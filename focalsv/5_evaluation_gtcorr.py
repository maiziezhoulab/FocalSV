import os
from argparse import ArgumentParser
from subprocess import check_call, CalledProcessError
from utils import setup_logging  # Assuming setup_logging exists in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/5_evaluation/"

parser = ArgumentParser(description="Evaluate final results:")
# General inputs
parser.add_argument('--bam_file', '-bam', help="BAM file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--reference', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--num_threads', '-t_chr', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--data_type', '-d', type=int, help="HIFI = 0 or CLR = 1 data", default=0)

args = parser.parse_args()

def run_command(command, logger, step):
    logger.info(f"Executing {step}: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"{step} completed successfully.")
    except CalledProcessError as e:
        logger.error(f"Error during {step}: {e}")
        raise

def clean(out_dir, cpu, datatype, logger):
    use_cmd = (f"python3.7 {code_path}clean.py "
               f"-o {out_dir} "
               f"-t {cpu} "
               f"-d {datatype}")
    
    run_command(use_cmd, logger, "Clean")

def fullchr_eval_gtCorrected(chr_num, out_dir, logger):
    target_sv = os.path.join(out_dir, "target_sv.vcf")
    input_regions = os.path.join(out_dir, "regions")
    results_dir = out_dir

    use_cmd = (f"python3.7 {code_path}final_eval_newremove_gtCorrected.py "
               f"-r {results_dir} "
               f"-chr {chr_num}")
    
    run_command(use_cmd, logger, "Full Chromosome Evaluation (GT Corrected)")

if __name__ == "__main__":
    bam_file = args.bam_file
    chr_num = args.chr_num
    reference = args.reference
    out_dir = args.out_dir
    threads = args.num_threads
    cpu = args.num_cpus
    datatype = args.data_type

    # Initialize logger
    logger = setup_logging("5_EVALUTAION", out_dir)

    logger.info(f"Starting evaluation for chromosome {chr_num} with data type {datatype}")
    
    # Step 1: Full chromosome evaluation with GT correction
    fullchr_eval_gtCorrected(chr_num, out_dir, logger)
    
    # Step 2: Clean up the results
    clean(out_dir, cpu, datatype, logger)

    logger.info("Evaluation and cleanup completed successfully.")
