import os
from argparse import ArgumentParser
from subprocess import check_call, CalledProcessError
from utils import setup_logging

# Paths
script_path = os.path.dirname(os.path.abspath(__file__))
code_path = os.path.join(script_path, "evaluation")

# Argument Parser
parser = ArgumentParser(description="Evaluate final results")
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--bench_vcf', '-vcf', help="Benchmark VCF file", required=True)
parser.add_argument('--reference', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', help="Output directory with results from the previous step", required=True)
parser.add_argument('--num_threads', '-t_chr', type=int, help="Number of threads (default = 8)", default=8)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs (default = 10)", default=10)
parser.add_argument('--data_type', '-d', type=int, help="Data type: HIFI (0) or CLR (1), default = 0", default=0)

args = parser.parse_args()

def run_command(command, logger, step):
    logger.info(f"Executing {step}: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"{step} completed successfully.")
    except CalledProcessError as e:
        logger.error(f"Error during {step}: {e}")
        raise

def clean(out_dir, num_cpus, data_type, logger):
    command = (
        f"python3.7 {os.path.join(code_path, 'clean.py')} "
        f"-o {out_dir} "
        f"-t {num_cpus} "
        f"-d {data_type}"
    )
    run_command(command, logger, "Cleanup")

def truvari_eval(chr_num, results_dir, vcf, logger):
    eval_dir = os.path.join(results_dir, 'eval')
    input_dir = os.path.join(results_dir, 'results')
    command = (
        f"{os.path.join(code_path, 'truvari_eval.sh')} "
        f"{chr_num} {input_dir} {eval_dir} FocalSV_variant_no_redundancy "
        f"500 0.5 0.5 30 0.01"
    )
    run_command(command, logger, "Truvari Evaluation")

if __name__ == "__main__":
    chr_num = args.chr_num
    bench_vcf = args.bench_vcf
    reference = args.reference
    out_dir = args.out_dir
    num_threads = args.num_threads
    num_cpus = args.num_cpus
    data_type = args.data_type

    logger = setup_logging("EVALUATION", out_dir)
    logger.info(f"Starting evaluation for chromosome {chr_num} with data type {data_type}")

    # Step 1: Truvari evaluation
    truvari_eval(chr_num, out_dir, bench_vcf, logger)

    # Step 2: Clean up the results
    clean(out_dir, num_cpus, data_type, logger)

    logger.info("Evaluation and cleanup completed successfully.")
