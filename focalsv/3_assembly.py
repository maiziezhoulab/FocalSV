import os
from argparse import ArgumentParser
from subprocess import check_call, CalledProcessError
from utils import setup_logging  # Assuming setup_logging is in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/3_assembly/"

parser = ArgumentParser(description="Assemble sequences and call SVs:")
# General inputs
parser.add_argument('--bam_file', '-bam', help="BAM file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', required=True, help="Output directory")
parser.add_argument('--num_threads', '-t_chr', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--data_type', '-d', type=int, help="HIFI = 0 or CLR = 1 data", default=0)

def run_command(command, logger):
    logger.info(f"Executing: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"Command finished successfully: {command}")
    except CalledProcessError as e:
        logger.error(f"Error occurred during execution: {e}")
        raise

def assembly(out_dir, cpu, threads, data_type, logger):
    # Run assembly step
    run_command(f"python3 {code_path}run_assembly.py -o {out_dir} -t {cpu} -tc {threads} -d {data_type}", logger)
    
    # Run post-assembly step
    run_command(f"python3 {code_path}post_assembly.py -o {out_dir} -t {cpu} -tc {threads} -d {data_type}", logger)
    
    # Combine FASTA files
    run_command(f"python3 {code_path}combine_fas.py -o {out_dir} -d {data_type}", logger)


if __name__ == "__main__":
    args = parser.parse_args()
    bam_file = args.bam_file
    chr_num = args.chr_num
    ref_file = args.ref_file
    out_dir = args.out_dir
    threads = args.num_threads
    cpus = args.num_cpus
    data_type = args.data_type
    
    # Initialize logger
    logger = setup_logging("3_ASSEMBLY", out_dir)

    # Perform assembly steps
    logger.info("Starting assembly process")
    assembly(out_dir, cpus, threads, data_type, logger)

    
