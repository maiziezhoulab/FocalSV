import os
from argparse import ArgumentParser
from subprocess import check_call, CalledProcessError
from utils import setup_logging  # Assuming setup_logging exists in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/4_sv_calling/"

parser = ArgumentParser(description="Assemble sequences:")
# General inputs
parser.add_argument('--bam_file', '-bam', help="BAM file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--reference', '-r', help="Reference fasta file", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", default="./RegionBased_results")
parser.add_argument('--num_threads', '-t_chr', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--log_dir', '-log', help="Position of log directory", required=True)
parser.add_argument('--data_type', '-d', type=int, help="HIFI = 0 or CLR = 1 data", default=0)

def run_command(command, logger):
    logger.info(f"Executing: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"Command finished successfully: {command}")
    except CalledProcessError as e:
        logger.error(f"Error occurred during execution: {e}")
        raise

def sv_calling(chr_num, out_dir, cpus, data_type, logger):
    # Run SV calling
    run_command(f"python3 {code_path}run_dippav.py -chr {chr_num} -o {out_dir} -t {cpus} -d {data_type}", logger)

def remove_redundancy(out_dir, logger):
    target_sv = os.path.join(out_dir, "target_sv.vcf")
    input_regions = os.path.join(out_dir, "regions")
    run_command(f"python3 {code_path}remove_redundancy_region_based.py -t {target_sv} -rg {input_regions} -o {out_dir} -vcf dippav_variant_no_redundancy.vcf", logger)

if __name__ == "__main__":
    args = parser.parse_args()
    bam_file = args.bam_file
    chr_num = args.chr_num
    reference = args.reference
    out_dir = args.out_dir
    threads = args.num_threads
    cpus = args.num_cpus
    log_dir = args.log_dir
    datatype = args.data_type

    # Initialize logger
    logger = setup_logging("4_SV_CALLING", out_dir)

    # Perform SV calling
    logger.info("Starting SV calling process")
    sv_calling(chr_num, out_dir, cpus, datatype, logger)
    
    # Remove redundancy
    logger.info("Removing redundancy")
    remove_redundancy(out_dir, logger)
