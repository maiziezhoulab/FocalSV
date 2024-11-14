#!/usr/bin/env python3
from subprocess import check_call, CalledProcessError
from argparse import ArgumentParser
import os
import sys
from utils import setup_logging  # Assuming setup_logging exists in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/"
__author__ = "Maizie&Jamie&Can@Vandy"

parser = ArgumentParser(description="Author: maiziezhoulab@gmail.com\n", usage='use "python3 %(prog)s --help" for more information')

# General inputs
parser.add_argument('--bam_file', '-bam', help="BAM file", required=True)
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)

# Defaulted inputs
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results, default = ./RegionBased_results", default="./RegionBased_results")
parser.add_argument('--data_type', '-d', help="HIFI/CLR/ONT", default="HIFI")

# Process information
parser.add_argument('--num_cpus', '-cpu', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--num_threads', '-thread', type=int, help="Number of threads, default = 8 (recommended)", default=8)

args = parser.parse_args()

def run_command(command, logger, step):
    logger.info(f"Executing {step}: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"{step} completed successfully.")
    except CalledProcessError as e:
        logger.error(f"Error during {step}: {e}")
        raise

def Evaluation(bam_file, chr_num, out_dir, reference, threads, cpus, datatype, logger):
    use_cmd = (f"python3 {code_path}5_evaluation_gtcorr.py "
               f"--bam_file {bam_file} "
               f"--chr_num {chr_num} "
               f"--reference {reference} "
               f"-o {out_dir} "
               f"--num_threads {threads} "
               f"--num_cpus {cpus} "
               f"--data_type {datatype}")
    
    run_command(use_cmd, logger, "Evaluation")

def main():
    if len(sys.argv) == 1:
        print("Displaying help for evaluation.py")
        check_call("python3 evaluation.py -h", shell=True)
    else:
        # Parse command line arguments
        bam_file = args.bam_file
        ref_file = args.ref_file
        chr_num = args.chr_num
        out_dir = args.out_dir
        data_type = args.data_type
        num_threads = args.num_threads
        num_cpus = args.num_cpus
        
        # Create output directory if it doesn't exist
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            print(f"Created output folder: {out_dir}")
        else:
            print(f"Using existing output folder: {out_dir}")
            
        # Convert data_type to integer
        data_type_map = {"HIFI": 0, "CLR": 1, "ONT": 2}
        data_type = data_type_map.get(data_type, 0)  # Default to HIFI if not recognized

        # Initialize logging
        logger = setup_logging("EVALUATION", out_dir)

        try:
            # Clean up and evaluation
            logger.info("Starting evaluation step...")
            Evaluation(bam_file, chr_num, out_dir, ref_file, num_threads, num_cpus, data_type, logger)

            logger.info("All steps completed successfully.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()
