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

# For single region
parser.add_argument('--region_start', '-S', type=int, help="Target region starting index", required=False)
parser.add_argument('--region_end', '-E', type=int, help="Target region ending index", required=False)

# For multiple regions
parser.add_argument('--target_bed', '-target_bed', help="BED file with multiple target regions", required=False)

# Defaulted inputs
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results, default = ./RegionBased_results", default="./RegionBased_results")
parser.add_argument('--data_type', '-d', help="HIFI/CLR/ONT", default="HIFI")

# Process information
parser.add_argument('--num_cpus', '-cpu', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--num_threads', '-thread', type=int, help="Number of threads, default = 8 (recommended)", default=8)

# Whole Chromosome Evaluation
parser.add_argument('--eval', '-e', action='store_true', help="Whole Chromsome evaluation")
parser.add_argument('--bed_file', '-bed', help="Bed file, used to filter SVs (-whole)")
parser.add_argument('--svlen_threshold', '-sv_thresh', type=int, help="Threshold number for initially filtering SVs (-whole)", default=50)
parser.add_argument('--vcf_file', '-v', help="VCF file, called by FreeBayes (-whole)", required=False)
parser.add_argument('--flanking', '-flank', type=int, help="Length of flanking region", default=50000)

args = parser.parse_args()

def run_command(command, logger, step):
    logger.info(f"Executing {step}: {command}")
    try:
        check_call(command, shell=True)
        logger.info(f"{step} completed successfully.")
    except CalledProcessError as e:
        logger.error(f"Error during {step}: {e}")
        raise

def Filter_SVs(vcf_file, chr_num, out_dir, bed_file, SV_THRESH, FLANK, logger):
    use_cmd = (f"python3 {code_path}0_filter_sv.py "
               f"--vcf_file {vcf_file} "
               f"--bed_file {bed_file} "
               f"--svlen_threshold {SV_THRESH} "
               f"--flanking {FLANK} "
               f"--chr_num {chr_num} "
               f"--out_dir {out_dir}")
    run_command(use_cmd, logger, "Filter SVs")

def Crop_Bam(bam_file, chr_num, target_bed, region_start, region_end, out_dir, logger):
    if target_bed:
        use_cmd = (f"python3 {code_path}1_crop_bam.py "
                   f"--bam_file {bam_file} "
                   f"--chr_num {chr_num} "
                   f"--target_bed {target_bed} "
                   f"--out_dir {out_dir}")
    else:
        use_cmd = (f"python3 {code_path}1_crop_bam.py "
                   f"--bam_file {bam_file} "
                   f"--chr_num {chr_num} "
                   f"--region_start {region_start} "
                   f"--region_end {region_end} "
                   f"--out_dir {out_dir}")
    
    run_command(use_cmd, logger, "BAM Cropping")

def Phase_Bam(ref_file, out_dir, logger):
    use_cmd = (f"python3 {code_path}2_phasing.py "
               f"--ref_file {ref_file} "
               f"--out_dir {out_dir}")
    
    run_command(use_cmd, logger, "Phasing")

def Assembly(bam_file, chr_num, out_dir, reference, threads, cpus, log_dir, datatype, logger):
    use_cmd = (f"python3 {code_path}3_assembly.py "
               f"--bam_file {bam_file} "
               f"--chr_num {chr_num} "
               f"--ref_file {reference} "
               f"-o {out_dir} "
               f"--num_threads {threads} "
               f"--num_cpus {cpus} "
               f"--data_type {datatype}")
    
    run_command(use_cmd, logger, "Assembly")

def SVCalling(bam_file, chr_num, out_dir, reference, threads, cpus, log_dir, datatype, logger):
    use_cmd = (f"python3 {code_path}4_sv_calling.py "
               f"--bam_file {bam_file} "
               f"--chr_num {chr_num} "
               f"--reference {reference} "
               f"-o {out_dir} "
               f"--num_threads {threads} "
               f"--num_cpus {cpus} "
               f"--log_dir {log_dir} "
               f"--data_type {datatype}")
    
    run_command(use_cmd, logger, "SV Calling")

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
        print("Displaying help for main.py")
        check_call("python3 main.py -h", shell=True)
    else:
        # Parse command line arguments
        bam_file = args.bam_file
        ref_file = args.ref_file
        chr_num = args.chr_num
        region_start = args.region_start
        region_end = args.region_end
        out_dir = args.out_dir
        data_type = args.data_type
        num_threads = args.num_threads
        num_cpus = args.num_cpus
        evaluation = args.eval
        target_bed = args.target_bed
        
        if evaluation:
            if target_bed or region_start is not None or region_end is not None:
                parser.error("For whole chromosome analysis (--eval), do not provide --target_bed or --region_start/--region_end.")
        elif target_bed and (region_start is not None or region_end is not None):
            parser.error("Specify either --target_bed or --region_start/--region_end, but not both.")
        elif not target_bed and (region_start is None or region_end is None):
            parser.error("You must specify both --region_start and --region_end if --target_bed is not provided, unless using --eval for whole chromosome analysis.")

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
        logger = setup_logging("MAIN", out_dir)

        try:
            # Step 0: Filter SVs (ran with -whole set)
            if evaluation:
                logger.info("Starting Filter SV step...")
                bed_file = args.bed_file
                SV_THRESH = args.svlen_threshold
                FLANK = args.flanking
                vcf_file = args.vcf_file
                
                if not bed_file or not vcf_file:
                    logger.error("Both --bed_file and --vcf_file are required for whole chromosome analysis.")
                    sys.exit(1)
                
                # Filter_SVs(vcf_file, chr_num, out_dir, bed_file, SV_THRESH, FLANK, logger)
            
            # Step 1: Crop target region
            logger.info("Starting BAM cropping step...")
            if eval:
                # When evaluation mode is on, use the generated BED file for cropping
                output_bed = os.path.join(out_dir, "target_regions.bed")
                # Crop_Bam(bam_file, chr_num, output_bed, None, None, out_dir, logger)
            elif target_bed:
                # If the user provided a BED file, use it for cropping
                Crop_Bam(bam_file, chr_num, target_bed, None, None, out_dir, logger)
            else:
                # Otherwise, use the start and end region for cropping
                Crop_Bam(bam_file, chr_num, None, region_start, region_end, out_dir, logger)


            # Step 2: Phase Bams
            logger.info("Starting BAM phasing step...")
            # Phase_Bam(ref_file, out_dir, logger)

            # Step 3: Assemble
            logger.info("Starting assembly step...")
            # Assembly(bam_file, chr_num, out_dir, ref_file, num_threads, num_cpus, out_dir, data_type, logger)

            # Step 4: SV Calling
            logger.info("Starting SV calling step...")
            # SVCalling(bam_file, chr_num, out_dir, ref_file, num_threads, num_cpus, out_dir, data_type, logger)

            # Step 5: Clean up and evaluation
            if evaluation:
                logger.info("Starting evaluation step...")
                Evaluation(bam_file, chr_num, out_dir, ref_file, num_threads, num_cpus, data_type, logger)

            logger.info("All steps completed successfully.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()
