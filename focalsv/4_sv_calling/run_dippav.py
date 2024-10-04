import os
import sys
from argparse import ArgumentParser
from joblib import Parallel, delayed
from utils import setup_logging  # Assuming setup_logging exists in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, script_path + '/Dippav')
from DipPAV_variant_call import dippav_variant_call

parser = ArgumentParser(description="Run DipPAV variant calling:")
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--data_type', '-d', type=int, help="HIFI = 0 or CLR = 1 or ONT = 2 data", default=0)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)

def dippav_run(chr_num, out_dir, cpu, data_type, logger):
    # Get all folders starting with "Region"
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Processing regions: {fds}")

    bam_fd, hp1_fd, hp2_fd, output_fd = [], [], [], []
    for fd in fds:
        bam_file = os.path.join(fd, "region.bam")
        hp1_contig = os.path.join(fd, "HP1.fa")
        hp2_contig = os.path.join(fd, "HP2.fa")
        outputdir = os.path.join(fd, "results")

        # Validate the existence of BAM and contig files
        if not all([os.path.exists(bam_file), os.path.exists(hp1_contig), os.path.exists(hp2_contig)]):
            logger.error(f"Missing required files in {fd}: region.bam, HP1.fa, or HP2.fa")
            continue

        bam_fd.append(bam_file)
        hp1_fd.append(hp1_contig)
        hp2_fd.append(hp2_contig)
        output_fd.append(outputdir)

        # Logging the command
        cmd = (f"python3 DipPAV_variant_call.py -rbam {bam_file} -hp1 {hp1_contig} -hp2 {hp2_contig} "
               f"-o {outputdir} -chr {chr_num} -dtype {data_type}")
        logger.info(f"Executing: {cmd}")

    n = len(bam_fd)
    if n == 0:
        logger.error("No valid regions found for processing.")
        return

    # Run the variant calling in parallel
    try:
        Parallel(n_jobs=cpu)(
            delayed(dippav_variant_call)(data_type, bam_fd[i], hp1_fd[i], hp2_fd[i], output_fd[i], chr_num) for i in range(n)
        )
        logger.info("DipPAV variant calling completed successfully.")
    except Exception as e:
        logger.error(f"Error during DipPAV variant calling: {e}")

if __name__ == "__main__":
    args = parser.parse_args()
    chr_num = args.chr_num
    cpu = args.num_cpus
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    data_type_map = {0: "CCS", 1: "CLR", 2: "ONT"}
    data_type = data_type_map.get(args.data_type, "UNKNOWN")

    # Initialize logger
    logger = setup_logging("4_run_dippav", out_dir_general)

    logger.info(f"Starting DipPAV variant calling for chromosome {chr_num} with data type {data_type}")
    dippav_run(chr_num, out_dir, cpu, data_type, logger)
