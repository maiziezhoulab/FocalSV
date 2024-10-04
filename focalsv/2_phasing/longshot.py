import os
from argparse import ArgumentParser
from subprocess import check_call
from utils import setup_logging  # Import the logging function

parser = ArgumentParser(description="Run longshot for phasing BAM.")
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)

args = parser.parse_args()

if __name__ == "__main__":
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    ref_file = args.ref_file

    # Initialize logger
    logger = setup_logging("2_longshot", out_dir_general)

    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]

    for fd in fds:
        bam = os.path.join(fd, "region.bam")
        phased_bam = os.path.join(fd, "region_phased.bam")
        phased_vcf = os.path.join(fd, "region_phased.vcf")

        cmd = f"longshot -O {phased_bam} -b {bam} -f {ref_file} -o {phased_vcf} -F"
        logger.info(f"Running longshot on {bam}")

        try:
            check_call(cmd, shell=True)
            logger.info(f"Longshot successfully run on {bam}")
            with open(os.path.join(out_dir, "2_phasing.txt"), "a") as log_file:
                log_file.write(f"{cmd}\n")
        except Exception as e:
            logger.error(f"Error running longshot on {bam}: {e}")
