import os
from argparse import ArgumentParser
from joblib import Parallel, delayed
from subprocess import check_call
from utils import setup_logging  # Import the logging function

parser = ArgumentParser(description="Run longshot for phasing BAM.")
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--n_threads', '-t', help="number of threads",type = int, default=10, required=False)

args = parser.parse_args()


def run_longshot(fd, ref_file, logger):
    bam = os.path.join(fd, "region.bam")
    phased_bam = os.path.join(fd, "region_phased.bam")
    phased_vcf = os.path.join(fd, "region_phased.vcf")

    cmd = f"longshot -O {phased_bam} -b {bam} -f {ref_file} -o {phased_vcf} -F"
    logger.info(f"Running longshot on {bam}")

    try:
        check_call(cmd, shell=True)
        logger.info(f"Longshot successfully run on {bam}")
    except Exception as e:
        logger.error(f"Error running longshot on {bam}: {e}")

def run_longshot_parallel(out_dir, ref_file,logger, n_jobs=4):

    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    
    Parallel(n_jobs=n_jobs)(delayed(run_longshot)(fd, ref_file, logger) for fd in fds)



if __name__ == "__main__":
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    ref_file = args.ref_file
    n_threads = args.n_threads
    logger = setup_logging("2_longshot", out_dir_general)
    run_longshot_parallel(out_dir, ref_file, logger, n_jobs=n_threads)
