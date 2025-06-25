import os
from argparse import ArgumentParser
from subprocess import check_call
from utils import setup_logging  # Import the logging function

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/2_phasing/"

parser = ArgumentParser(description="Phasing step for longshot.")
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--n_threads', '-t', help="number of threads",type = int, default=10, required=False)

args = parser.parse_args()

def run_longshot(out_dir, ref_file, logger, n_threads):
    use_cmd = f"python3 {code_path}longshot.py --ref_file {ref_file} --out_dir {out_dir} -t {n_threads}"
    logger.info(f"Executing: {use_cmd}")
    try:
        check_call(use_cmd, shell=True)
        logger.info("Longshot completed successfully.")
    except Exception as e:
        logger.error(f"Error running longshot: {e}")

def putback(out_dir, logger, n_threads):
    use_cmd = f"python3 {code_path}output_fas.py --out_dir {out_dir} -t {n_threads}"
    logger.info(f"Executing: {use_cmd}")
    try:
        check_call(use_cmd, shell=True)
        logger.info("Output FASTA files successfully.")
    except Exception as e:
        logger.error(f"Error generating FASTA output: {e}")

if __name__ == "__main__":
    out_dir = args.out_dir
    ref_file = args.ref_file
    n_threads = args.n_threads
    # Initialize logger
    logger = setup_logging("2_PHASING", out_dir)

    run_longshot(out_dir, ref_file, logger, n_threads)
    putback(out_dir, logger, n_threads)
