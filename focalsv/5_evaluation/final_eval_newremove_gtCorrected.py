import os
from argparse import ArgumentParser
from subprocess import Popen, CalledProcessError
from utils import setup_logging  # Assuming `setup_logging` is available in utils.py

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path

parser = ArgumentParser(description="Combine VCFs and evaluate:")
parser.add_argument('--results_dir', '-r', help="Results directory")
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number")

def run_command(cmd, logger, step):
    logger.info(f"Executing {step}: {cmd}")
    try:
        Popen(cmd, shell=True).wait()
        logger.info(f"{step} completed successfully.")
    except CalledProcessError as e:
        logger.error(f"Error during {step}: {e}")
        raise
    
def truvari_eval(chr_num, results_dir, logger):
    eval_dir = os.path.join(results_dir, 'eval/')
    cmd = f"{code_path}/truvari_eval.sh {chr_num} {results_dir} {eval_dir} FocalSV_variant_no_redundancy 500 0.5 0.5 30 0.01"
    run_command(cmd, logger, "Truvari Evaluation")

if __name__ == "__main__":
    args = parser.parse_args()

    results_dir = os.path.join(args.results_dir, "results_gtCorrected")
    chr_num = args.chr_num

    # Initialize logger
    logger = setup_logging("5_final_eval", args.results_dir)

    truvari_eval(chr_num, results_dir, logger)
