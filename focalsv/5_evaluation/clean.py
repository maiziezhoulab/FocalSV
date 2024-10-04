import os
from argparse import ArgumentParser
from joblib import Parallel, delayed
from utils import setup_logging  # Assuming `setup_logging` is available in utils.py

parser = ArgumentParser(description="Clean data:")
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", default="./RegionBased_results")
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--data_type', '-d', type=int, help="Data type: CLR (1), HIFI (0), ONT (2)", default=0)

def clean(fd, datatype, logger):
    logger.info(f"Cleaning folder: {fd}")
    
    # Delete all non `.fa` files under 'PS'
    if datatype == 1:  # CLR lib
        delete = [os.path.join(fd, f) for f in os.listdir(fd) if f.startswith("PS") and not f.endswith(".fa")]
        for folder in delete:
            in_pb = [os.path.join(folder, f) for f in os.listdir(folder) if f != "assembly.fasta"]
            for f in in_pb:
                logger.info(f"Removing: {f}")
                os.system(f"rm -r {f}")

    elif datatype == 0:  # HIFI
        delete = [os.path.join(fd, f) for f in os.listdir(fd) if f.startswith("PS") and not f.endswith(".fa")]
        for file in delete:
            logger.info(f"Removing: {file}")
            os.system(f"rm {file}")

    elif datatype == 2:  # ONT
        delete = [os.path.join(fd, f) for f in os.listdir(fd) if f.startswith("PS") and not f.endswith(".fa")]
        for folder in delete:
            in_pb = [os.path.join(folder, f) for f in os.listdir(folder) if not f.endswith(".fasta")]
            for f in in_pb:
                logger.info(f"Removing: {f}")
                os.system(f"rm -r {f}")

if __name__ == "__main__":
    args = parser.parse_args()

    out_dir_region = os.path.join(args.out_dir, "regions")
    cpu = args.num_cpus
    datatype = args.data_type

    # Initialize logger
    logger = setup_logging("5_clean", args.out_dir)

    logger.info("***START Clean up***")
    fds = [os.path.join(out_dir_region, fd) for fd in os.listdir(out_dir_region) if fd.startswith("Region")]
    Parallel(n_jobs=cpu)(delayed(clean)(fds[i], datatype, logger) for i in range(len(fds)))
    logger.info("***FINISHED Clean up***")
