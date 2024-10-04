import os
import glob
from argparse import ArgumentParser
from utils import setup_logging  # Assuming setup_logging exists in utils.py

parser = ArgumentParser(description="Assemble sequences:")
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--data_type', '-d', type=int, help="hifi = 0 or clr = 1 or ont = 2 data", default=0)

def combine_fas_hifi(out_dir, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    for fd in fds:
        fas1 = [os.path.join(fd, f) for f in os.listdir(fd) if f.endswith("hp1.asm.p_ctg.gfa.fa")]
        fas2 = [os.path.join(fd, f) for f in os.listdir(fd) if f.endswith("hp2.asm.p_ctg.gfa.fa")]

        new_hp1 = os.path.join(fd, "HP1.fa")
        new_hp2 = os.path.join(fd, "HP2.fa")

        try:
            with open(new_hp1, 'w') as outfile1:
                for fa1 in fas1:
                    with open(fa1, 'r') as infile1:
                        outfile1.writelines(infile1.read())
            logger.info(f"*** Finished combining HP1 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP1 for {fd}: {e}")

        try:
            with open(new_hp2, 'w') as outfile2:
                for fa2 in fas2:
                    with open(fa2, 'r') as infile2:
                        outfile2.writelines(infile2.read())
            logger.info(f"*** Finished combining HP2 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP2 for {fd}: {e}")

def combine_fas_clr(out_dir, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    for fd in fds:
        fas1 = [os.path.join(fd, f, "assembly.fasta") for f in os.listdir(fd) if f.endswith("hp1_flye")]
        fas2 = [os.path.join(fd, f, "assembly.fasta") for f in os.listdir(fd) if f.endswith("hp2_flye")]

        new_hp1 = os.path.join(fd, "HP1.fa")
        new_hp2 = os.path.join(fd, "HP2.fa")

        try:
            with open(new_hp1, 'w') as outfile1:
                for fa1 in fas1:
                    try:
                        with open(fa1, 'r') as infile1:
                            outfile1.writelines(infile1.read())
                    except FileNotFoundError:
                        logger.warning(f"{fa1} does not exist!")
            logger.info(f"*** Finished combining HP1 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP1 for {fd}: {e}")

        try:
            with open(new_hp2, 'w') as outfile2:
                for fa2 in fas2:
                    try:
                        with open(fa2, 'r') as infile2:
                            outfile2.writelines(infile2.read())
                    except FileNotFoundError:
                        logger.warning(f"{fa2} does not exist!")
            logger.info(f"*** Finished combining HP2 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP2 for {fd}: {e}")

def combine_fas_ont(out_dir, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    for fd in fds:
        fas1 = [os.path.join(fd, f, "/*ssembly.fasta") for f in os.listdir(fd) if f.endswith("hp1_flye")]
        fas2 = [os.path.join(fd, f, "/*ssembly.fasta") for f in os.listdir(fd) if f.endswith("hp2_flye")]

        new_hp1 = os.path.join(fd, "HP1.fa")
        new_hp2 = os.path.join(fd, "HP2.fa")

        try:
            with open(new_hp1, 'w') as outfile1:
                for fa1 in fas1:
                    try:
                        fa1 = glob.glob(fa1)[0]
                        with open(fa1, 'r') as infile1:
                            outfile1.writelines(infile1.read())
                    except IndexError:
                        logger.warning(f"{fa1} does not exist!")
            logger.info(f"*** Finished combining HP1 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP1 for {fd}: {e}")

        try:
            with open(new_hp2, 'w') as outfile2:
                for fa2 in fas2:
                    try:
                        fa2 = glob.glob(fa2)[0]
                        with open(fa2, 'r') as infile2:
                            outfile2.writelines(infile2.read())
                    except IndexError:
                        logger.warning(f"{fa2} does not exist!")
            logger.info(f"*** Finished combining HP2 for {fd} ***")
        except Exception as e:
            logger.error(f"Error combining HP2 for {fd}: {e}")

if __name__ == "__main__":
    args = parser.parse_args()
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    data_type = args.data_type

    # Initialize logger
    logger = setup_logging("3_combine_fas", out_dir_general)

    if data_type == 0:
        logger.info("Combining FASTA files for HIFI")
        combine_fas_hifi(out_dir, logger)
    elif data_type == 1:
        logger.info("Combining FASTA files for CLR")
        combine_fas_clr(out_dir, logger)
    elif data_type == 2:
        logger.info("Combining FASTA files for ONT")
        combine_fas_ont(out_dir, logger)
