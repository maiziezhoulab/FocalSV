import os
from joblib import Parallel, delayed
from argparse import ArgumentParser
from utils import setup_logging  # Import the logging setup function
import logging

parser = ArgumentParser(description="Assemble sequences:")
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--num_threads', '-tc', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--data_type', '-d', type=int, help="CLR = 1 or ONT = 2 data", default=0)

def run_clr(out_dir, threads, input_file, logger):
    cmd = f"flye --pacbio-raw {input_file} -o {out_dir} -t {threads}"
    logger.info(f"Executing: {cmd}")
    try:
        os.system(cmd)
        logger.info(f"FLYE assembly completed for {input_file}")
    except Exception as e:
        logger.error(f"Error during FLYE assembly for {input_file}: {e}")

def rerun_clr(out_dir, threads, cpu, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Regions found: {fds}")

    unvalid = []
    for fd in fds:
        fas = [os.path.join(fd, f, "assembly.fasta") for f in os.listdir(fd) if f.endswith("flye")]
        for fa in fas:
            if not os.path.exists(fa):
                unvalid.append(fa)
    
    inputs = []
    outputs = []
    for pbhp in unvalid:
        inputs.append(pbhp.replace("_flye/assembly.fasta", ".fa"))
        outputs.append(pbhp.replace("/assembly.fasta", ""))
    
    n = len(inputs)
    logger.info(f"Re-running CLR assembly for {n} invalid regions")
    Parallel(n_jobs=cpu)(delayed(run_clr)(outputs[i], threads, inputs[i], logger) for i in range(n))

def run_ont(out_dir, input_file, logger):
    try:
        os.system(f"rm -r {out_dir}")
        os.system(f"mkdir {out_dir}")
        cmd = f"../software/Shasta/shasta-Linux-0.10.0 --input {input_file} --assemblyDirectory {out_dir} --config Nanopore-UL-Dec2019"
        logger.info(f"Executing: {cmd}")
        os.system(cmd)
        cmd_cleanup = f"../software/Shasta/shasta-Linux-0.10.0 --command cleanupBinaryData --assemblyDirectory {out_dir}"
        os.system(cmd_cleanup)
        logger.info(f"ONT assembly completed for {input_file}")
    except Exception as e:
        logger.error(f"Error during ONT assembly for {input_file}: {e}")

def rerun_ont(out_dir, threads, cpu, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Regions found: {fds}")

    unvalid = []
    for fd in fds:
        fas = [os.path.join(fd, f, "assembly.fasta") for f in os.listdir(fd) if f.endswith("flye")]
        for fa in fas:
            if not os.path.exists(fa):
                unvalid.append(fa)

    inputs = []
    outputs = []
    for pbhp in unvalid:
        inputs.append(pbhp.replace("_flye/assembly.fasta", ".fa"))
        outputs.append(pbhp.replace("/assembly.fasta", ""))
    
    n = len(inputs)
    logger.info(f"Re-running ONT assembly for {n} invalid regions")
    Parallel(n_jobs=cpu)(delayed(run_ont)(outputs[i], inputs[i], logger) for i in range(n))

# For HIFI: index data
def fa_index(out_dir, logger):
    logger.info(f"HIFI: fa index {out_dir}")
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    file_count = 0

    for fd in fds:
        fas = [os.path.join(fd, f) for f in os.listdir(fd) if f.endswith("asm.p_ctg.gfa")]
        file_count += 1

        # call on each fa
        for fa in fas:
            input_file = fa
            output_file = f"{fa}.fa"
            cmd = f"awk '/^S/{{print \">{output_file}\\n\"$3}}' {input_file} | fold > {output_file}"
            logger.info(f"Executing: {cmd}")
            os.system(cmd)
            logger.info(f"FA index completed for {output_file}")

if __name__ == "__main__":
    args = parser.parse_args()
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    threads = args.num_threads
    cpu = args.num_cpus
    data_type = args.data_type

    # Initialize logger
    logger = setup_logging("3_post_assembly", out_dir_general)

    if data_type == 0:
        logger.info("Starting HIFI FA index")
        fa_index(out_dir, logger)
    elif data_type == 1:
        for i in range(3):
            logger.info(f"Starting CLR assembly, iteration {i+1}")
            rerun_clr(out_dir, threads, cpu, logger)
    elif data_type == 2:
        logger.info("Starting ONT assembly")
        rerun_ont(out_dir, threads, cpu, logger)
