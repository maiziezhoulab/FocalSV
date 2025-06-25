import os
from joblib import Parallel, delayed
from argparse import ArgumentParser
from utils import setup_logging  # Import the logging setup function

parser = ArgumentParser(description="Assemble sequences:")
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", required=True)
parser.add_argument('--num_cpus', '-t', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--num_threads', '-tc', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--data_type', '-d', type=int, help="HIFI = 0, CLR = 1, ONT = 2 data", default=0)




def run_hifiasm(output, threads, input_file, logger):
    code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
    if 'unphased' in input_file:
        v = "0.16.1"
    else:
        v = "0.14"
    cmd = f"{code_dir}/../../software/hifiasm-{v}/hifiasm -o {output} -t {threads} {input_file}"
    print(cmd)
    logger.info(f"Executing: {cmd}")
    try:
        os.system(cmd)
        logger.info(f"HIFI assembly completed for {input_file}")
    except Exception as e:
        logger.error(f"Error during HIFI assembly for {input_file}: {e}")

def assembly_hifi(out_dir, threads, cpu, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Regions found: {fds}")
    
    inputs, outputs = [], []
    for fd in fds:
        fas = [os.path.join(fd, f) for f in os.listdir(fd) if f.endswith(".fa")]
        for fa in fas:
            inputs.append(fa)
            outputs.append(fa[:-3] + ".asm")
    
    n = len(inputs)
    logger.info(f"Starting HIFI assembly for {n} regions.")
    Parallel(n_jobs=cpu)(delayed(run_hifiasm)(outputs[i], threads, inputs[i], logger) for i in range(n))
    logger.info("HIFI assembly finished.")

def run_flye(out_dir, threads, input_file, logger):
    cmd = f"flye --pacbio-raw {input_file} -o {out_dir} -t {threads}"
    logger.info(f"Executing: {cmd}")
    try:
        os.system(cmd)
        logger.info(f"CLR assembly completed for {input_file}")
    except Exception as e:
        logger.error(f"Error during CLR assembly for {input_file}: {e}")

def assembly_clr(out_dir, threads, cpu, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Regions found: {fds}")
    
    inputs, outputs = [], []
    for fd in fds:
        fas = [os.path.join(fd, f) for f in os.listdir(fd) if (f.startswith("PS") and f.endswith(".fa"))]
        for fa in fas:
            outdir_flye = fa[:-3] + "_flye"
            os.makedirs(outdir_flye, exist_ok=True)
            if not os.path.exists(os.path.join(outdir_flye, "assembly.fasta")):
                inputs.append(fa)
                outputs.append(outdir_flye)
    
    n = len(inputs)
    logger.info(f"Starting CLR assembly for {n} regions.")
    Parallel(n_jobs=cpu)(delayed(run_flye)(outputs[i], threads, inputs[i], logger) for i in range(n))
    logger.info("CLR assembly finished.")

def run_flye_ont(out_dir, threads, input_file, logger):
    cmd = f"flye --nano-raw {input_file} -o {out_dir} -t {threads}"
    logger.info(f"Executing: {cmd}")
    try:
        os.system(cmd)
        logger.info(f"ONT assembly completed for {input_file}")
    except Exception as e:
        logger.error(f"Error during ONT assembly for {input_file}: {e}")

def assembly_ont(out_dir, threads, cpu, logger):
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    logger.info(f"Regions found: {fds}")
    
    inputs, outputs = [], []
    for fd in fds:
        fas = [os.path.join(fd, f) for f in os.listdir(fd) if (f.startswith("PS") and f.endswith(".fa"))]
        for fa in fas:
            outdir_flye = fa[:-3] + "_flye"
            os.makedirs(outdir_flye, exist_ok=True)
            if not os.path.exists(os.path.join(outdir_flye, "assembly.fasta")):
                inputs.append(fa)
                outputs.append(outdir_flye)
    
    n = len(inputs)
    logger.info(f"Starting ONT assembly for {n} regions.")
    Parallel(n_jobs=cpu)(delayed(run_flye_ont)(outputs[i], threads, inputs[i], logger) for i in range(n))
    logger.info("ONT assembly finished.")

if __name__ == "__main__":
    args = parser.parse_args()
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")
    threads = args.num_threads
    cpu = args.num_cpus
    data_type = args.data_type

    # Initialize logger
    logger = setup_logging("3_assembly", out_dir_general)

    # Perform assembly based on data type
    if data_type == 0:
        logger.info("Starting HIFI assembly")
        assembly_hifi(out_dir, threads, cpu, logger)
    elif data_type == 1:
        logger.info("Starting CLR assembly")
        assembly_clr(out_dir, threads, cpu, logger)
    elif data_type == 2:
        logger.info("Starting ONT assembly")
        assembly_ont(out_dir, threads, cpu, logger)
