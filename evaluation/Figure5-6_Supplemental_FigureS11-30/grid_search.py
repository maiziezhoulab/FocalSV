import os
import subprocess
from itertools import product
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import glob
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


def run_truvari(prefix, vcf, ref, bed, bench_del, bench_ins, output_root, p=None, r=None, O=None):
    # Construct output directory

    if p is not None and r is not None:
        subdir = f"p_{p}_r_{r}"
        truvari_cmds = [
            f"truvari bench -b {bench_del} -c {output_root}/{prefix}_DEL_noXY_sorted.vcf.gz -f {ref} -o {output_root}/{subdir}/DEL_50_ --includebed {bed} -p {p} -r {r} --passonly --sizemin 50",
            f"truvari bench -b {bench_ins} -c {output_root}/{prefix}_INS_noXY_sorted.vcf.gz -f {ref} -o {output_root}/{subdir}/INS_50_ --includebed {bed} -p {p} -r {r} --passonly --sizemin 50"
        ]
    elif p is not None and O is not None:
        subdir = f"p_{p}_O_{O}"
        truvari_cmds = [
            f"truvari bench -b {bench_del} -c {output_root}/{prefix}_DEL_noXY_sorted.vcf.gz -f {ref} -o {output_root}/{subdir}/DEL_50_ --includebed {bed} -p {p} -O {O} --passonly --sizemin 50",
            f"truvari bench -b {bench_ins} -c {output_root}/{prefix}_INS_noXY_sorted.vcf.gz -f {ref} -o {output_root}/{subdir}/INS_50_ --includebed {bed} -p {p} -O {O} --passonly --sizemin 50"
        ]
    else:
        raise ValueError("Unsupported parameter combination")

    os.makedirs(os.path.join(output_root, subdir, ), exist_ok=True)
    # os.makedirs(os.path.join(output_root, subdir, 'INS_50_'), exist_ok=True)

    for cmd in truvari_cmds:
        subprocess.run(cmd, shell=True, check=True)
    
    # clean
    for dtype in ['DEL','INS']:
        cmd = f"find {output_root}/{subdir}/{dtype}_50_ -type f ! -name 'summary.json' ! -name 'log.txt' -delete"
        subprocess.run(cmd, shell=True, check=True)

def process_one_vcf(vcf, prefix, work_dir, n_thread, ref, bench_dir):

    # prefix = "ontLib6"
    # vcf = "/path/to/your_input.vcf"
    # ref = "/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa"
    bed = f"{bench_dir}/HG002_SVs_Tier1_v0.6_chr_noXY.bed"
    bench_del = f"{bench_dir}/HG002_SVs_Tier1_v0.6_chr_del_noXY.vcf.gz"
    bench_ins = f"{bench_dir}/HG002_SVs_Tier1_v0.6_chr_ins_noXY.vcf.gz"
    output_root = os.path.abspath(f"{work_dir}/Truvari_{prefix}_hm")

    if os.path.exists(output_root):
        os.system("rm -r " + output_root)
    os.makedirs(output_root, exist_ok=True)
    os.chdir(output_root)

    # Run the VCF filtering and sorting outside the loop once
    subprocess.run(f"python3 {code_dir}/vcf_filter.py -v {vcf} -p {prefix} -o_dir . --remove_small_sv", shell=True, check=True)
    for typ in ["DEL", "INS"]:
        subprocess.run(f"vcf-sort {prefix}_{typ}_noXY.vcf > {prefix}_{typ}_noXY_sorted.vcf", shell=True, check=True)
        subprocess.run(f"bgzip -c {prefix}_{typ}_noXY_sorted.vcf > {prefix}_{typ}_noXY_sorted.vcf.gz", shell=True, check=True)
        subprocess.run(f"tabix -p vcf {prefix}_{typ}_noXY_sorted.vcf.gz", shell=True, check=True)

    # Define grid
    p_values = [round(x * 0.1, 1) for x in range(11)]
    r_values = list(range(0, 1001, 100))
    O_values = [round(x * 0.1, 1) for x in range(11)]

    # p-r grid
    tasks = [delayed(run_truvari)(prefix, vcf, ref, bed, bench_del, bench_ins, output_root, p=p, r=r) for p, r in product(p_values, r_values)]
    Parallel(n_jobs=n_thread)(tqdm(tasks, total=len(tasks), desc="Running p-r grid"))

    # p-O grid
    tasks = [delayed(run_truvari)(prefix, vcf, ref, bed, bench_del, bench_ins, output_root, p=p, O=O) for p, O in product(p_values, O_values)]
    Parallel(n_jobs=n_thread)(tqdm(tasks, total=len(tasks), desc="Running p-O grid"))

def main():
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--config','-conf', help = "each row is prefix,vcfpath ")
    parser.add_argument('--input_dir','-i')
    parser.add_argument('--tool','-tool')
    parser.add_argument('--work_dir','-w')
    parser.add_argument('--reference','-ref')
    parser.add_argument('--bench_dir','-bench')
    parser.add_argument('--n_thread','-t', type = int, default = 22 )
    # parser.add_argument('--delete_temp_file','-d', action='store_true')
    args = parser.parse_args()
    # config = args.config
    tool = args.tool
    input_dir = args.input_dir
    work_dir = os.path.abspath(args.work_dir)
    n_thread = args.n_thread
    ref = args.reference
    bench_dir = args.bench_dir
    # with open(config,'r') as f:
    #     for line in f:
    #         prefix, vcf = line.strip().split(',')
    #         if not os.path.exists(vcf):
    #             print(f"{vcf} does not exist!")
    #             exit()

    for vcf in glob.glob(f"{input_dir}/*_{tool}.vcf"):
        prefix = '_'.join(os.path.basename(vcf).split('_')[:2])
        process_one_vcf(vcf, prefix, work_dir, n_thread, ref, bench_dir)


if __name__ == "__main__":
    main()
