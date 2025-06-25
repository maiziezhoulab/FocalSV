import os
import subprocess
import sys
from subprocess import Popen
import argparse
from argparse import ArgumentParser
from GT_impute import *
from ONT_var_process import *

# Argument parser
parser = ArgumentParser(description="Filters structural variants (SVs) from a VCF file, calculates signature support, and performs genotype correction.",
	usage='use "python3 %(prog)s --help" for more information',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Adding arguments
#-------required
parser.add_argument("--bam_file", '-bam', type=str, help="Path to the input whole genome reads BAM file.")
parser.add_argument("--vcf_file", '-vcf',type=str, help="Path to the input FocalSV Candidate SV VCF file.")
parser.add_argument("--data_type", '-d',type=str, choices=["ONT", "Hifi", "CLR"], help="Data type (e.g., ONT, Hifi, CLR).")
parser.add_argument("--ref_file", '-r',type=str, help="Path to the reference genome file in FASTA format.")
parser.add_argument(
    "--chr_num",'-chr',
    type=str,
    choices = [str(i) for i in range(1,23)]+['wgs'],
    help="Chromosome number for targeted processing (use 'wgs' to enable whole genome scale post-processing).",
)
#------optional
parser.add_argument("--out_dir",'-o', type=str, help="Working directory for signature and corrected VCF outputs.", default= "./FocalSV_Final_VCF")
parser.add_argument("--num_threads",'-thread', type=int, help="Number of threads to use.", default = 10)
# parser.add_argument(
#     "--sigdir",'-sig',
#     type=str,
#     help="Pre-extracted reads signature (use 'None' to enable auto-extraction).",
# )
parser.add_argument(
    "--sigdir",'-sig',
    type=str,
    help="Pre-extracted reads signature.",
)


# Parse arguments
args = parser.parse_args()

# Validate inputs
if not os.path.isfile(args.bam_file):
    print(f"Error: BAM file '{args.bam_file}' does not exist.")
    sys.exit(1)
if not os.path.isfile(args.vcf_file):
    print(f"Error: VCF file '{args.vcf_file}' does not exist.")
    sys.exit(1)
if not os.path.isfile(args.ref_file):
    print(f"Error: Reference genome file '{args.ref_file}' does not exist.")
    sys.exit(1)

# Process arguments
bamfile = os.path.realpath(args.bam_file)
vcffile = os.path.realpath(args.vcf_file)
dtype = args.data_type
if dtype == 'HIFI':
    dtype = "Hifi"
wdir = os.path.realpath(args.out_dir)
t = args.num_threads
reference = os.path.realpath(args.ref_file)
chr_num = args.chr_num
sigdir = args.sigdir
# sigdir = None # aoto-extract enabled
#---------------test only-----------
# sigdir = "/data/maiziezhou_lab/CanLuo/FocalSV/LargeINDEL/HiFi_l2exp_test/out/reads_sig/"
# Print argument values for debugging or verification
print(f"BAM file: {bamfile}")
print(f"VCF file: {vcffile}")
print(f"Data type: {dtype}")
print(f"Working directory: {wdir}")
print(f"Threads: {t}")
print(f"Reference genome: {reference}")
print(f"Chromosome number: {chr_num}")
print(f"Signature directory: {sigdir}")



if chr_num!='wgs':
    chr_num = int(chr_num)
ft_vtype = "DEL"


# Derived variables
filtered_vcf = f"{wdir}/{os.path.basename(vcffile).replace('.vcf', '')}_filter_{ft_vtype}.vcf"
dtype_lc = dtype.lower()
gtdir = f"{wdir}/GT_Correction"
candidate_vcf = f"{wdir}/FocalSV_Candidate_SV.vcf"
candidate_vcf_ins = f"{wdir}/FocalSV_Candidate_SV_INS.vcf"
candidate_vcf_del = f"{wdir}/FocalSV_Candidate_SV_DEL.vcf"
cand_vcf_ont = filtered_vcf.replace(".vcf","_updated_GT.vcf")


final_vcf = f"{wdir}/FocalSV_Final_SV.vcf"
workdir = wdir
# Ensure necessary directories exist
os.makedirs(wdir, exist_ok=True)


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


# call reads signature

def call_sig(dtype,bamfile, sigdir, reference, chromosome):

    os.system("mkdir -p " + sigdir)

    if  dtype == "Hifi":
        para ='''--max_cluster_bias_INS      1000 \
        --diff_ratio_merging_INS    0.9 \
        --max_cluster_bias_DEL  1000 \
        --diff_ratio_merging_DEL    0.5 '''
    elif dtype == 'CLR':
        para = '''--max_cluster_bias_INS      100 \
        --diff_ratio_merging_INS    0.3 \
        --max_cluster_bias_DEL  200 \
        --diff_ratio_merging_DEL    0.5 '''
    else:
        para = '''--max_cluster_bias_INS      100 \
        --diff_ratio_merging_INS    0.3 \
        --max_cluster_bias_DEL  100 \
        --diff_ratio_merging_DEL    0.3 '''
    cmd = f'''python3 {code_dir}/Reads_Based_Scan/Reads_Based_Scan.py \
    {bamfile} \
    {reference} \
    {sigdir}/reads_draft_variants.vcf \
    {sigdir} \
    {para} \
    -chr {chromosome} -t {t} --genotype --retain_work_dir '''
    print(cmd)
    Popen(cmd, shell = True).wait()


# call sig
if sigdir is None:
	sigdir = f"{wdir}/reads_sig"
	call_sig(dtype,bamfile, sigdir, reference, chr_num)


reads_draft_vcf = sigdir + "/reads_draft_variants.vcf"

# Ensure necessary directories exist
os.makedirs(gtdir, exist_ok=True)

def run_command(command):
    """Helper function to run shell commands."""
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

# Step 1: Calculate signature support for each variant
run_command(f"""
    python3 {code_dir}/calculate_signature_support.py \
    -v {vcffile} \
    -ct {sigdir} \
    -w {wdir}
""")

# Step 2: Filter VCF by empirical signature coverage

run_command(f"""
    python3 {code_dir}/filter_vcf_by_sig_cov_insdel.py \
    -i {vcffile} \
    -d {dtype_lc} \
    -a volcano \
    -v {ft_vtype} \
    -w {wdir}
""")


# Step 3: Perform GT correction; extra process for ONT
if dtype == 'Hifi':
     
    # Step 3: Perform GT correction for deletions
    run_command(f"""
        python3 {code_dir}/correct_gt_del_real_data.py \
        -i {filtered_vcf} \
        -o {gtdir}/bnd_del_real.tsv \
        -bam {bamfile} \
        -sig {sigdir}/DEL.sigs \
        -t 7 \
        -d {dtype} -v DEL
    """)

    # Step 4: Perform GT correction for insertions
    run_command(f"""
        python3 {code_dir}/correct_gt_ins_real_data.py \
        -i {filtered_vcf} \
        -o {gtdir}/bnd_ins_real.tsv \
        -bam {bamfile} \
        -sig {sigdir}/INS.sigs \
        -t 7 \
        -d {dtype} -v INS
    """)

    # Step 5: Combine corrected VCFs and generate final VCF
    header_file = f"{workdir}/header"
    run_command(f"grep '#' {filtered_vcf} > {header_file}")

    # Concatenate the header and corrected VCF files, then sort the result
    run_command(f"""
        cat {header_file} {filtered_vcf}.newgt.DEL {filtered_vcf}.newgt.INS | vcf-sort > {final_vcf}
    """)

    print(f"Process completed. Final VCF saved at: {final_vcf}")

elif dtype == 'CLR':
    # -------- use reads inferred GT instead
    gt_impute(filtered_vcf, reads_draft_vcf, final_vcf, 1000, 0.5)
else:
    # -------- use reads inferred GT instead
    gt_impute(filtered_vcf, reads_draft_vcf, cand_vcf_ont, 1000, 0.5)

    # filter DEL & union INS
    final_process_ont(cand_vcf_ont,reads_draft_vcf, final_vcf )



