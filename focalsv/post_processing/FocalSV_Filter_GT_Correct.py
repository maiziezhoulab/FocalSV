import os
import subprocess
import sys
from subprocess import Popen
import argparse
from argparse import ArgumentParser

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
#     "sigdir",
#     type=str,
#     help="Pre-extracted reads signature (use 'None' to enable auto-extraction).",
# )

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
# sigdir = args.sigdir
sigdir = None # aoto-extract enabled
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
final_vcf = f"{wdir}/FocalSV_Final_SV.vcf"
workdir = wdir
# Ensure necessary directories exist
os.makedirs(wdir, exist_ok=True)


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


# call reads signature

def call_sig(dtype,bamfile, sigdir, reference, chr_num):
    if chr_num == 'wgs':
        bed_para = " "
    else:
        # get chromosome length
        fai_file = reference + ".fai"
        with open(fai_file,'r') as f:
            chr_len = int(f.readlines()[chr_num -1 ].split()[1])
        bed_file = sigdir + "/sample.bed"
        with open(bed_file, 'w') as f:
            f.write("chr"+str(chr_num)+"\t1\t"+str(chr_len)+'\n')
        
        bed_para = " -include_bed " + bed_file

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
    cmd = f'''python3 {code_dir}/sig_extract.py \
    {bamfile} \
    {reference} \
    {sigdir} \
    {para} \
    {bed_para} -t {t}'''
    print(cmd)
    Popen(cmd, shell = True).wait()


# call sig
if sigdir == 'None':
	sigdir = f"{wdir}/reads_sig"
	call_sig(dtype,bamfile, sigdir, reference, chr_num)


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
if dtype == 'ONT':
	run_command(f"""
		python3 {code_dir}/filter_vcf_by_sig_cov_insdel.py \
		-i {vcffile} \
		-d {dtype_lc} \
		-a volcano \
		-v {ft_vtype} \
        -w {wdir}
	""")
else:
    cmd = f"cp {vcffile} {filtered_vcf}"
    os.system(cmd)

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
