import os
import subprocess
import sys
from subprocess import Popen

def print_help():
    # Get the current script's file name
	script_name = os.path.basename(__file__)
	# print(f"The script name is: {script_name}")

	"""Prints help message for the script usage."""
	help_message = f"""
Usage:
    python3 {script_name} <bamfile> <vcffile> <dtype> <wdir> <t> <reference> <chr_num> <sigdir>

Description:
    This script filters structural variants (SVs) from a VCF file, calculates signature support,
    and performs genotype correction based on BAM alignments and signatures.

Arguments:
    bamfile    : Path to the input BAM file.
    vcffile    : Path to the input VCF file.
    dtype      : Data type (e.g., ONT, Hifi, CLR).
    wdir       : Working directory for signature and corrected VCF outputs.
    t          : Number of threads to use.
    reference  : Path to the reference genome file in FASTA format.
    chr_num    : Chromosome number for targeted processing (use None to enable whole genome scale post processing).
    sigdir	   : preextracted reads signature (use None to enable auto-extract).

Example:
    python3 {script_name} input.bam input.vcf ONT /path/to/workdir 5 reference.fasta 1 None

Notes:
    - Ensure all required Python dependencies and tools (samtools, bcftools, etc.) are installed.
    - The script generates a corrected VCF file and combines all processed VCFs into a final version.
    """
	print(help_message)

# Check if the user requested help
if len(sys.argv) < 9 or sys.argv[1] in ["-h", "--help"]:
    print_help()
    sys.exit(1)
    
# Input arguments
bamfile = os.path.realpath(sys.argv[1])
vcffile = os.path.realpath(sys.argv[2])
dtype = sys.argv[3]  # ONT HiFi CLR
# sigdir = os.path.realpath(sys.argv[4])
wdir = os.path.realpath(sys.argv[4])
t = sys.argv[5]
reference = os.path.realpath(sys.argv[6])
chr_num = sys.argv[7]
sigdir = sys.argv[8]

ft_vtype = "DEL"

# Derived variables
filtered_vcf = f"{wdir}/{os.path.basename(vcffile).replace('.vcf', '')}_filter_{ft_vtype}.vcf"
dtype_lc = dtype.lower()
gtdir = f"{wdir}/GT_Correction"


final_vcf = f"{wdir}/FocalSV_final_variants.vcf"
workdir = wdir
# Ensure necessary directories exist
os.makedirs(wdir, exist_ok=True)


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


# call reads signature

def call_sig(dtype,bamfile, sigdir, reference, chr_num):
    if chr_num == 'None':
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
