from subprocess import Popen

import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
# parser.add_argument('--hp1fa','-hp1')
# parser.add_argument('--hp2fa','-hp2')
# parser.add_argument('--indelvcf','-vcf')
parser.add_argument('--bam_file','-bam')
parser.add_argument('--bed_file','-bed')
# parser.add_argument('--bed_file','-bed')
parser.add_argument('--ref_file','-ref')
parser.add_argument('--data_type','-d', choices=['HIFI','CLR','ONT'])
parser.add_argument('--out_dir','-o')
parser.add_argument('--num_threads','-thread', type = int, default = 8 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_dir = args.input_dir
# hp1fa = args.hp1fa
# hp2fa = args.hp2fa
# indelvcf = args.indelvcf
bamfile = args.bam_file
out_dir = args.out_dir
bed_file = args.bed_file
reference = args.ref_file
n_thread = args.num_threads
datatype = args.data_type

import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

def split_bed(bed_file, out_dir):
    if not os.path.exists(out_dir):
        os.system("mkdir -p " + out_dir)

    for svtype in ['DUP','INV','TRA']:
        cmd = f"grep {svtype} {bed_file} > {out_dir}/{svtype}_regions.bed"
        os.system(cmd)

bed_file_tra = f"{out_dir}/TRA_regions.bed"
bed_file_inv = f"{out_dir}/INV_regions.bed"
bed_file_dup = f"{out_dir}/DUP_regions.bed"

split_bed(bed_file, out_dir)
cmd = f'''python3 {code_dir}/call_DUP.py \
    -i {input_dir} \
    -bam {bamfile} \
    -bed {bed_file_dup} \
    -ref {reference} \
    -d {datatype} \
    -o {out_dir}/DUP \
    -t {n_thread}'''
Popen(cmd, shell = True).wait()

cmd = f'''python3 {code_dir}/Reads_Based_INV_Call.py \
    -bam {bamfile} \
    -bed {bed_file_inv} \
    -o {out_dir}/INV \
    -t {n_thread}'''
print(cmd)
Popen(cmd, shell = True).wait()

cmd = f'''python3 {code_dir}/Reads_Based_TRA_Call.py \
    -bam {bamfile} \
    -bed {bed_file_tra} \
    -o {out_dir}/TRA \
    -t {n_thread}'''
Popen(cmd, shell = True).wait()

out_vcf = out_dir+"/FocalSV_TRA_INV_DUP.vcf"
cmd = f'''cat {out_dir}/DUP/DUP.vcf |grep '#' > {out_vcf};
cat {out_dir}/DUP/DUP.vcf {out_dir}/INV/INV.vcf {out_dir}/TRA/TRA.vcf|grep -v '#'|vcf-sort >> {out_vcf}; sed -i 's/svim_asm/focalsv/g' {out_vcf}
'''
Popen(cmd, shell = True).wait()