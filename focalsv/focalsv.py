#!/usr/bin/env python3
from subprocess import check_call, CalledProcessError
from argparse import ArgumentParser
import os
import sys
from utils import setup_logging  # Assuming setup_logging exists in utils.py
from main import wrapper
import glob
script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/"
__author__ = "Maizie&Jamie&Can@Vandy"

parser = ArgumentParser(description="Author: maiziezhoulab@gmail.com\n", usage='use "python3 %(prog)s --help" for more information')

# General inputs
parser.add_argument('--bam_file', '-bam', help="BAM file", required=True)
parser.add_argument('--ref_file', '-r', help="Reference FASTA file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region, use 0 to select all autosomal chromsomes", required=True)

# For single region
parser.add_argument('--region_start', '-S', type=int, help="Target region starting index", required=False)
parser.add_argument('--region_end', '-E', type=int, help="Target region ending index", required=False)

# For multiple regions
parser.add_argument('--target_bed', '-target_bed', help="BED file with multiple target regions", required=False)

# Defaulted inputs
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results, default = ./RegionBased_results", default="./RegionBased_results")
parser.add_argument('--data_type', '-d', help="HIFI/CLR/ONT")

# Process information
parser.add_argument('--num_cpus', '-cpu', type=int, help="Number of CPUs, default = 10", default=10)
parser.add_argument('--num_threads', '-thread', type=int, help="Number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--early_threads', '-ethread', type=int, help="Number of threads for cropping bam and phasing, default = 8 (recommended)", default=8)

args = parser.parse_args()


bam_file = args.bam_file
ref_file = args.ref_file
chr_num = args.chr_num
region_start = args.region_start
region_end = args.region_end
out_dir = args.out_dir
data_type = args.data_type
num_threads = args.num_threads
early_threads = args.early_threads
num_cpus = args.num_cpus
target_bed = args.target_bed


def split_bed(bed_file, out_dir):
    if not os.path.exists(out_dir):
        os.system("mkdir -p "+out_dir)
    cmd = f"""awk '{{print > "{out_dir}/" $1 ".bed"}}' {bed_file}"""
    print(cmd)
    os.system(cmd)

def extract_bed(bed_file, out_dir, chrom):
    if not os.path.exists(out_dir):
        os.system("mkdir -p "+out_dir)
    cmd = f"""grep -w {chrom} {bed_file} > {out_dir}/{chrom}.bed"""
    print(cmd)
    os.system(cmd)

def merge_vcf(in_dir, chr_num_list):

    cmd = f'''grep '#' {in_dir}/chr{chr_num_list[0]}/Final_SV/FocalSV_Final_SV.vcf > {in_dir}/FocalSV_Final_SV.vcf
    grep -v '#' {in_dir}/chr*/Final_SV/FocalSV_Final_SV.vcf|vcf-sort >> {in_dir}/FocalSV_Final_SV.vcf'''
    os.system(cmd)


if target_bed is None:
    wrapper(bam_file,ref_file,chr_num,out_dir,
                data_type,target_bed,
                region_start,region_end,num_threads,early_threads,num_cpus )

else:
    if chr_num == 0:
        split_bed(target_bed, out_dir)
        chr_num_list = [int(os.path.basename(bed).split('.')[0][3:]) for bed in glob.glob(f"{out_dir}/chr*.bed")]
    else:
        extract_bed(target_bed, out_dir, f"chr{chr_num}")
        chr_num_list = [chr_num]

    for use_chr_num in chr_num_list:
        wrapper(bam_file,ref_file,use_chr_num,out_dir,
                data_type,f"{out_dir}/chr{use_chr_num}.bed",
                region_start,region_end,num_threads,early_threads,num_cpus )

    if chr_num == 0:
        merge_vcf(out_dir, chr_num_list)

    

