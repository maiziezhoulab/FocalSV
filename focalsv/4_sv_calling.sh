#!/bin/bash
set -euo pipefail

in_dir=$1       # Input base directory
out_dir=$2      # Output base directory
bam_file=$3     # BAM file
ref_file=$4 
dtype=$5        # Data type string
code_dir=$6     # Code directory (contains rename_fa.py and DipPAV_variant_call.py)
chr=$7          # Chromosome number, e.g., 21
t=$8

# Ensure output dirs exist
mkdir -p "${out_dir}/SV/chr${chr}"

# Step 1: Concatenate HP1.fa and HP2.fa
cat ${in_dir}/regions/*/HP1.fa > ${out_dir}/chr${chr}_HP1.fa
cat ${in_dir}/regions/*/HP2.fa > ${out_dir}/chr${chr}_HP2.fa

# Step 2: Rename fasta files
python3 ${code_dir}/rename_fa.py \
    -i ${out_dir}/chr${chr}_HP1.fa \
    -o ${out_dir}/chr${chr}_HP1_new.fa \
    -hp 1

python3 ${code_dir}/rename_fa.py \
    -i ${out_dir}/chr${chr}_HP2.fa \
    -o ${out_dir}/chr${chr}_HP2_new.fa \
    -hp 2

# Step 3: Run DipPAV variant caller
python3 ${code_dir}/4_sv_calling/Dippav/DipPAV_variant_call.py \
    -rbam ${bam_file} \
    -ref ${ref_file} \
    -hp1 ${out_dir}/chr${chr}_HP1_new.fa \
    -hp2 ${out_dir}/chr${chr}_HP2_new.fa \
    -o ${out_dir}/SV/chr${chr} \
    -d ${dtype} \
    -chr ${chr} \
    -t $t

python3 ${code_dir}/5_post_processing/FocalSV_Filter_GT_Correct.py \
        -bam ${bam_file} \
        -r ${ref_file} \
        -chr ${chr} \
        -o ${out_dir} -thread $t \
        -d ${dtype}

# cp ${out_dir}/SV/chr${chr}/final_vcf/dippav_variant_no_redundancy.vcf  ${out_dir}/FocalSV_Candidate_SV.vcf
