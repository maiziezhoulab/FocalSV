import os
from argparse import ArgumentParser
import sys
from utils import setup_logging
import allel

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/"

# Argument Parsing
parser = ArgumentParser(description="Find target SVs and produce BED regions:")
parser.add_argument('--vcf_file', '-v', help="VCF file called by FreeBayes", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target variant or region", required=True)
parser.add_argument('--out_dir', '-o', help="Directory to store assembly results", default="./RegionBased_results")
parser.add_argument('--bed_file', '-bed', help="BED file used to filter SVs", required=True)
parser.add_argument('--svlen_threshold', '-sv_thresh', type=int, help="Threshold number for initially filtering SVs", default=50)
parser.add_argument('--flanking', '-flank', type=int, help="Length of flanking region", default=50000)
args = parser.parse_args()

def get_chr_length(chr_num, chr_length_file=code_path+"0_chr_length.txt"):
    """
    Function to fetch the length of a chromosome from a file.
    """
    try:
        with open(chr_length_file, "r") as f:
            for line in f:
                chrom, length = line.strip().split()
                if chrom == f"chr{chr_num}":
                    return int(length)
    except FileNotFoundError:
        raise FileNotFoundError(f"Chromosome length file {chr_length_file} not found.")
    
    raise ValueError(f"Chromosome chr{chr_num} not found in {chr_length_file}.")

def produce_bed_file(out_dir, chr_num, FLANK, logger):
    """
    Function to produce BED file for regions around the target SVs.
    """
    output_bed = os.path.join(out_dir, "target_regions.bed")
    
    # Ensure the region directory exists
    out_dir_region = os.path.join(out_dir, "regions")
    os.makedirs(out_dir_region, exist_ok=True)
    
    # Get chromosome length
    LEN = get_chr_length(chr_num)
    
    # Read the VCF file
    callset = allel.read_vcf(vcf_file, fields=['variants/CHROM', 'variants/POS', 'variants/SVTYPE', 'variants/SVLEN'])
    total_snv = len(callset['variants/CHROM'])
    
    bed_lines = []
    
    # Process each variant
    for i in range(total_snv):
        chr = callset['variants/CHROM'][i]
        pos = callset['variants/POS'][i]
        sv_type = callset['variants/SVTYPE'][i]
        
        # Calculate the start and end based on SV type
        if sv_type == "INS":
            start = max(pos - FLANK, 0)
            end = min(pos + FLANK, LEN)
        else:
            sv_len = abs(callset['variants/SVLEN'][i])
            start = max(pos - FLANK, 0)
            end = min(pos + sv_len + FLANK, LEN)
        
        # Append the BED format line
        bed_lines.append(f"{chr}\t{start}\t{end}\n")
    
    # Write the BED file
    with open(output_bed, 'w') as output:
        output.writelines(bed_lines)
    
    logger.info(f"BED file written to {output_bed}")

def filter_target_sv(vcf_file, chr_num, out_dir, bed_file, svlen_threshold, flank, logger):
    """
    Filter target SVs from the VCF file based on chromosome number, BED file regions, and SV length threshold.
    """
    target_svs = os.path.join(out_dir, "target_sv.vcf")
    logger.info(f"Filtering target SVs from VCF: {vcf_file} for chr{chr_num} with SV length threshold: {svlen_threshold}")
    
    bed_positions = []
    try:
        with open(bed_file, "r") as bed:
            for line in bed:
                bed_data = line.split()
                if bed_data[0] == f"chr{chr_num}":
                    bed_positions.append((int(bed_data[1]), int(bed_data[2])))
    except FileNotFoundError:
        logger.error(f"BED file {bed_file} not found.")
        sys.exit(1)

    insertions, deletions = 0, 0

    try:
        with open(vcf_file, "r") as vcf, open(target_svs, "w") as target_vcf:
            for line in vcf:
                if line.startswith("#"):
                    target_vcf.write(line)
                else:
                    snv = line.split()
                    sv_chr, sv_pos = snv[0], int(snv[1])
                    sv_ref, sv_alt = snv[3], snv[4]
                    svlen = len(sv_alt) - len(sv_ref)
                    
                    if sv_chr == f"chr{chr_num}" and abs(svlen) >= svlen_threshold:
                        for bed_start, bed_end in bed_positions:
                            if bed_start <= sv_pos <= bed_end and snv[6] == "PASS":
                                if svlen > 0:
                                    insertions += 1
                                else:
                                    deletions += 1
                                target_vcf.write(line)
                                break

        logger.info(f"Filtering complete. Insertions: {insertions}, Deletions: {deletions}")
    except FileNotFoundError:
        logger.error(f"VCF file {vcf_file} not found.")
        sys.exit(1)

    produce_bed_file(out_dir, chr_num, flank, logger)

if __name__ == "__main__":
    vcf_file = args.vcf_file
    chr_num = args.chr_num
    out_dir = args.out_dir
    bed_file = args.bed_file
    svlen_threshold = args.svlen_threshold
    flank = args.flanking
    
    # Initialize logger
    logger = setup_logging("0_FILTER_SV", out_dir)

    try:
        filter_target_sv(vcf_file, chr_num, out_dir, bed_file, svlen_threshold, flank, logger)
        logger.info("Process completed successfully")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)