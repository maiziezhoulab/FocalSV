import os
import sys
from argparse import ArgumentParser
from utils import setup_logging  # Import the setup_logging function

script_path = os.path.dirname(os.path.abspath(__file__))
code_path = script_path + "/"

# Set up argument parsing
parser = ArgumentParser(description="Crop target region from BAM file:")
parser.add_argument('--bam_file', '-bam', help="Input BAM file", required=True)
parser.add_argument('--chr_num', '-chr', type=int, help="Chromosome number for target region", required=True)

# Accept either start/end region or a BED file
parser.add_argument('--region_start', '-S', type=int, help="Start index of target region", required=False)
parser.add_argument('--region_end', '-E', type=int, help="End index of target region", required=False)
parser.add_argument('--target_bed', '-target_bed', help="BED file with multiple target regions", required=False)

parser.add_argument('--out_dir', '-o', help="Output directory", required=True)

args = parser.parse_args()

# Function to crop the BAM file for a specific region
def crop_single_region(bam_file, chr_num, region_start, region_end, out_dir, logger):
    out_dir_region = os.path.join(out_dir, "regions")
    os.makedirs(out_dir_region, exist_ok=True)
    
    region = f"chr{chr_num}:{region_start}-{region_end}"
    output_dir = f"{out_dir_region}/Region_chr{chr_num}_S{region_start}_E{region_end}"
    os.makedirs(output_dir, exist_ok=True)

    cropped_file = f"{output_dir}/region.bam"

    logger.info(f"Starting to crop BAM file for region {region}")
    
    # Call samtools to crop and index the BAM file
    cmd1 = f"samtools view {bam_file} {region} -Sb > {cropped_file}"
    cmd2 = f"samtools index {cropped_file}"
    
    try:
        os.system(cmd1)
        os.system(cmd2)
        logger.info(f"*** BAM region cropped and indexed: {region} ***")
        logger.info(f"{cmd1}\n{cmd2}\n")
    except Exception as e:
        logger.error(f"Error during cropping or indexing: {e}")

# Function to crop the BAM file for multiple regions from a BED file
def crop_multiple_regions(bam_file, chr_num, bed_file, out_dir, logger):
    out_dir_region = os.path.join(out_dir, "regions")
    os.makedirs(out_dir_region, exist_ok=True)

    logger.info(f"Starting to crop BAM file for multiple regions from {bed_file}")
    
    # Process each line from the BED file
    with open(bed_file, "r") as bed:
        for line in bed:
            cols = line.strip().split()
            bed_chr, start, end = cols[0], int(cols[1]), int(cols[2])
            if bed_chr != f"chr{chr_num}":
                continue

            region = f"{bed_chr}:{start}-{end}"
            output_dir = f"{out_dir_region}/Region_{bed_chr}_S{start}_E{end}"
            os.makedirs(output_dir, exist_ok=True)
            
            cropped_file = f"{output_dir}/region.bam"

            # Call samtools to crop and index the BAM file
            cmd1 = f"samtools view {bam_file} {region} -Sb > {cropped_file}"
            cmd2 = f"samtools index {cropped_file}"
            
            try:
                os.system(cmd1)
                os.system(cmd2)
                logger.info(f"*** BAM region cropped and indexed: {region} ***")
                logger.info(f"{cmd1}\n{cmd2}\n")
            except Exception as e:
                logger.error(f"Error during cropping or indexing for region {region}: {e}")
                

if __name__ == "__main__":
    # Get arguments
    bam_file = args.bam_file
    chr_num = args.chr_num
    region_start = args.region_start
    region_end = args.region_end
    target_bed = args.target_bed
    out_dir = args.out_dir

    # Initialize the logger
    logger = setup_logging("1_CROP_BAM", out_dir)
    
    # Validate input: Ensure that either start/end or target_bed is provided, but not both
    if (region_start is not None and region_end is not None) and target_bed:
        logger.error("Specify either --region_start/--region_end or --target_bed, but not both.")
        sys.exit(1)
    elif not target_bed and (region_start is None or region_end is None):
        logger.error("You must specify both --region_start and --region_end if --target_bed is not provided.")
        sys.exit(1)

    # Run crop_bam with logging
    if target_bed:
        crop_multiple_regions(bam_file, chr_num, target_bed, out_dir, logger)
    else:
        crop_single_region(bam_file, chr_num, region_start, region_end, out_dir, logger)
