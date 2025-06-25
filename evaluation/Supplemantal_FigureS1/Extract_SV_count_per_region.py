import pandas as pd
import os
import gzip
from joblib import Parallel, delayed
from intervaltree import IntervalTree

# Function to build an interval tree for fast region lookup
def build_interval_tree(regions):
    tree = IntervalTree()
    for _, row in regions.iterrows():
        tree[row["start"]:row["end"]] = row["chr"]
    return tree

# Function to count SVs in a VCF file within sampled regions and aggregate by chromosome
def count_sv_in_vcf(vcf_file, regions, tool):
    """Counts the number of SV entries in a given VCF file that fall within the sampled regions and aggregates per chromosome."""
    if not os.path.exists(vcf_file):
        print(f"Warning: {vcf_file} not found. Skipping...")
        return None

    print(f"Processing {vcf_file}...")

    # Read the VCF file correctly, skipping comment lines
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                break
    
    # Read VCF data (ensuring correct columns)
    df = pd.read_csv(vcf_file, comment="#", sep="\t", names=header)

    if "Dipcall" not in tool:
        # Non-Dipcall tools: Extract SVTYPE and filter INS/DEL only
        df["svtype"] = df["INFO"].str.extract(r"SVTYPE=(INS|DEL)")[0]
        df = df.dropna(subset=["svtype"])
        df["svlen"] = df["INFO"].str.extract(r"SVLEN=(-?\d+)")[0].astype(float).abs()
        df = df[df["svlen"] >= 50]  # Filter SVs with abs(SVLEN) >= 50
    else:
        # Dipcall: Infer SV length from REF and ALT
        df["svlen"] = (df["ALT"].str.len() - df["REF"].str.len()).abs()
        df = df[df["svlen"] >= 50]  # Keep only SVs where |len(REF) - len(ALT)| >= 50

    # Count SVs per region
    sv_counts = []
    for _, region in regions.iterrows():
        count = df[(df["#CHROM"] == region["chr"]) & (df["POS"].between(region["start"], region["end"]))].shape[0]
        sv_counts.append([region["chr"], region["start"], region["end"], count])

    sv_counts_df = pd.DataFrame(sv_counts, columns=["chr", "start", "end", "sv_count"])

    return sv_counts_df

# Function to count SVs for FocalSV directly from the provided file
def count_sv_for_focalsv(focal_vcf_file, regions):
    """Counts SVs per region from the FocalSV VCF file."""
    print(f"Processing FocalSV file: {focal_vcf_file}")

    if not os.path.exists(focal_vcf_file):
        print(f"Warning: {focal_vcf_file} not found. Skipping FocalSV processing...")
        return None

    records = []
    with gzip.open(focal_vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            records.append((chrom, pos))

    # Convert to DataFrame
    df = pd.DataFrame(records, columns=["chr", "POS"])

    # Count per region
    sv_counts = []
    for _, region in regions.iterrows():
        count = df[(df["chr"] == region["chr"]) & (df["POS"].between(region["start"], region["end"]))].shape[0]
        sv_counts.append([region["chr"], region["start"], region["end"], count])

    return pd.DataFrame(sv_counts, columns=["chr", "start", "end", "sv_count"])

import argparse
parser = argparse.ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
# parser.add_argument('--output_path','-o')
# parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
vcf_dir = args.input_dir
# output_path = args.output_path
# n_thread = args.n_thread


# Load sampled regions BED file
bed_file = "1k_SV_negative_regions.bed"
print("Loading filtered regions...")
bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "end"])
print(f"Loaded {len(bed_df)} sampled regions.")

# Directory containing VCF files
# vcf_dir = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/CanLuo/FocalSV/GR_Revision/VCF_Collection"

# Output directories
output_dirs = {
    "clr": "clr_all_tools",
    "hifi": "hifi_all_tools",
    "ont": "ont_all_tools"
}
for dir_path in output_dirs.values():
    os.makedirs(dir_path, exist_ok=True)

# List of tools to process
other_tools = ['FocalSV',"cuteSV", "Dipcall",
               "PAV",  "PBSV", "Sniffles2", "SVIM-asm", "SVIM",'sawfish','SKSV']

def process_tool(data_type, prefix, tool):
    if tool == "FocalSV":
        focal_vcf_file = f"FocalSV_in_sampled_region/{data_type}_FocalSV.vcf.gz"
        region_df = count_sv_for_focalsv(focal_vcf_file, bed_df)
    else:
        vcf_file = f"{prefix}{tool}.vcf"
        vcf_path = os.path.join(vcf_dir, vcf_file)
        if not os.path.exists(vcf_path):
            print(f"Warning: {vcf_path} not found. Skipping...")
            return
        region_df = count_sv_in_vcf(vcf_path, bed_df, tool)

    if region_df is None:
        return

    output_csv = os.path.join(output_dirs[data_type], f"sv_counts_by_region_{tool}.csv")
    region_df.to_csv(output_csv, sep=",", index=False)
    print(f"Per-region SV counts saved to {output_csv}")

# Run parallel processing
cpu = os.cpu_count()
Parallel(n_jobs=cpu)(
    delayed(process_tool)(data_type, prefix, tool)
    for data_type, prefix in zip(["hifi","clr",'ont'], ["Hifi_L1_",'CLR_L1_','ONT_L1_'])
    for tool in other_tools  # Ensure FocalSV is processed first
)

print("All processing completed.")
