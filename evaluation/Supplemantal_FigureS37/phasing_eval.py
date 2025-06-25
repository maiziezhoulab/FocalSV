import os
import subprocess
import re
import logging
from concurrent.futures import ProcessPoolExecutor

import argparse
parser = argparse.ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
parser.add_argument('--bench_dir','-b')
args = parser.parse_args()
BASE_PATH = args.input_dir
bench_dir = args.bench_dir



# === Logging Setup ===
logging.basicConfig(
    filename="phasing_eval.log",
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# === Updated Constants ===
# BASE_PATH = "/data/maiziezhou_lab/Jamiezhou/region_based/Review/ps_hp"
GOLD_VCF = f"{bench_dir}/HG002_GRCh37_GIAB_highconf_phased_SNPs.vcf.gz"
EVAL_SCRIPT = "./measure_phasing_performance.pl"
TECHS = ["hifi_l1", "clr_l1", "ont_l1"]

def run(cmd):
    logger.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.cmd}\nExit code: {e.returncode}\nStderr: {e.stderr.decode()}")
        raise

def extract_coords(region_folder):
    match = re.search(r"S(\d+)_E(\d+)", region_folder)
    if match:
        return int(match.group(1)), int(match.group(2))
    else:
        raise ValueError(f"Cannot extract coordinates from {region_folder}")

def vcf_has_variants(vcf_path):
    try:
        output = subprocess.check_output(f"bcftools view {vcf_path} | grep -v '^#' | head -n 1", shell=True)
        return bool(output.strip())
    except subprocess.CalledProcessError:
        return False

def evaluate_region(args):
    tech, chrom, region_folder = args
    region_path = os.path.join(BASE_PATH, tech, chrom, "regions", region_folder)
    vcf_path = os.path.join(region_path, "region_phased.vcf")
    if not os.path.exists(vcf_path):
        logger.warning(f"Missing VCF file: {vcf_path}")
        return

    try:
        start, end = extract_coords(region_folder)
        compressed_vcf = vcf_path + ".gz"
        cropped_gold = os.path.join(region_path, "gold_cropped.vcf.gz")
        output_dir = os.path.join(region_path, "phasing_eval")

        # Step 1: Compress and index Longshot VCF if needed
        if not os.path.exists(compressed_vcf):
            run(f"bgzip -c {vcf_path} > {compressed_vcf}")
        run(f"tabix -p vcf {compressed_vcf} -f")

        # Step 2: Crop gold standard VCF
        run(f"bcftools view -r {chrom}:{start}-{end} {GOLD_VCF} -O z -o {cropped_gold}")
        run(f"tabix -p vcf {cropped_gold} -f")

        # Step 3: Filter both to **biallelic SNPs**
        filtered_longshot = os.path.join(region_path, "region_phased_snps.vcf.gz")
        filtered_gold = os.path.join(region_path, "gold_cropped_snps.vcf.gz")

        run(f"bcftools view -v snps -m2 -M2 {compressed_vcf} -O z -o {filtered_longshot}")
        run(f"tabix -p vcf {filtered_longshot} -f")

        run(f"bcftools view -v snps -m2 -M2 {cropped_gold} -O z -o {filtered_gold}")
        run(f"tabix -p vcf {filtered_gold} -f")

        # Step 4: Check for content
        if vcf_has_variants(filtered_longshot) and vcf_has_variants(filtered_gold):
            run(f"{EVAL_SCRIPT} -i {filtered_longshot} -r {filtered_gold} -o {output_dir}")
            logger.info(f"✅ Completed evaluation: {output_dir}")
        else:
            logger.warning(f"⚠️ Skipping {region_path} — one or both SNP VCFs are empty.")

    except Exception as e:
        logger.error(f"[ERROR] Failed processing {region_path}: {e}", exc_info=True)

def gather_all_tasks():
    tasks = []
    for tech in TECHS:
        tech_path = os.path.join(BASE_PATH, tech)
        for chrom in [f"chr{i}" for i in range(1, 23)]:
            region_base = os.path.join(tech_path, chrom, "regions")
            if not os.path.exists(region_base):
                logger.warning(f"Missing: {region_base}")
                continue
            for region_folder in os.listdir(region_base):
                region_path = os.path.join(region_base, region_folder)
                vcf_file = os.path.join(region_path, "region_phased.vcf")
                if os.path.isfile(vcf_file):
                    tasks.append((tech, chrom, region_folder))
                    logger.info(f"Found: {vcf_file}")
                else:
                    logger.warning(f"No VCF in: {region_path}")
    return tasks

# === Main ===
if __name__ == "__main__":
    logger.info("🚀 Starting phasing evaluation pipeline...")

    all_tasks = gather_all_tasks()
    logger.info(f"📦 Total regions to evaluate: {len(all_tasks)}")

    max_workers = os.cpu_count() // 2 or 4
    logger.info(f"🧵 Using {max_workers} parallel workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(evaluate_region, all_tasks)

    logger.info("✅ All evaluations complete.")
