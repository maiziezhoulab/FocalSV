import os
import pysam
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

import argparse
parser = argparse.ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--auto_dir','-auto')
parser.add_argument('--target_dir','-target')
args = parser.parse_args()
auto_base = args.auto_dir
target_base = args.target_dir



# === Configuration ===
auto_paths = {
    "hifi_l1": f"{auto_base}/Hifi_L1",
    "clr_l1": f"{auto_base}/clr_L1",
    "ont_l1": f"{auto_base}/ont_L1"
}

target_paths = {
    "hifi_l1": f"{target_base}/hifi_l1",
    "clr_l1": f"{target_base}/clr_l1",
    "ont_l1": f"{target_base}/ont_l1"
}

bin_edges = np.linspace(0, 100, 30)
bar_colors = {
    "hifi_l1": "#89A8B2",
    "clr_l1": "#B3C8CF",
    "ont_l1": "#E5E1DA"
}
mean_color = "#d62728"
median_color = "#1f77b4"

# === BAM processing ===
def process_bam(region_bam):
    try:
        if not os.path.exists(region_bam):
            return None
        with pysam.AlignmentFile(region_bam, "rb") as bam:
            total_reads = 0
            phased_reads = 0
            unique_tags = set()
            for read in bam.fetch(until_eof=True):
                total_reads += 1
                try:
                    ps = read.get_tag("PS")
                    hp = read.get_tag("HP")
                    phased_reads += 1
                    unique_tags.add((ps, hp))
                except KeyError:
                    continue
            if total_reads > 0:
                percent_phased = (phased_reads / total_reads) * 100
                return (percent_phased, len(unique_tags))
    except Exception as e:
        print(f"❌ Error reading {region_bam}: {e}")
    return None

# === BAM file collector ===
def collect_bam_paths(root_dir, platform, mode):
    bam_paths = []
    for chr_num in range(1, 23):
        if mode == "auto":
            if platform == "hifi_l1":
                chr_path = os.path.join(root_dir, f"FocalSV_results_HIFI_L1_chr{chr_num}", "regions")
            else:
                chr_path = os.path.join(root_dir, f"chr{chr_num}_output", "regions")
        else:  # target mode
            chr_path = os.path.join(root_dir, f"chr{chr_num}", "regions")

        if not os.path.exists(chr_path):
            continue
        for region in os.listdir(chr_path):
            region_bam = os.path.join(chr_path, region, "region_phased.bam")
            bam_paths.append(region_bam)
    return bam_paths

# === Load data ===
def load_phased_data(paths_dict, mode):
    data = {}
    for platform, path in paths_dict.items():
        print(f"\n🔍 Processing {mode.upper()} - {platform.upper()}")
        bam_paths = collect_bam_paths(path, platform, mode)
        print(f"🧬 Found {len(bam_paths)} BAM files")

        with Pool(processes=cpu_count() // 2 or 4) as pool:
            results = pool.map(process_bam, bam_paths)

        results = [r for r in results if r is not None]
        if not results:
            print(f"⚠️ No valid results for {platform}")
            continue

        phased_percentages, _ = zip(*results)
        data[platform] = phased_percentages
    return data

# === Combined plot ===
def plot_combined_histogram(target_data, auto_data):
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), sharex='col', constrained_layout=True)
    platforms = ["hifi_l1", "clr_l1", "ont_l1"]
    row_labels = ["a", "b", "c"]
    col_titles = ["Target", "Auto"]

    for row_idx, platform in enumerate(platforms):
        for col_idx, data in enumerate([target_data, auto_data]):
            ax = axes[row_idx, col_idx]
            values = data.get(platform, [])
            if not values:
                ax.set_visible(False)
                continue

            mean_val = np.mean(values)
            median_val = np.median(values)

            ax.hist(values, bins=bin_edges, color=bar_colors[platform], edgecolor='black', alpha=0.9)
            ax.axvline(mean_val, color=mean_color, linestyle='--', linewidth=2)
            ax.axvline(median_val, color=median_color, linestyle=':', linewidth=2)

            ax.text(0.97, 0.90, f"Mean: {mean_val:.1f}%", transform=ax.transAxes,
                    ha='right', va='top', fontsize=11, color=mean_color)
            ax.text(0.97, 0.82, f"Median: {median_val:.1f}%", transform=ax.transAxes,
                    ha='right', va='top', fontsize=11, color=median_color)

            if col_idx == 0:
                ax.set_ylabel(platform.upper(), fontsize=13)
                ax.annotate(row_labels[row_idx], xy=(-0.15, 1.00), xycoords='axes fraction',
                            fontsize=15, fontweight='bold', ha='left', va='top')

            ax.tick_params(axis='both', which='major', labelsize=11)
            ax.grid(axis='y', linestyle='--', alpha=0.5)

    for ax in axes[2]:
        ax.set_xlabel("Phased Read Percentage (%)", fontsize=13)

    for col_idx, title in enumerate(col_titles):
        axes[0, col_idx].set_title(title, fontsize=15, fontweight='bold')

    fig.suptitle("Phased Read Percentage by Platform and Mode", fontsize=16, y=1.03)
    plt.savefig("FocalSV_combined_histogram_both.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Saved: FocalSV_combined_histogram_both.png")

# === Main ===
if __name__ == "__main__":
    target_data = load_phased_data(target_paths, mode="target")
    auto_data = load_phased_data(auto_paths, mode="auto")
    plot_combined_histogram(target_data, auto_data)
