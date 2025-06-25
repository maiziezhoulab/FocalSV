import os
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
parser = argparse.ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
# parser.add_argument('--bench_dir','-o')
args = parser.parse_args()
BASE_PATH = args.input_dir
# bench_dir = args.bench_dir

# === Config ===
# base_dir = "/data/maiziezhou_lab/Jamiezhou/region_based/Review/ps_hp"
libraries = ["hifi_l1", "clr_l1", "ont_l1"]

metric_keys = {
    "SER(%)": "SER(%)",
    "long_switch_error": "ref.het.phased.correct_call_in_pred.phased_in_pred.long_switch_error",
    "point_switch_error": "ref.het.phased.correct_call_in_pred.phased_in_pred.point_switch_error",
    "switch_error": "ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error",
    "correct_phased_ratio": ["ref.het.phased.correct_call_in_pred.phased_in_pred", "ref.het.phased"]
}

# Custom titles for plots
metric_titles = {
    "SER(%)": "Switch Error Rate (SER)",
    "long_switch_error": "Number of Long Switch Errors",
    "point_switch_error": "Number of Point Switch Errors",
    "switch_error": "Number of Total Switch Errors",
    "correct_phased_ratio": "Fraction of Correctly Phased Heterozygous Variants"
}

metric_names = list(metric_keys.keys())
all_metrics_by_library = {}  # Global structure for grid plot

bar_colors = {
    "hifi_l1": "#89A8B2",  # pastel blue
    "clr_l1": "#B3C8CF",   # pastel green
    "ont_l1": "#E5E1DA"    # pastel pink
}
mean_color = "#d62728"     # red
median_color = "#1f77b4"   # blue

# === Metric Parser ===
def extract_metrics(eval_file):
    data = {}
    try:
        with open(eval_file) as f:
            for line in f:
                if ":" not in line or line.startswith("#"):
                    continue
                key, value = line.strip().split(":", 1)
                key = key.strip().split()[0]
                value = value.strip().split("#")[0].strip()
                try:
                    data[key] = float(value)
                except ValueError:
                    continue

        # Compute correct_phased_ratio
        ratio_keys = metric_keys["correct_phased_ratio"]
        if all(k in data for k in ratio_keys):
            num = data[ratio_keys[0]]
            den = data[ratio_keys[1]]
            data["correct_phased_ratio"] = num / den if den > 0 else 0.0
        else:
            data["correct_phased_ratio"] = 0.0

        return data
    except Exception as e:
        print(f"❌ Error reading {eval_file}: {e}")
        return None

# === File Collector ===
def gather_eval_files(lib):
    eval_paths = []
    missing_regions = []
    total_regions = 0

    root = os.path.join(base_dir, lib)
    for chr_num in range(1, 23):
        region_dir = os.path.join(root, f"chr{chr_num}", "regions")
        if not os.path.exists(region_dir):
            continue
        for region in os.listdir(region_dir):
            total_regions += 1
            eval_path = os.path.join(region_dir, region, "phasing_eval.out")
            if os.path.isfile(eval_path):
                eval_paths.append(eval_path)
            else:
                missing_regions.append(os.path.join(region_dir, region))

    return eval_paths, total_regions, missing_regions

# === Parallel Extract ===
def parallel_extract(eval_paths, n_workers=8):
    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(extract_metrics, path): path for path in eval_paths}
        for future in as_completed(futures):
            res = future.result()
            if res:
                results.append(res)
    return results

# === Main Execution ===
if __name__ == "__main__":
    for lib in libraries:
        print(f"\n📦 Processing library: {lib}")
        lib_metrics = {key: [] for key in metric_keys}

        eval_files, total_regions, missing_regions = gather_eval_files(lib)

        print(f"📊 Total regions in {lib}: {total_regions}")
        print(f"✅ Regions with phasing_eval.out: {len(eval_files)}")
        print(f"⚠️ Missing phasing_eval.out in {len(missing_regions)} regions")

        if missing_regions:
            missing_file = f"{lib}_missing_regions.txt"
            with open(missing_file, "w") as f:
                for path in missing_regions:
                    f.write(path + "\n")
            print(f"📄 Missing region paths saved to: {missing_file}")

        parsed_results = parallel_extract(eval_files, n_workers=os.cpu_count() // 2 or 4)

        for metrics in parsed_results:
            for key in metric_keys:
                if key == "correct_phased_ratio":
                    lib_metrics[key].append(metrics.get("correct_phased_ratio", 0.0))
                else:
                    lib_metrics[key].append(metrics.get(metric_keys[key], 0.0))

        all_metrics_by_library[lib] = lib_metrics

    # === 3x5 Grid Plot with a/b/c in row corners and y-axis labels ===
    fig, axes = plt.subplots(3, 5, figsize=(24, 12), constrained_layout=True)
    row_labels = ["a", "b", "c"]
    y_axis_labels = ["HIFI_L1", "CLR_L1", "ONT_L1"]

    for row_idx, lib in enumerate(libraries):
        for col_idx, metric in enumerate(metric_names):
            ax = axes[row_idx, col_idx]
            values = all_metrics_by_library[lib][metric]

            # Convert ratio to percentage
            if metric == "correct_phased_ratio":
                values = [v * 100 for v in values]

            # Title for the column
            metric_display = metric_titles.get(metric, metric)

            mean_val = np.mean(values)
            median_val = np.median(values)
            bins = np.arange(0, max(values) + 2, 1) if "switch" in metric else 20

            # Histogram
            ax.hist(values, bins=bins, color=bar_colors[lib], edgecolor='black', alpha=0.9, linewidth=1.5)

            # Mean and median
            ax.axvline(mean_val, color=mean_color, linestyle='--', linewidth=2)
            ax.axvline(median_val, color=median_color, linestyle=':', linewidth=2)

            # Annotations
            ax.text(0.97, 0.93, f"Mean: {mean_val:.2f}", transform=ax.transAxes,
                    ha='right', va='top', fontsize=10, color=mean_color)
            ax.text(0.97, 0.85, f"Median: {median_val:.2f}", transform=ax.transAxes,
                    ha='right', va='top', fontsize=10, color=median_color)

            # Top row titles
            if row_idx == 0:
                ax.set_title(metric_display, fontsize=13)
            # Y-axis label for first column
            if col_idx == 0:
                ax.set_ylabel(y_axis_labels[row_idx], fontsize=13)

            ax.tick_params(axis='both', which='major', labelsize=10)
            ax.grid(axis='y', linestyle='--', alpha=0.4)

            if "switch" in metric:
                ax.set_yscale("log")
                ax.set_ylim(1, None)

    # Add a/b/c labels at the top-left of each row
    for row_idx, label in enumerate(["a", "b", "c"]):
        ax = axes[row_idx, 0]
        ax.annotate(label, xy=(-0.25, 1.00), xycoords='axes fraction',
                    fontsize=15, fontweight='bold', ha='center', va='top')

    # Save the plot
    final_colored_path = "phasing_metrics_grid_3x5_rowlabels_v2.png"
    fig.suptitle("Phasing Evaluation Metrics Across Sequencing Technologies", fontsize=16, y=1.03)
    plt.savefig(final_colored_path, dpi=300, bbox_inches='tight')
    plt.close()
