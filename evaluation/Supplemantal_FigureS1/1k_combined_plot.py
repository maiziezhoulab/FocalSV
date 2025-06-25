import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42  # TrueType fonts (editable)
mpl.rcParams['ps.fonttype'] = 42


# === Paths and constants ===
csv_dir = "./"
bed_raw_path = os.path.join(csv_dir, "1k_SV_negative_regions.bed")
bed_norm_path = os.path.join(csv_dir, "1k_SV_negative_regions.bed")
# f1_path = os.path.join(csv_dir, "Recall_F1.csv")
f1_path = os.path.join("./", "Recall_F1.csv")



output_dirs = {
    "clr": "clr_all_tools",
    "hifi": "hifi_all_tools",
    "ont": "ont_all_tools"
}

bar_colors = {
    "hifi": "#89A8B2",
    "clr": "#B3C8CF",
    "ont": "#E5E1DA"
}

tools = ["FocalSV", "PAV", "SVIM-asm", "Dipcall",'sawfish', "cuteSV", "SVIM", "PBSV", "Sniffles2", "SKSV"]
chrom_order = [f"chr{i}" for i in range(1, 23)]
title_map = {"hifi": "HIFI_L1", "clr": "CLR_L1", "ont": "ONT_L1"}

# === Load BEDs ===
bed_raw_df = pd.read_csv(bed_raw_path, sep="\t", header=None, names=["chr", "start", "end"])
bed_raw_regions = set(zip(bed_raw_df["chr"], bed_raw_df["start"], bed_raw_df["end"]))

bed_norm_df = pd.read_csv(bed_norm_path, sep="\t", header=None, names=["chr", "start", "end"])
bed_norm_regions = set(zip(bed_norm_df["chr"], bed_norm_df["start"], bed_norm_df["end"]))

# === Load F1 Table ===
f1_data_raw = pd.read_csv(f1_path, header=0)
f1_data_raw = f1_data_raw.rename(columns={f1_data_raw.columns[1]: "Metric"})
f1_data = f1_data_raw.set_index(f1_data_raw.columns[0])

# === Helper: F1 averages ===
def get_avg_f1_from_consolidated(data_type):
    del_row = f1_data.loc[f"{data_type.upper()} DEL"]
    ins_row = f1_data.loc[f"{data_type.upper()} INS"] if f"{data_type.upper()} INS" in f1_data.index else None
    del_f1 = del_row[del_row["Metric"] == "F1"]
    ins_f1 = ins_row[ins_row["Metric"] == "F1"] if ins_row is not None else None

    avg_f1 = {}
    for tool in tools:
        try:
            del_val = float(str(del_f1[tool].values[0]).replace(",", "")) if tool in del_f1 else 1.0
            ins_val = float(str(ins_f1[tool].values[0]).replace(",", "")) if ins_f1 is not None and tool in ins_f1 else del_val
            avg_f1[tool] = (del_val + ins_val) / 2
        except:
            avg_f1[tool] = 1.0
    return avg_f1

# === Count total FP over BED regions ===
def count_fp(data_type, bed_regions, normalize=False):
    raw_fp = {}
    f1_scores = get_avg_f1_from_consolidated(data_type) if normalize else None

    for tool in tools:
        # csv_path = os.path.join(csv_dir,output_dirs[data_type], f"sv_counts_by_region_{tool}.csv")
        csv_path = os.path.join("./",output_dirs[data_type], f"sv_counts_by_region_{tool}.csv")
        if not os.path.exists(csv_path):
            print(f"[Missing] {csv_path}")
            continue

        df = pd.read_csv(csv_path)
        df = df[df["chr"].isin(chrom_order)]
        df_bed = df[df.apply(lambda row: (row["chr"], row["start"], row["end"]) in bed_regions, axis=1)]
        total_fp = df_bed["sv_count"].sum()

        if normalize:
            adjusted_fp = total_fp * (1 - f1_scores.get(tool, 1.0))
            raw_fp[tool] = adjusted_fp
        else:
            raw_fp[tool] = total_fp

    return pd.Series(raw_fp).sort_values(ascending=False)

# === Plot all in one figure ===
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11, 11))

column_titles = ["Raw False Positive Counts", "Normalized False Positive Counts"]
row_labels = ['a', 'b', 'c']
row_names = ["HIFI_L1", "CLR_L1", "ONT_L1"]
bar_width = 0.8  # full-width bars now

for row_idx, dt in enumerate(["hifi", "clr", "ont"]):
    raw = count_fp(dt, bed_raw_regions, normalize=False)
    norm = count_fp(dt, bed_norm_regions, normalize=True)

    for col_idx, (data, ax) in enumerate(zip([raw, norm], axes[row_idx])):
        # Plot directly with tool names (categorical axis)
        bars = ax.bar(
            data.index,
            data.values,
            width=bar_width,
            color=bar_colors[dt],
            edgecolor='black'
        )

        y_max = max(data.values) * 1.15
        ax.set_ylim(0, y_max)

        for bar in bars:
            height = bar.get_height()
            label = f"{height:.2f}" if col_idx == 1 else f"{int(height)}"
            ax.annotate(label,
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=10)

        # ax.set_xticklabels(data.index, rotation=45, ha='right', fontsize=8)
        ax.tick_params(axis='x', labelrotation=45, labelsize=10)
        ax.set_yticks([])
        ax.set_ylabel("")
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if col_idx == 0:
            ax.annotate(row_names[row_idx], xy=(-0.02, 0.5), xycoords='axes fraction',
                        fontsize=11, fontweight='bold', ha='right', va='center', rotation=90)
            ax.annotate(row_labels[row_idx], xy=(-0.09, 1.0), xycoords='axes fraction',
                        fontsize=13, fontweight='bold', ha='center', va='top')

# === Manual top titles for compact layout ===
# === Adjust layout to leave space for top titles ===
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave top 5% space

# === Add column titles within the canvas ===
fig.text(0.27, 0.96, "Raw False Positive Counts", ha='center', fontsize=14, fontweight='bold')
fig.text(0.75, 0.96, "Normalized False Positive Counts", ha='center', fontsize=14, fontweight='bold')
fig.suptitle("False Positive Structural Variant Counts", fontsize=15, weight='bold', y=1.02)


# === Save final plot ===
# plt.savefig("final_combined_fp_plot_ultracompact.png", dpi=300, bbox_inches='tight')
plt.savefig("final_combined_fp_plot_ultracompact.pdf",  bbox_inches='tight')
plt.close()
print("✅ Final compact plot saved as PDF: final_combined_fp_plot_ultracompact.pdf")