from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dirs','-i', nargs='+')
parser.add_argument('--output_dir','-o')
parser.add_argument('--libnames','-libs', nargs= '+')
parser.add_argument('--color_dict_path','-cl', help = "json file")
parser.add_argument('--target_tool','-tg', help = "target tool")
# parser.add_argument('--target_tool_alias','-tga', help = "target tool alias")
parser.add_argument('--tool_order','-od',help = "separated by comma; optional, if None, order by alphabet ")
#parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_dirs = args.input_dirs
libnames = args.libnames
out_dir = args.output_dir
tool_order = args.tool_order
target_tool = args.target_tool
# target_tool_alias = args.target_tool_alias
color_dict_path = args.color_dict_path
#output_path = args.output_path
if tool_order is not None:
    tool_order = tool_order.split(',')[::-1]



import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import os
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Create figure and axis #1
# colors=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#800000', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
# tls=['focalsv','nanosv','smartie-sv_aln','sniffles','svim','cutesv','nanovar','pbsv','sksv','sniffles2','mamnet','debreak','dipcall', 'pav','smartie-sv_asm','svim-asm']

# color_dict = dict()
# for i in range(len(tls)):
#     color_dict[tls[i]]=colors[i]
with open(color_dict_path, 'r') as f:
    color_dict = eval(f.read())


import os
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


import os
import math
import string
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def draw_plots_grid(work_dirs, out_dir, y_column, libnames,
                    tool_order=None, target_tool=None, 
                    tool_dict=None, color_dict=None,
                    legend=True):
    """
    work_dirs: list of up to 9 directories, each containing tools.csv
    out_dir: where to write the combined PDF
    y_column: e.g. "F1"
    libnames: list of titles for each subplot
    tool_order: optional list to fix the order of tools
    target_tool: name of the tool you want to highlight
    tool_dict: mapping tool->'asm' or other, to pick linestyles
    color_dict: mapping tool->color
    """
    os.makedirs(out_dir, exist_ok=True)
    n = len(work_dirs)
    cols = 3
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(28, 15), sharex=True)
    axes = axes.flatten()

    # prepare subplot labels 'a'.. 'i'
    labels = list(string.ascii_lowercase[:n])

    line_handles = None
    line_labels = None

    for i, work_dir in enumerate(work_dirs):
        ax1 = axes[i]
        ax2 = ax1.twinx()

        # --- load and prep data ---
        df = pd.read_csv(os.path.join(work_dir, "tools.csv"))
        tool_set = {x.split("_")[-1] for x in df.columns[2:-2]}
        if tool_order is not None:
            tool_list = [t for t in tool_order if t in tool_set]
        else:
            tool_list = sorted(tool_set)

        # bar‐color by deletion/insertion
        df["color"] = df["lenth_range"].apply(
            lambda x: "lightgray" if eval(x)[0] < 0 else "darkgray"
        )
        df["lenth_range"] = df["lenth_range"].str[1:-1]

        # --- plot bars ---
        ax1.bar(df["lenth_range"], df["SV_number"],
                width=0.5, alpha=0.7, color=df["color"])
        ax1.set_ylabel("Number of SVs", size=14)
        if i > 5:
            ax1.set_xlabel("SV size", size=14)
        ax1.tick_params(axis="x", rotation=90, labelsize=12)
        ax1.tick_params(axis="y", labelsize=12)

        # --- add bolded subplot label ---
        ax1.text(
            -0.15, 1.08,
            labels[i],
            transform=ax1.transAxes,
            fontsize=45, fontweight="regular",
            va="top", ha="left"
        )

        # --- plot lines ---
        for tool in tool_list:
            ls = "solid" if tool_dict.get(tool) == "asm" else "dashed"
            lw = 2.0 if tool == target_tool else 1.6
            alpha = 1.0 if tool == target_tool else 0.9

            ax2.plot(df["lenth_range"], df[f"{y_column}_{tool}"],
                     label=tool, linestyle=ls,
                     linewidth=lw, alpha=alpha,
                     color=color_dict[tool])

        ax2.set_ylabel(y_column, size=14)
        ax2.tick_params(axis="y", labelsize=12)

        # title = user‐provided libname
        ax1.set_title(libnames[i], size=20, pad = 10)

        # capture line‐plot legend from the first axes
        if legend and line_handles is None:
            line_handles, line_labels = ax2.get_legend_handles_labels()

    # remove any unused subplots
    for j in range(n, rows * cols):
        fig.delaxes(axes[j])

    if legend:
        # 1) bar legend (deletion/insertion) above
        bar_handles = [
            Patch(facecolor="lightgray", edgecolor=None, alpha=0.7, label="Deletion"),
            Patch(facecolor="darkgray",  edgecolor=None, alpha=0.7, label="Insertion"),
        ]
        fig.legend(bar_handles, ["Deletion", "Insertion"],
                   loc="upper left",
                   bbox_to_anchor=(0.9, 0.98),
                   frameon=False,
                   prop={"size": 15})

        # 2) line legend (tools) below
        if line_handles:
            fig.legend(line_handles, line_labels,
                       loc="upper left",
                       bbox_to_anchor=(0.9, 0.94),
                       frameon=False,
                       prop={"size": 15})

    plt.tight_layout(rect=[0, 0, 0.90, 1.00])
    outfile = os.path.join(out_dir, "3x3_diffsvlen_f1.pdf")
    fig.savefig(outfile, bbox_inches="tight")
    plt.show()

tool_list = ["PBHoney",
"NanoSV",
"Smartie-sv_aln",
"Sniffles",
"SVIM",
"NanoVar",
"PBSV",
"cuteSV",
"SKSV",
"DeBreak",
"MAMnet",
"Sniffles2",
"Dipcall",
"Smartie-sv_asm",
"PAV",
"SVIM-asm",
"FocalSV(target)",
"FocalSV(auto)",
"sawfish"
]
tool_type = ['aln'] * 12 + ['asm']*7
tool_dict = dict(zip(tool_list,tool_type))
# df=pd.read_csv('line_barplot/p0.1_O0.1_r900/'+"/tools.csv")
#data_dir = "config/"
#out_dir = 'figure/'
# draw_plot('F1',input_dir,input_dir,legend=True)

draw_plots_grid(input_dirs, out_dir, 'F1',libnames,
                    tool_order, target_tool,
                    tool_dict, color_dict,
                    legend=True)
#draw_plot('F1',data_dir + '/p0.3_O0.3_r600/',out_dir,legend=True)
#draw_plot('F1',data_dir + '/p0.6_O0.6_r300/',out_dir,legend=True)
#draw_plot('F1',data_dir + '/p0.9_O0.9_r100/',out_dir,legend=True)
