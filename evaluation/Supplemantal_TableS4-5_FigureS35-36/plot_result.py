import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import glob
from collections import defaultdict

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 # TrueType fonts (editable)
mpl.rcParams['ps.fonttype'] = 42


def organize_df(eval_params, dfs, ):
    use_metrices = ['bench','Recall','Precision','F1']
    df = pd.concat(dfs,axis = 0)
    filtered_df = df[df['metric'].isin(use_metrices)].reset_index(drop = True)
    n_list = []
    for n in eval_params:
        n_list.extend([n]*len(use_metrices))
    
    filtered_df.insert(0,'number of supporting tools',n_list)
    return filtered_df
    
def write_excel(dc_stats, eval_params, outfile):
    with pd.ExcelWriter(outfile) as writer:
        for var, dfs in dc_stats.items():
            df = organize_df(eval_params, dfs)
            df.to_excel(writer, sheet_name=str(var), index=False)



def plot_eval_results(excel_files, eval_params, outdir,sample):
    """
    For each variant type (each sheet), plot a figure with 4 subplots (Recall, Precision, F1, and bench).
    Each subplot will contain line plots of evaluation results vs. evaluation parameter for each tool.
    
    Parameters:
      excel_files (list of str): List of Excel file paths. Each file corresponds to one evaluation parameter.
      eval_params (list): List of evaluation parameter values corresponding to each Excel file.
      outdir (str): Output directory for saving the figures.
    """
    # Define the metrics we are interested in
    metrics_to_plot = ['Recall', 'Precision', 'F1', 'bench']
    
    # Dictionary to hold the data for each variant type. 
    # Structure: {variant_type: {metric: {tool: [list of eval results across files]}}}
    variant_results = {}

    dc_stats = defaultdict(list)
    
    # Process each Excel file and its corresponding evaluation parameter.
    for idx, file in enumerate(excel_files):
        eval_value = eval_params[idx]
        # Read all sheets from the Excel file; each sheet corresponds to one variant type.
        sheets_dict = pd.read_excel(file, sheet_name=None)
        for variant, df in sheets_dict.items():
            tools =['FocalSV','PAV','SVIM-asm','Dipcall','sawfish','cuteSV','SVIM','PBSV','Sniffles2','SKSV']
            dc_stats[variant].append(df[['metric']+tools])
            # Initialize the dictionary structure if this is the first time encountering the variant.
            if variant not in variant_results:
                variant_results[variant] = {metric: {} for metric in metrics_to_plot}
            
            # Ensure that the first column is treated as a string for metric comparison.
            # Assume the first column contains the metric names.
            metric_names = df.iloc[:, 0].astype(str)
            
            # Process each desired metric.
            for metric in metrics_to_plot:
                # Find the row where the first column matches the metric (case-insensitive).
                mask = metric_names.str.lower() == metric.lower()
                row = df[mask]
                if not row.empty:
                    # We take the first matching row.
                    row = row.iloc[0]
                    # For each tool (columns other than the first), record the result.
                    for tool in df.columns[1:]:
                        result_value = row[tool]
                        # Initialize the list for this tool if necessary.
                        if tool not in variant_results[variant][metric]:
                            variant_results[variant][metric][tool] = []
                        variant_results[variant][metric][tool].append(result_value)
                else:
                    print(f"Warning: Metric '{metric}' not found in sheet '{variant}' in file '{file}'")
    
    # Ensure the output directory exists.
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # write excel
    outfile = outdir+"/All_tool_eval_summary.xlsx"
    write_excel(dc_stats, eval_params, outfile)

    # Define the mapping of your data labels to line styles.
    line_style_dict = {
        'Hifi': '-',
        'CLR': '--',
        'ONT': ':'  # Change 'Other' to the appropriate label if needed.
    }
    # Define a color dictionary for 5 tools


    color_dict ={

    "FocalSV":"magenta",
    "FocalSV-auto":"magenta",
    "PAV":  "darkgreen",
    "SVIM-asm":"springgreen",
    "Dipcall":"skyblue",
    "cuteSV":"royalblue",
    "SVIM":"deeppink",
    "pbsv":"coral",
    "PBSV":"coral",
    "Sniffles2":"peru",
    "SKSV":"yellow",
    "sawfish": "red"
    }


    # For each variant type, create a figure with 4 subplots.
    variants = list(variant_results.keys())[::-1]
    # print(variants)
    # for variant, metrics_data in variant_results.items():

    subplot_labels = ['a', 'b', 'c', 'd']

    for variant in variants:
        metrics_data = variant_results[variant]
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        axs = axs.flatten()  # Flatten so we can iterate over them
        
        for i, metric in enumerate(metrics_to_plot):
            ax = axs[i]
            # Plot a line for each tool for the current metric.
            tools = list(metrics_data[metric].keys())[::-1]
            tools = ['SKSV','Sniffles2','PBSV','SVIM','cuteSV','sawfish','Dipcall','SVIM-asm','PAV','FocalSV']
            # for tool, results in metrics_data[metric].items():
            for tool in tools:
                results = metrics_data[metric][tool]
                # dat, tool = dat_tool.split('@')
                dat = "Hifi"
                if tool in ['FocalSV','Dipcall','PAV','sawfish','SVIM-asm']:
                    linestyle = 'solid'
                else:
                    linestyle = 'dashed'


                if metric == 'bench':
                    line_color = 'darkgray'
                else:
                    line_color = color_dict[tool]
                ax.plot(eval_params, results, marker='o', label = tool.replace("FocalSV",'FocalSV(auto)'), color = line_color, linewidth = 1, linestyle = linestyle)
                # Add subplot label (e.g., (a), (b), ...)
                ax.text(-0.1, 1.05, subplot_labels[i], transform=ax.transAxes,
                        fontsize=25, fontweight='normal', va='top', ha='right')

                            
            # ax.set_title(metric.replace('bench','Number of bench'))
            ax.set_xlabel('Number of supporting tools')
            ax.set_ylabel(metric.replace('bench','Number of bench'))
            # Removed individual legends from each subplot

        # Create one common legend for all subplots and place it in the top right corner outside the axes.
        handles_labels = {}
        for ax in axs:
            handles, labels = ax.get_legend_handles_labels()
            for handle, label in zip(handles, labels):
                if label not in handles_labels:
                    handles_labels[label] = handle
        if handles_labels:
            fig.legend(list(handles_labels.values()), list(handles_labels.keys()), 
                       loc='upper right', bbox_to_anchor=(1.15, 0.95))
            
        from matplotlib.lines import Line2D

        
        # # Create custom legend handles for line styles.
        # line_style_handles = [Line2D([0], [0], linestyle=ls, color='black', label=label)
        #                     for label, ls in line_style_dict.items()]

        # # Add the custom legend for line styles.
        # # Adjust loc and bbox_to_anchor as necessary to avoid overlapping with other legends.
        # fig.legend(handles=line_style_handles, title="Line Style", loc='upper right', bbox_to_anchor=(1.08, 0.75))

        


        if suffix:
            title = f'{sample} Variant: {variant}  ({suffix})'
            filename = f"{sample}_{variant}_{suffix}"
        else:
            title = f'{sample} {variant}'
            filename = f"{sample}_{variant}"

        fig.suptitle(title, fontsize=16)
        # Adjust layout to prevent overlapping titles/subplots
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.2, wspace=0.2)
        outfile = os.path.join(outdir, f"{filename}.pdf")
        fig.savefig(outfile, bbox_inches = 'tight')
        outfile = os.path.join(outdir, f"{filename}.png")
        fig.savefig(outfile, bbox_inches = 'tight')
        plt.close(fig)

    
    print("Plotting complete. Figures saved in", outdir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Draw line plots of eval_para vs eval_result for Recall, Precision, F1 and bench.")

    parser.add_argument('--work_dir', '-w',required=True, help='input dir and output directory to save figures')
    parser.add_argument('--suffix', '-suffix', default='')
    args = parser.parse_args()


    work_dir = args.work_dir
    suffix = args.suffix
    excels = glob.glob(work_dir+"/All*union*xlsx")
    patient = excels[0].split('_')[-3]
    source_list = set([ excel.split('_')[-2] for excel in excels])
    eval_params = sorted( list(set([ int(excel.split('_')[-1].split('.')[0]) for excel in excels])))
    
    for source in source_list:
        excels = [ f"{work_dir}/All_tool_eval_union_{patient}_{source}_{n}.xlsx" for n in eval_params ]
        sample = f"{patient}"
        print(sample)
        plot_eval_results(excels, eval_params, work_dir, sample)
        
        # break
