import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import pandas as pd
import json

color_map={
"FocalSV(target)":"magenta",
"FocalSV(auto)":"magenta",
"PAV":	"darkgreen",
"SVIM-asm":"springgreen",
"Dipcall":"skyblue",
"cuteSV":"royalblue",
"SVIM":"deeppink",
"PBSV":"coral",
"Sniffles2":"peru",
"SKSV":"green",
"sawfish":"darkorange",
}

asm=["FocalSV(target)", "sawfish", "PAV", "SVIM-asm", "Dipcall"]
aln=["cuteSV","SVIM","PBSV","Sniffles2","SKSV"]
hybrid=["FocalSV(auto)",]


def read_truvari_summary(summary_file):
    with open(summary_file,"r") as sumryf:
        summary_dict=json.load(sumryf)
    return summary_dict

def read_eval_results(tool_list, data_type):
    data_dict = dict()
    with open(tool_list, "r") as f:
        for line in f:
            tool, main_dir = line.rstrip("\n").split("\t")
            config_list = [i for i in os.listdir(main_dir) if os.path.isdir(main_dir+'/'+i)] #DEL_P0.3_50_
            data_list = []
            for i in config_list:
                d = read_truvari_summary(main_dir+"/"+i+"/summary.json")[data_type]
                if d is None:
                    d = 0
                data_list.append(d)
            eval_data = [i.split("_")[:2] + [d] for i, d in zip(config_list, data_list)] #[DEL, P0.3, d]
            eval_data = [(i[0]+"-"+i[1][0], float(i[1][1:]), i[2]) for i in eval_data] #(DEL-P, 0.3, d)
            eval_data.sort(key=lambda x: (x[0], x[1]))

            temp_dict = defaultdict(list)
            for para, para_value, d in eval_data:
                temp_dict[para].append((para_value, d))

            for para, values in temp_dict.items():
                d_list = [i[1] for i in values]
                if para not in data_dict:
                    data_dict[para] = pd.DataFrame({para.split("-")[-1]:[i[0] for i in values]})
                data_dict[para][tool] = d_list

    return data_dict

def plot_data(data_dict, save_dir, data_type):

    mpl.rcParams['pdf.fonttype'] = 42
    
    for para, df in data_dict.items():
        fig, ax = plt.subplots(figsize=(10, 6))
        asm_tools = [tool for tool in df.columns[1:] if tool in asm]
        aln_tools = [tool for tool in df.columns[1:] if tool in aln]
        hybrid_tools = [tool for tool in df.columns[1:] if tool in hybrid]

        df.plot(x=para.split("-")[-1],
                y=asm_tools,
                marker='o',
                ax=ax,
                ylim=(0,1.0),
                color=[color_map[tool] for tool in asm_tools])

        df.plot(x=para.split("-")[-1],
                y=aln_tools,
                marker='o',
                ax=ax,
                linestyle='--',
                ylim=(0,1.0),
                color=[color_map[tool] for tool in aln_tools])
        
        df.plot(x=para.split("-")[-1],
                y=hybrid_tools,
                marker='^',
                ax=ax,
                linestyle='-',
                ylim=(0,1.0),
                color=[color_map[tool] for tool in hybrid_tools])

        
        ax.set_title("Truvari "+data_type+" for "+para.split("-")[-1])
        ax.set_xlabel(para.split("-")[-1])
        ax.set_ylabel(data_type)
        ax.grid(True)
        plt.legend(loc=3, bbox_to_anchor=(1.05, 0), borderaxespad=0)

        plt.savefig(save_dir+"/Truvari_"+data_type+"_"+para+'.pdf',transparent=False, bbox_inches='tight')
        plt.close()

def main(tool_list, data_type, save_dir):

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    data_dict = read_eval_results(tool_list, data_type)
    plot_data(data_dict, save_dir, data_type)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--tool_list', type=str,)
    parser.add_argument('--data_type', type=str,default="f1")
    parser.add_argument('--save_dir', type=str,default=".")
    args = parser.parse_args()

    main(args.tool_list, args.data_type, args.save_dir)
