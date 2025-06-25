import os
import numpy as np
import pickle

# plot heatmap #
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

from collections import OrderedDict
import math

# def plot_heatmap(plot_info,save_dir,data_key,xtks,ytks):

#     for SV_dir in plot_info.keys():
#         fig, ax = plt.subplots(figsize=(10,8))
#         #sns.heatmap(plot_info[SV_dir], annot=False, ax=ax,cmap="viridis",vmin=0, vmax=1)
#         sns.heatmap(plot_info[SV_dir], annot=True, ax=ax,cmap="viridis",vmin=0, vmax=1,fmt='.3f')
#         ax.set_xlabel("r")
#         ax.set_ylabel("p",rotation=0)
#         ax.set_xticklabels(xtks)
#         ax.set_yticklabels(ytks)
#         ax.set_title(SV_dir+data_key)
#         plt.yticks(rotation=0)
#         mpl.rcParams['pdf.fonttype'] = 42
#         plt.savefig(save_dir+'/'+SV_dir+data_key+'_truvari_heatmap.pdf',bbox_inches='tight')
#         plt.close()
#     return 0

def read_truvari_summary(summary_file):
    with open(summary_file,"r") as sumryf:
        summary_dict=eval(sumryf.read().replace("null",'"NaN"'))
    return summary_dict

def heatmap_truvari_results_parser(truvari_results_dir,data_key="f1"): #data_key can be: "TP-call","FP","FN","recall","precision","f1","TP-call_TP-gt","TP-call_FP-gt","gt_recall","gt_precision","gt_f1"
    '''p_0.1_r_100/DEL_50_'''
    Truvari_config_list = [i for i in os.listdir(truvari_results_dir) if os.path.isdir(truvari_results_dir+'/'+i)]
    args = OrderedDict()

    for config in Truvari_config_list:
        config = config.split('_')
        if config[0] not in args:
            args[config[0]] = list()
        if config[2] not in args:
            args[config[2]] = list()
        args[config[0]].append(config[1])
        args[config[2]].append(config[3])


    #args[config[0]] = sorted(set(args[config[0]]),key=lambda x: float(x))
    #args[config[2]] = sorted(set(args[config[2]]),key=lambda x: float(x))
    try:
        for key, value in args.items():
            args[key] = sorted(set(value),key=lambda x: float(x))
    except:
        print(truvari_results_dir)

    heat_map_data = dict()
    
    for config_dir in Truvari_config_list:
        SV_dirs = [i for i in os.listdir(truvari_results_dir+'/'+config_dir) if os.path.isdir(truvari_results_dir+'/'+config_dir+'/'+i)]

        config = config_dir.split('_')
        arg_pair = (config[0],config[2])
        arg1_index = args[config[0]].index(config[1])
        arg2_index = args[config[2]].index(config[3])

        for SV_dir in SV_dirs:
            if SV_dir not in heat_map_data:
                heat_map_data[SV_dir] = dict()
            if arg_pair not in heat_map_data[SV_dir]:
                heat_map_data[SV_dir][arg_pair] = np.zeros((len(args[config[0]]),len(args[config[2]])))
                #heat_map_data[SV_dir] = np.zeros((len(args[config[0]]),len(args[config[2]])))

            summary_dict = read_truvari_summary(truvari_results_dir+'/'+config_dir+'/'+SV_dir+'/summary.json')
            if summary_dict[data_key] == 'NaN':
                summary_dict[data_key] = 0

            heat_map_data[SV_dir][arg_pair][arg1_index,arg2_index] = float(summary_dict[data_key])

    return args, heat_map_data
    #plot_heatmap(heat_map_data,truvari_results_dir,data_key,args[config[2]],args[config[0]])

def truvari_results_heatmap(input_dir, save_dir, data_key="f1", suffix = None):

    os.system("mkdir -p " + save_dir)
    #sv_tool_name\ttruvari_results_dir

    heat_map_info = dict()
    tool_num = 0

    tools_truvari_dict = {'_'.join(subdir.split('_')[1:3]): os.path.join(input_dir,subdir) for subdir in os.listdir(input_dir) if "Truvari" in subdir}
    tools_truvari_list = []
    # sort by lib order
    for dtype in ['Hifi','CLR','ONT']:
        for i in range(3):
            lib = f'{dtype}_{i+1}'
            if lib in tools_truvari_dict:
                tools_truvari_list.append((lib, tools_truvari_dict[lib]))
    # heatmap_data_
    # with open(tools_truvari_list,'r') as ttl:
    #     for line in ttl:
    #         if line[0]!='#':

    for line in tools_truvari_list:
        tool_num+=1
        line = line.rstrip('\n').split('\t')
        args, heat_map_data = heatmap_truvari_results_parser(line[1],data_key=data_key)
        for sv_type in heat_map_data.keys():
            if sv_type not in heat_map_info:
                heat_map_info[sv_type] = dict()
            for arg_pair in heat_map_data[sv_type].keys():
                if arg_pair not in heat_map_info[sv_type]:
                    heat_map_info[sv_type][arg_pair] = list()
                heat_map_info[sv_type][arg_pair].append([line[0], args, heat_map_data[sv_type][arg_pair]])

    with open(save_dir+'/heatmap_data.pkl','wb') as hmpk:
        pickle.dump(heat_map_info,hmpk)


    # print("2")
    with open(save_dir+'/heatmap_data.pkl', 'rb') as f:
        heat_map_info = pickle.load(f)
    # print("2")


    col_hm = max(math.ceil(tool_num**0.5) ,3)

    col = col_hm+1 #additional axis for colorbar
    row = math.ceil(tool_num/col_hm)

    # --- Customizable parameters ---
    annot_fontsize = 15
    tick_fontsize = 21
    label_fontsize = 28
    cbar_labelsize = 28
    fig_labelsize = 45
    # -------------------------------
    
    for sv_type in heat_map_info.keys():
        for arg_pair in heat_map_info[sv_type].keys():
            if sv_type == 'INS_50_':
                if arg_pair!=('p','r'):
                    continue
            else:
                if arg_pair!=('p','O'):
                    continue
            

            fig, axes = plt.subplots(
                row, col,
                figsize=(col_hm * 10 + 3, row * 10),
                gridspec_kw={'width_ratios': [10] * col_hm + [0.5]},
                squeeze=False
            )

            mpl.rcParams['pdf.fonttype'] = 42

            for i, hm_info in enumerate(heat_map_info[sv_type][arg_pair]):
                ax_row = i // col_hm
                ax_col = i % col_hm
                heatmap_data = hm_info[-1]
                xticks = hm_info[1][arg_pair[1]]
                yticks = hm_info[1][arg_pair[0]]
                title = hm_info[0]

                # Choose cbar axis only for the last column
                if ax_col == col_hm - 1:
                    sns.heatmap(
                        heatmap_data, annot=True, fmt='.3f',
                        ax=axes[ax_row][ax_col],
                        cmap="viridis", vmin=0, vmax=1,
                        cbar=True, cbar_ax=axes[ax_row][ax_col + 1],
                        annot_kws={"size": annot_fontsize}
                    )
                    axes[ax_row][ax_col + 1].tick_params(labelsize=cbar_labelsize)
                else:
                    sns.heatmap(
                        heatmap_data, annot=True, fmt='.3f',
                        ax=axes[ax_row][ax_col],
                        cmap="viridis", vmin=0, vmax=1,
                        cbar=False, annot_kws={"size": annot_fontsize}
                    )

                # Add subplot labels like (a), (b), ...
                subplot_label = f"{chr(97 + i)}"  # 97 is ASCII for 'a'
                axes[ax_row][ax_col].text(
                    -0.1, 1.08, subplot_label,
                    transform=axes[ax_row][ax_col].transAxes,
                    fontsize=fig_labelsize, fontweight='regular',
                    va='top', ha='left'
                )

                # Customize ticks and labels
                axes[ax_row][ax_col].set_xticklabels(xticks, fontsize=tick_fontsize)
                axes[ax_row][ax_col].set_yticklabels(yticks, fontsize=tick_fontsize, rotation=0)
                axes[ax_row][ax_col].set_xlabel(arg_pair[1], fontsize=label_fontsize)
                axes[ax_row][ax_col].set_ylabel(arg_pair[0], fontsize=label_fontsize, rotation=0, labelpad = 20)
                axes[ax_row][ax_col].set_title(title, fontsize=label_fontsize)

            # axes = axes.flatten()
            # for i in range(len(heat_map_info[sv_type][arg_pair]), row * col):
            #     axes[i].axis("off")
            # plt.suptitle(f"{sv_type} {data_key} {translate(arg_pair)}", y=0.95, fontsize=label_fontsize + 2)
            if suffix is None:
                outfile = f"{save_dir}/{sv_type}{data_key}_{translate(arg_pair)}_truvari_heatmap.pdf"
            else:
                outfile = f"{save_dir}/{sv_type}{data_key}_{translate(arg_pair)}_truvari_heatmap_{suffix}.pdf"
            plt.savefig(outfile, bbox_inches='tight')
            plt.close()


            print(outfile)
            # exit()


def translate(arg_pair):
    dc = {"p": "pseq",
          "P":"Psize",
          'r':'r',
          'O':'O'}
    return dc[arg_pair[0]]+'_'+dc[arg_pair[1]]
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_dir', type=str,)
    parser.add_argument('--data_type', type=str,default="f1")
    parser.add_argument('--save_dir', type=str,default=".")
    parser.add_argument('--suffix', type=str,default=None)
    args = parser.parse_args()

    # tools_truvari_list = args.tools_truvari_list
    input_dir = args.input_dir
    data_type = args.data_type
    save_dir = args.save_dir
    suffix = args.suffix

    truvari_results_heatmap(input_dir, save_dir, data_key=data_type, suffix = suffix)

