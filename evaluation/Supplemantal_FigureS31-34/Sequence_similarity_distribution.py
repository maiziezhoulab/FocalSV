import re
import matplotlib
import matplotlib.pyplot as plt
import gzip
import math
import os
matplotlib.use('agg')

def plot_sequence_similarity_distribution(seqsim_list, save_dir_list,sv_type):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    num_dirs = len(seqsim_list)
    rows = cols = math.ceil(math.sqrt(num_dirs))

    fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 8 * rows))
    axes = axes.flatten()  # Flatten for easier iteration
    fontsize = 24
    titlesize = 28

    for idx, (seqsim, save_dir) in enumerate(zip(seqsim_list, save_dir_list)):
        if seqsim[sv_type]:
            if sv_type == 'INS':
                color = 'red'
            else:
                color = 'blue'
            axes[idx].hist(seqsim[sv_type], bins=200, range=None, histtype='bar', log=False, 
                        color= color, label= sv_type )
            axes[idx].set_xlim(0.5, 1.0)

            # Add axis labels
            axes[idx].set_xlabel("Sequence Similarity", fontsize=fontsize)
            axes[idx].set_ylabel("SV Count", fontsize=fontsize)

            # Adjust tick size
            axes[idx].tick_params(axis='both', which='major', labelsize=fontsize)

            # Adjust legend size and position
            axes[idx].legend(prop={'size': fontsize}, loc='upper right')

        else:
            axes[idx].text(0.5, 0.5, f'No {sv_type} data', ha='center', va='center', fontsize=12)

        # Add subplot labels (e.g., 'a', 'b', 'c', etc.)
        label = chr(97 + idx)  # Convert index to letter (97 is 'a')
        axes[idx].text(-0.05, 1.07, label, transform=axes[idx].transAxes, 
                    fontsize=titlesize+6, fontweight='regular', va='top', ha='right')

        axes[idx].set_title(f"{libs[idx]}", fontsize=titlesize)

    # Hide any unused axes
    for ax in axes[len(seqsim_list):]:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'{out_dir}/sequence_similarity_distribution_all_{sv_type}.pdf', bbox_inches='tight')
    plt.close(fig)


def sequence_similarity_analysis(truv_rslt_dir, svtype):
    tp_call_vcf = f"{truv_rslt_dir}/{svtype}_50_/tp-comp.vcf.gz"

    seqsim_dict = {'INS': [], 'DEL':[]}

    with gzip.open(tp_call_vcf, 'rt') as vf:
        for line in vf:
            if line[0] != '#':
                seqsim = float(re.findall(r"PctSeqSimilarity=(\d*\.?\d*)", line)[0])
                fields = line.rstrip('\n').split('\t')

                try:
                    svtype = re.findall(r"SVTYPE=(\w+)", line)[0]
                except:
                    svtype = 'INS' if len(fields[3]) < len(fields[4]) else 'DEL'

                if svtype == 'INS':
                    seqsim_dict['INS'].append(seqsim)
                else:
                    seqsim_dict['DEL'].append(seqsim)

                

    return seqsim_dict


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--truv_rslt_dirs', '-trds', nargs='+', help="List of directories containing Truvari results")
    parser.add_argument('--libs', '-libs', nargs='+', help="List of library names")
    parser.add_argument('--out_dir', '-o', help="Output directory to save plots", required=True)

    args = parser.parse_args()

    truv_rslt_dirs = args.truv_rslt_dirs
    libs = args.libs
    out_dir = args.out_dir

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    assert len(truv_rslt_dirs) == len(libs)

    for svtype in ['INS','DEL']:
        print(f"Processing {svtype}...")
        all_seqsim = []
        for dir_path in truv_rslt_dirs:
            seqsim = sequence_similarity_analysis(dir_path,svtype)
            all_seqsim.append(seqsim)

        plot_sequence_similarity_distribution(all_seqsim, truv_rslt_dirs,svtype)

