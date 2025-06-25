import re
import matplotlib
import matplotlib.pyplot as plt
import gzip
import math
import os
matplotlib.use('agg')

def plot_breakpoint_shift_distribution(breakpoint_shifts, save_dir_list, sv_type):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    num_dirs = len(breakpoint_shifts)
    rows = cols = math.ceil(math.sqrt(num_dirs))

    fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 8 * rows))
    axes = axes.flatten()  # Flatten the axes array for easier iteration
    fontsize = 24
    titlesize = 28

    for idx, (shift, save_dir) in enumerate(zip(breakpoint_shifts, save_dir_list)):
        if shift[sv_type]:
            color = 'blue' if sv_type == 'DEL' else 'red'
            axes[idx].hist(shift[sv_type], bins=200, range=None, histtype='bar', log=False, 
                           color=color, label=sv_type)
            axes[idx].set_xlim(-210, 210)

            # Set x and y labels
            axes[idx].set_xlabel("Breakpoint Shift", fontsize=fontsize)
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
                       fontsize=titlesize + 6, fontweight='regular', va='top', ha='right')

        axes[idx].set_title(f"{libs[idx]}", fontsize=titlesize)

    # Hide any unused axes
    for ax in axes[len(breakpoint_shifts):]:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'{out_dir}/breakpoint_shift_distribution_all_{sv_type}.pdf', bbox_inches='tight')
    plt.close(fig)


def group_by_matchid(vcf):
    brkpt_dict = dict()

    with gzip.open(vcf, 'rt') as vf:
        for line in vf:
            if line[0] != '#':
                matchid = re.findall("MatchId=(\d+)", line)[0]
                fields = line.rstrip('\n').split('\t')

                try:
                    svtype = re.findall("SVTYPE=(\w+)", line)[0]
                    svlen = float(re.findall("SVLEN=-?(\d+)", line)[0])
                except:
                    svlen = abs(len(fields[3]) - len(fields[4]))
                    svtype = 'DEL' if len(fields[3]) > len(fields[4]) else 'INS'

                start = float(fields[1])
                end = start + svlen if svtype == 'DEL' else start

                brkpt_dict[matchid] = [svtype, start, end]

    return brkpt_dict

def breakpoint_shift_analysis(truv_rslt_dir, svtype):
    tp_base_vcf = f"{truv_rslt_dir}/{svtype}_50_/tp-base.vcf.gz"
    tp_call_vcf = f"{truv_rslt_dir}/{svtype}_50_/tp-comp.vcf.gz"

    base_brkpt = group_by_matchid(tp_base_vcf)
    call_brkpt = group_by_matchid(tp_call_vcf)

    breakpoint_shift = {'DEL': [], 'INS': []}

    for matchid in base_brkpt.keys():
        if matchid not in call_brkpt:
            continue

        base_bp = base_brkpt[matchid]
        call_bp = call_brkpt[matchid]

        # if base_bp[0] != 'DEL':
        #     continue

        bp_shift = max(abs(call_bp[1] - base_bp[1]), abs(call_bp[2] - base_bp[2]))
        bp_shift = max(-201, min(bp_shift, 201))  # Clip to [-201, 201]

        breakpoint_shift[svtype].append(bp_shift)

    return breakpoint_shift

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

    for sv_type in ['DEL', 'INS']:
        print(f"Processing {sv_type}...")

        all_shifts = []
        for dir_path in truv_rslt_dirs:
            shifts = breakpoint_shift_analysis(dir_path, sv_type)
            all_shifts.append(shifts)

        plot_breakpoint_shift_distribution(all_shifts, truv_rslt_dirs, sv_type)

