import pickle
import pysam
from collections import defaultdict

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def parse_single_cell_bam(cell_bam, cell_strand, min_mapq=20):

    strand_count = defaultdict(lambda : {"P":0, "M":0}) # NOTE: pickle cannot handle defualtdict with lambda
    if cell_strand == "P":
        cell_strand_r = "M"
    elif cell_strand == "M":
        cell_strand_r = "P"
    else:
        raise ValueError("Invalid cell strand label:", cell_strand, "should be either P or M")
    
    bam = pysam.AlignmentFile(cell_bam, "rb")
    for contig in bam.references:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if not read.is_read1:
                continue  # Strand-seq: only use R1

            if read.is_reverse:
                strand_count[contig][cell_strand_r] += 1
            else:
                strand_count[contig][cell_strand] += 1

    return strand_count

def calculate_phasing_quality(informative_cell_pkl, input_dir ,chroms=None, min_mapq=20):

    if chroms is None:
        chroms = ["chr"+str(i) for i in range(1, 23)] # default to all autosomes

    with open(informative_cell_pkl, "rb") as f:
        informative_cell = pickle.load(f)

    for chrom in chroms:
        chrom_strand_count = defaultdict(lambda : {"P":0, "M":0})
        for cell, cell_strand in informative_cell[chrom]:

            cell_chrom_strand_count = parse_single_cell_bam(f"{input_dir}/{chrom}/{cell}.bam", cell_strand, min_mapq)
            for contig, s_count in cell_chrom_strand_count.items():
                for s, c in s_count.items():
                    chrom_strand_count[contig][s] += c

        PhaseQ_dict = dict()

        phase_blocks = set([contig.replace("HP1_", "").replace("HP2_", "") for contig in chrom_strand_count.keys()])
        
        for pb in phase_blocks:
            valid_pb = True
            if "HP1_"+pb not in chrom_strand_count:
                print("HP1_"+pb, "has 0 count for both P and M reads, probably not covered by any valid strand-seq reads")
                valid_pb = False
            if "HP2_"+pb not in chrom_strand_count:
                print("HP2_"+pb, "has 0 count for both P and M reads, probably not covered by any valid strand-seq reads")
                valid_pb = False
            if "unphased" in pb: # exclude unphased blocks
                valid_pb = False

            if not valid_pb:
                continue

            # if "HP1_"+pb not in chrom_strand_count:
            #     print("HP1_"+pb, "has 0 count for both P and M reads, probably not covered by any valid strand-seq reads")
            #     print("hap score default to 0")
            #     hp1_score = 0
            # else:
            hp1_p = chrom_strand_count["HP1_"+pb]["P"]
            hp1_m = chrom_strand_count["HP1_"+pb]["M"]
            hp1_score = (hp1_p-hp1_m)/(hp1_p+hp1_m)

            # if "HP2_"+pb not in chrom_strand_count:
            #     print("HP2_"+pb, "has 0 count for both P and M reads, probably not covered by any valid strand-seq reads")
            #     print("hap score default to 0")
            #     hp2_score = 0
            # else:
            hp2_p = chrom_strand_count["HP2_"+pb]["P"]
            hp2_m = chrom_strand_count["HP2_"+pb]["M"]
            hp2_score = (hp2_p-hp2_m)/(hp2_p+hp2_m)
            
            # PhaseQ_dict[pb] = abs(hp1_score - hp2_score) - 1 # not sure if this make sense
            PhaseQ_dict[pb] = abs(hp1_score - hp2_score)/2 # not sure if this make sense

        with open(f"{input_dir}/{chrom}/phasing_qual.pkl", "wb") as phaseQf:
            pickle.dump(PhaseQ_dict, phaseQf)

    
    #for contig in list(chrom_strand_count.keys())[:200]:
    # for contig in chrom_strand_count.keys():
    #     print(contig)
    #     print(chrom_strand_count[contig])


def plot_phaseQ_distribution(input_dir ,chroms=None):

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    if chroms is None:
        chroms = ["chr"+str(i) for i in range(1, 23)] # default to all autosomes

    total_phaseQ = list()
    
    for chrom in chroms:
        chrom_phaseQ = list()
        with open(f"{input_dir}/{chrom}/phasing_qual.pkl", "rb") as phaseQf:
            PhaseQ_dict = pickle.load(phaseQf)
            for PhaseQ in PhaseQ_dict.values():
                chrom_phaseQ.append(PhaseQ)
                total_phaseQ.append(PhaseQ)

        fig, ax = plt.subplots()
        counts, bins, patches = ax.hist(chrom_phaseQ, bins=20)
        for count, bin_left, bin_right in zip(counts, bins[:-1], bins[1:]):
            percentage = 100 * count / len(chrom_phaseQ)
            x = (bin_left + bin_right) / 2
            y = count
            ax.text(x, y + 0.2, f'{percentage:.1f}%', ha='center',fontsize=6)

        ax.set_title(f'Phasing Quality - {chrom}')
        ax.set_xlabel('PhaseQ')
        ax.set_ylabel('Phase block count')

        plt.savefig(f"{input_dir}/{chrom}/PhaseQ_distribution.pdf",bbox_inches='tight')
        plt.close(fig)

    fig, ax = plt.subplots()
    counts, bins, patches = ax.hist(total_phaseQ, bins=20)
    for count, bin_left, bin_right in zip(counts, bins[:-1], bins[1:]):
        percentage = 100 * count / len(total_phaseQ)
        x = (bin_left + bin_right) / 2
        y = count
        ax.text(x, y + 0.2, f'{percentage:.1f}%', ha='center',fontsize=6)

    ax.set_title('Phasing Quality')
    ax.set_xlabel('PhaseQ')
    ax.set_ylabel('Phase block count')

    plt.savefig(f"{input_dir}/aggregated_PhaseQ_distribution.pdf",bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--informative_cell','-c', required=True)
    parser.add_argument('--input_dir','-i', required=True)
    parser.add_argument('--chroms', nargs="+",required=False, default=None)
    parser.add_argument('--mapq', type=int, required=False, default=20)

    args = parser.parse_args()

    calculate_phasing_quality(args.informative_cell_pkl,
                              args.input_dir,
                              args.chroms,
                              args.mapq)

    plot_phaseQ_distribution(args.input_dir,args.chroms)

# calculate_phasing_quality("/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/2-align_to_T2T-HG002/informative_cell.pkl",
#                           "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/4-align_to_PhaseBlock/")

# plot_phaseQ_distribution("/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/4-align_to_PhaseBlock/")
