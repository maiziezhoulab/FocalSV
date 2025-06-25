# NOTE: this script is specifically designed for collecting cw informative cells from
# the strand-seq to HG002-T2T assembly alignment result. 
# DO NOT use it directly on strand-seq to reference alignment cw ratio results

import pickle

def collect_informative_cell(cell_cw_ratio_file, out_dir="./"):

    informative_cell = dict()
    for i in range(1, 23): # currently, we hard code the chromosome range (auto chromosome)
        informative_cell["chr"+str(i)] = list()

    with open(cell_cw_ratio_file, "rb") as f:
        cell_cw_ratio = pickle.load(f)

    for cell, contig_cw_ratio in cell_cw_ratio.items():
        for chrom in informative_cell.keys():
            mat_cw = contig_cw_ratio[chrom+"_MATERNAL"]["cw_ratio"]
            pat_cw = contig_cw_ratio[chrom+"_PATERNAL"]["cw_ratio"]

            if mat_cw is None or pat_cw is None:
                continue
            elif mat_cw * pat_cw < 0 and min(abs(mat_cw), abs(pat_cw)) > 0.75: # mat and pat has different W and C preference
                if mat_cw > 0:
                    informative_cell[chrom].append((cell, "M")) # forward strand (W) come from maternal hap
                elif pat_cw > 0:
                    informative_cell[chrom].append((cell, "P")) # forward strand (W) come from paternal hap
            else:
                continue

    for chrom, cells in informative_cell.items():
        if len(cells) == 0:
            print("[WARNING]: No informative cell found for", chrom)

    with open(out_dir+"/informative_cell.pkl", "wb") as icf:
        pickle.dump(informative_cell, icf)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--cell_cw_ratio','-i',)
    parser.add_argument('--out_dir','-o_dir', default="./")

    args = parser.parse_args()

    collect_informative_cell(args.cell_cw_ratio)