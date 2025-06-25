import pysam
import pickle

def compute_cw_ratio(bam_path, min_mapq=20):
    bam = pysam.AlignmentFile(bam_path, "rb")
    strand_counts = {}

    for contig in bam.references:
        W, C = 0, 0
        for read in bam.fetch(contig): #NOTE: We assume that the bam was indexed 
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if not read.is_read1:
                continue  # Strand-seq: only use R1

            if read.is_reverse:
                W += 1
            else:
                C += 1

        total = W + C
        if total == 0:
            ratio = None
        else:
            ratio = (C - W) / total

        strand_counts[contig] = {'W': W, 'C': C, 'cw_ratio': ratio}

    return strand_counts


def main(cell_bam_list, out_dir="./", min_mapq=20):

    cw_ratio_dict = dict()
    with open(cell_bam_list, "r") as f:
        for line in f:
            cell, bam = line.rstrip("\n").split("\t")
            cw_ratio_dict[cell] = compute_cw_ratio(bam, min_mapq)
            print(cell, "Done")

    
    with open(out_dir+"/cell_cw_ratio.pkl", "wb") as cwf:
        pickle.dump(cw_ratio_dict, cwf)

    print("All Done")



if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--cell_bam_list','-i',)
    parser.add_argument('--out_dir','-o_dir', default="./")
    parser.add_argument('--min_mapq', type=int, default=20)

    args = parser.parse_args()

    main(args.cell_bam_list, args.out_dir, args.min_mapq)