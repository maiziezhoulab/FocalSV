# NOTE: this script is for extracting HG002-T2T contig regions from the
# ref regions they were aligned to

import pysam
from collections import defaultdict
import pickle
import os
from Bio import SeqIO
from tqdm import tqdm

def rev_strd_coord(start, end, c_len):
    return c_len-end, c_len-start

def read_fai(fai):
    fai_dict = dict()
    with open(fai, "r") as f:
        for line in f:
            if line[0] != "#":
                line = line.rstrip("\n").split("\t")
                fai_dict[line[0]] = [int(i) for i in line[1:]]
    return fai_dict

def merge_intervals(intervals, max_gap):
    """
    intervals: List of (start, end) tuples or lists
    max_gap  : max distance
    Returns a merged list of intervals.
    """
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    cur_s, cur_e = intervals[0]
    
    for s, e in intervals[1:]:
        if s - cur_e <= max_gap:
            cur_e = max(cur_e, e)
        else:
            merged.append([cur_s, cur_e])
            cur_s, cur_e = s, e
    merged.append([cur_s, cur_e])
    return merged

def get_t2t_coord_in_region(segment, region_start, region_end):
    """
    For a given aligned segment and a target reference region,
    return the overlapping portion on both the reference and contig (query),
    in absolute coordinates relative to the full contig.

    Parameters:
        segment (pysam.AlignedSegment): The aligned segment
        region_start (int): Start of the reference region (0-based)
        region_end (int): End of the reference region (0-based, exclusive)
    """
    ref_pos = segment.reference_start
    query_pos = 0

    # contig_offset = segment.query_alignment_start  # alignment start on full contig
    # removed because this wil cause confusions and wrong shifs when hard-clips present
    # query_alignment_start : start index of the aligned query portion of the sequence (0-based, inclusive).
    #                         This the index of the first base in query_sequence that is not soft-clipped.

    query_overlap_start = None
    query_overlap_end = None
    ref_overlap_start = None
    ref_overlap_end = None

    seg_ref_span = 0

    for op, length in segment.cigartuples:
        if op in (0, 7, 8):  # M, =, X: consumes both
            ref_end = ref_pos + length
            query_end = query_pos + length

            if ref_end > region_start and ref_pos < region_end:
                # Overlaps with region
                if ref_overlap_start is None:
                    ref_overlap_start = max(ref_pos, region_start)
                    offset = ref_overlap_start - ref_pos
                    query_overlap_start = query_pos + offset

                ref_overlap_end = min(ref_end, region_end)

                query_overlap_end = query_pos + (ref_overlap_end - ref_pos)

            ref_pos = ref_end
            query_pos = query_end
            seg_ref_span += length

        elif op in (2, 3):  # D, N: consumes ref only
            ref_end = ref_pos + length

            if ref_end > region_start and ref_pos < region_end:
                # Overlaps with region
                if ref_overlap_start is None:
                    ref_overlap_start = max(ref_pos, region_start)
                    query_overlap_start = query_pos  # query does not advance

                ref_overlap_end = min(ref_end, region_end)
                query_overlap_end = query_pos  # still no change

            ref_pos = ref_end
            seg_ref_span += length

        elif op in (1, 4, 5):  # I, S, H: consumes query only
            query_pos += length

    if ref_overlap_start is None:
        return None  # No overlap with region

    # strand = "-" if segment.is_reverse else "+"
    # contig_overlap_start = contig_offset + query_overlap_start
    # contig_overlap_end = contig_offset + query_overlap_end
    contig_overlap_start = query_overlap_start
    contig_overlap_end = query_overlap_end

    # ref_overlap_span = ref_overlap_end - ref_overlap_start
    #contig_overlap_span = contig_overlap_end - contig_overlap_start

    return segment.query_name, [contig_overlap_start, contig_overlap_end], segment.is_reverse


def collect_t2t_contig_coords(t2t_to_ref_bam, ref_regions, t2t_fai):

    t2t_region_dict = defaultdict(list)

    bam = pysam.AlignmentFile(t2t_to_ref_bam, "rb")
    for chrom, r_start, r_end in tqdm(ref_regions):
        for segment in bam.fetch(chrom, r_start, r_end):
            if not segment.is_secondary: # exclude secondary alignment
                c_name, coords, is_reverse = get_t2t_coord_in_region(segment, r_start, r_end)

                if chrom not in c_name: # exclude contigs that from other chromosomes (not sure if make sense)
                    continue

                if is_reverse:
                    coords = rev_strd_coord(*coords, t2t_fai[c_name][0])

                t2t_region_dict[c_name].append(coords)

    for c_name, intervals in t2t_region_dict.items():
        t2t_region_dict[c_name] = merge_intervals(intervals, max_gap=100000) # 100kb hard coded merge thresh

    return t2t_region_dict


def collect_strand_seq_read_name(strand_to_t2t_bam, t2t_region_dict, chrom):

    read_names = set()

    bam = pysam.AlignmentFile(strand_to_t2t_bam, "rb")

    for contig in [chrom+"_MATERNAL", chrom+"_PATERNAL"]:
        for start, end in t2t_region_dict[contig]:
            for read in bam.fetch(contig, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < 20: # hard coded mapqual threshold
                    continue
                if not read.is_read1:
                    continue  # Strand-seq: only use R1

                read_names.add(read.query_name)

    return read_names


def main(t2t_fai_f,
         ref_regions_bed,
         informative_cell_pkl,
         t2t_to_ref_bam,
         strand_to_t2t_dir,
         strand_seq_dir,
         out_dir,
         rewrite_t2t_region = False):
    
    t2t_fai = read_fai(t2t_fai_f)

    ref_regions = list()
    with open(ref_regions_bed, "r") as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            start = int(start)
            end = int(end)
            ref_regions.append([chrom, start, end])

    with open(informative_cell_pkl, "rb") as f:
        informative_cell = pickle.load(f)

    print("Extracting T2T regions")
    if rewrite_t2t_region or not os.path.isfile(f"{out_dir}/t2t_regions.pkl"):
        t2t_region_dict = collect_t2t_contig_coords(t2t_to_ref_bam, ref_regions, t2t_fai)
        with open(f"{out_dir}/t2t_regions.pkl", "wb") as f:
            pickle.dump(t2t_region_dict, f)
        print(f"T2T regions saved to {out_dir}/t2t_regions.pkl")
    else:
        print(f"T2T regions loaded from existing file {out_dir}/t2t_regions.pkl")
        with open(f"{out_dir}/t2t_regions.pkl", "rb") as f:
            t2t_region_dict = pickle.load(f)

    print("Extracting Strand seq reads")
    for chrom, cells in informative_cell.items():
        print("Processing", chrom, "...")
        os.makedirs(f"{out_dir}/{chrom}/", exist_ok=True)
        for cell, _ in tqdm(cells):
            if os.path.isfile(f"{strand_seq_dir}/{cell}_1.idx") and os.path.isfile(f"{strand_seq_dir}/{cell}_2.idx"):
                cell_1_idx = SeqIO.index_db(f"{strand_seq_dir}/{cell}_1.idx")
                cell_2_idx = SeqIO.index_db(f"{strand_seq_dir}/{cell}_2.idx")
            else:
                cell_1_idx = SeqIO.index_db(f"{strand_seq_dir}/{cell}_1.idx",
                                            f"{strand_seq_dir}/{cell}_1_sequence.fastq.gz",
                                            "fastq")
                cell_2_idx = SeqIO.index_db(f"{strand_seq_dir}/{cell}_2.idx",
                                            f"{strand_seq_dir}/{cell}_2_sequence.fastq.gz",
                                            "fastq")
                
            read_names = collect_strand_seq_read_name(f"{strand_to_t2t_dir}/{cell}.bam", t2t_region_dict, chrom)
            with open(f"{out_dir}/{chrom}/{cell}_1.fastq", "w") as f1, open(f"{out_dir}/{chrom}/{cell}_2.fastq", "w") as f2:
                for read_name in read_names:
                    f1.write(cell_1_idx[read_name].format("fastq"))
                    f2.write(cell_2_idx[read_name].format("fastq"))


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--t2t_fai', required=True)
    parser.add_argument('--sv_region_bed', required=True)
    parser.add_argument('--informative_cell', required=True)
    parser.add_argument('--t2t_bam', required=True)
    parser.add_argument('--strand_to_t2t_dir', required=True)
    parser.add_argument('--strand_seq_dir', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--rewrite_t2t_region', action='store_true')

    args = parser.parse_args()
    
    main(args.t2t_fai,
        args.sv_region_bed,
        args.informative_cell,
        args.t2t_bam,
        args.strand_to_t2t_dir,
        args.strand_seq_dir,
        args.out_dir,
        args.rewrite_t2t_region)



# main("/data/maiziezhou_lab/Datasets/Assemblies/HG002-T2T.v1.1.fasta.fai",
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/4-align_to_PhaseBlock/SV_regions.bed", # we will do this later with multi-threading
#     # "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/4-align_to_PhaseBlock/chr1_SV_regions.bed", # for now we test on chr1 for speed consideration
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/2-align_to_T2T-HG002/informative_cell.pkl",
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/3-align_T2T-HG002_to_ref/HG002-T2T_hg19.bam",
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/2-align_to_T2T-HG002/",
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/0-strand-seq_reads",
#     "/data/maiziezhou_lab/Yichen/Projects/FocalSV/Strand-seq/4-align_to_PhaseBlock",)