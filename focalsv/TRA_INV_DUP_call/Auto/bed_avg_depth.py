#!/usr/bin/env python3
import argparse
import subprocess
from tqdm import tqdm
from joblib import Parallel, delayed

def compute_avg_depth(bam_path, chrom, start, end):
    """Run `samtools depth` on one interval and return (chrom, start, end, avg_depth)."""
    region = f"{chrom}:{start}-{end}"

    cmd = "samtools depth '{0}' -r {1} | awk '{{sum+=$3}} END {{ print sum/NR}}'".format(bam_path, region)
    try:
        output = subprocess.check_output(cmd, shell=True, text=True)
        avg_depth = float(output.strip())
        # print("Average depth:", avg_depth)
        return  avg_depth
    except subprocess.CalledProcessError as e:
        print("Command failed:", e)
        return  0


def main():
    parser = argparse.ArgumentParser(
        description="Compute avg read depth per BED interval (in parallel)."
    )
    parser.add_argument("--bed", '-bed',    required=True, help="Input BED file")
    parser.add_argument("--out_bed", '-o',    required=True, help="Output BED file")
    parser.add_argument("--bam", '-bam',    required=True, help="Indexed BAM file")
    # parser.add_argument("--out", '-o',    default="regions_with_depth.bed",
    #                     help="Output BED (col4=avg depth)")
    parser.add_argument("--threads",'-t', type=int, default=4,
                        help="Number of parallel workers")
    args = parser.parse_args()

    # 1) load intervals
    intervals = []
    with open(args.bed) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.split()[:3]
            intervals.append((chrom, start, end))

    # 2) compute in parallel
    results = Parallel(n_jobs=args.threads)(
        delayed(compute_avg_depth)(args.bam, chrom, start, end)
        for chrom, start, end in tqdm(intervals)
    )

    results_left = Parallel(n_jobs=args.threads)(
        delayed(compute_avg_depth)(args.bam, chrom, int(start)-1000, start)
        for chrom, start, end in tqdm(intervals)
    )

    results_right = Parallel(n_jobs=args.threads)(
        delayed(compute_avg_depth)(args.bam, chrom, end, int(end)+1000)
        for chrom, start, end in tqdm(intervals)
    )

    # 3) write output
    with open(args.out_bed, "w") as out:
        i = 0
        with open(args.bed) as f:
            for line in f:
                a,b,c = results[i], results_left[i], results_right[i]
                i+=1
                out.write(f"{line[:-1]}\t{a}\t{b}\t{c}\t{2*a-b-c}\n")

if __name__ == "__main__":
    main()
