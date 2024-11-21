from argparse import ArgumentParser
import os
from collections import defaultdict
import pysam
from utils import setup_logging  # Import the logging function

parser = ArgumentParser(description="Output FASTA files according to phase blocks.")
parser.add_argument('--out_dir', '-o', required=True)

def output_fa(fd, out_dir, logger):
    logger.info(f"Processing phased BAM: {fd}")
    total, unphased = 0, 0

    phased_bam = os.path.join(fd, "region_phased.bam")
    phase_blocks = defaultdict(list)
    unphased_block = []

    # Open BAM file and iterate through reads
    samfile = pysam.AlignmentFile(phased_bam, "rb")
    samiter = samfile.fetch(until_eof=True)

    for read in samiter:
        total += 1
        try:
            ps = read.get_tag('PS')
            hp = read.get_tag('HP')
            key = f"{ps}_{hp}"
            phase_blocks[key].append((read.qname, read.seq, read.reference_start, read.reference_end))
        except KeyError:
            unphased += 1
            unphased_block.append((read.qname, read.seq, read.reference_start, read.reference_end))

    # Calculate phase block boundaries
    start, end = defaultdict(int), defaultdict(int)
    for key, values in phase_blocks.items():
        pb = int(key.split("_")[0])
        start[pb] = min(start.get(pb, float('inf')), min([s[2] for s in values]))
        end[pb] = max(end.get(pb, float('-inf')), max([s[3] for s in values]))

    # Allocate unphased reads
    for readname, readseq, reads, reade in unphased_block:
        if len(phase_blocks) == 2:
            for key in phase_blocks:
                phase_blocks[key].append((readname, readseq, reads, reade))
        else:
            max_overlap, max_pb = -float('inf'), None
            for pb in end:
                overlap = min(reade, end[pb]) - max(reads, start[pb])
                if overlap > max_overlap:
                    max_overlap = overlap
                    max_pb = pb
            if max_pb is not None:
                phase_blocks[f"{max_pb}_1"].append((readname, readseq, reads, reade))
                phase_blocks[f"{max_pb}_2"].append((readname, readseq, reads, reade))

    # Write phased reads to FASTA files
    for key, values in phase_blocks.items():
        output_file = f"PS{key.split('_')[0]}_hp{key.split('_')[1]}.fa"
        output_path = os.path.join(fd, output_file)
        logger.info(f"Writing {output_path} with {len(values)} reads")

        with open(output_path, "w") as fw:
            dup_reads = set()
            for readname, readseq, _, _ in values:
                if readname not in dup_reads:
                    fw.write(f">{readname}\n{readseq}\n")
                    dup_reads.add(readname)

if __name__ == "__main__":
    args = parser.parse_args()
    out_dir_general = args.out_dir
    out_dir = os.path.join(out_dir_general, "regions")

    # Initialize logger
    logger = setup_logging("2_output_fas", out_dir_general)

    # Process each region folder
    fds = [os.path.join(out_dir, fd) for fd in os.listdir(out_dir) if fd.startswith("Region")]
    for fd in fds:
        output_fa(fd, out_dir, logger)
