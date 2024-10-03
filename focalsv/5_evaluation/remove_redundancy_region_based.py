import os
import numpy as np
import networkx as nx
import edlib
from argparse import ArgumentParser
from utils import setup_logging  # Assuming `setup_logging` is available in utils.py

parser = ArgumentParser(description="Remove redundancy based on region information")
parser.add_argument('--target_sv', '-t', required=True)
parser.add_argument('--regions', '-rg', required=True)
parser.add_argument('--output_dir', '-o', required=True)
parser.add_argument('--flanking', '-fl', type=int, default=50000, help='Flanking region size')
parser.add_argument('--dist_thresh', '-r', type=int, default=500)
parser.add_argument('--dist_thresh_del', '-rd', type=int, default=3000)
parser.add_argument('--overlap_thresh', '-O', type=float, default=0)
parser.add_argument('--size_sim_thresh', '-P', type=float, default=0.5)
parser.add_argument('--size_sim_thresh_del', '-Pd', type=float, default=0.1)
parser.add_argument('--seq_sim_thresh', '-p', type=float, default=0.5)
parser.add_argument('--vcf_prefix', '-vcf', default='dippav_variant_no_redundancy.vcf')

args = parser.parse_args()

def merge_vcf(target_sv, regions, flanking, outvcf, vcf_prefix, logger):
    """
    Merge VCF files based on flanking regions.
    """
    logger.info(f"Merging VCFs from target_sv: {target_sv}, regions: {regions}")
    rg_list, folder_list = [], []
    cnt = 0
    
    try:
        with open(target_sv, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data = line.split()
                pos = int(data[1])
                rg_start = pos - flanking
                rg_end = pos + flanking
                if 'SVTYPE=DEL' in data[7]:
                    svlen = abs(len(data[3]) - len(data[4]))
                    rg_end += svlen
                rg_list.append((rg_start, rg_end))
                folder_list.append(f'Out{cnt+1}_{pos}')
                cnt += 1
    except Exception as e:
        logger.error(f"Error processing target_sv: {e}")
        raise

    vcf_list = [os.path.join(regions, folder, 'results/final_vcf', vcf_prefix) for folder in folder_list]
    
    try:
        header = []
        with open(vcf_list[0], 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header.append(line)
                else:
                    break

        header.insert(-2, "##INFO=<ID=Region_start,Number=1,Type=Integer,Description=\"Region starting position\">\n")
        header.insert(-2, "##INFO=<ID=Region_end,Number=1,Type=Integer,Description=\"Region ending position\">\n")

        lines = []
        for i, vcf in enumerate(vcf_list):
            with open(vcf, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        data = line.split()
                        data[2] = f'AquilaSV{cnt}'
                        cnt += 1
                        data[7] += f';Region_start={rg_list[i][0]};Region_end={rg_list[i][1]}'
                        lines.append('\t'.join(data) + '\n')

        with open(outvcf, 'w') as f:
            f.writelines(header + lines)
        logger.info(f"VCF merged and written to {outvcf}")
    except Exception as e:
        logger.error(f"Error merging VCFs: {e}")
        raise

def sort_sig(sig_list):
    """
    Sort signal list by chromosome and position.
    """
    sorted_sig = []
    for i in range(1, 23):
        chr_name = 'chr' + str(i)
        sig_list_chr = [sig for sig in sig_list if sig[0] == chr_name]
        sorted_sig.extend(sort_sig_per_chr(sig_list_chr))
    return sorted_sig

def sort_sig_per_chr(sig_list):
    """
    Sort signal list within a single chromosome by position.
    """
    pos_list = [sig[1] for sig in sig_list]
    idx_list = np.argsort(pos_list)
    return [sig_list[idx] for idx in idx_list]

def vcf_to_sig(vcf_path):
    """
    Convert VCF data to signature for DEL/INS.
    """
    del_sig, ins_sig, dc = [], [], {}
    header = []
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header.append(line)
                else:
                    data = line.split()
                    data[1] = int(data[1])
                    data[3] = data[3].upper()
                    data[4] = data[4].upper()
                    dc[data[2]] = data
                    if 'SVTYPE=DEL' in line:
                        del_sig.append(data)
                    elif 'SVTYPE=INS' in line:
                        ins_sig.append(data)
        header.append("##INFO=<ID=CollapseId,Number=1,Type=Integer,Description=\"collapse match ID\">\n")
    except Exception as e:
        logger.error(f"Error reading VCF: {e}")
        raise

    return sort_sig(del_sig), sort_sig(ins_sig), dc, header

def edit_sim(seq1, seq2):
    """
    Calculate sequence similarity using edit distance.
    """
    scr = edlib.align(seq1, seq2)
    totlen = len(seq1) + len(seq2)
    return (totlen - scr["editDistance"]) / totlen

def get_size_sim(svlen1, svlen2):
    """
    Calculate size similarity between two SV lengths.
    """
    return min(abs(svlen1), abs(svlen2)) / max(abs(svlen1), abs(svlen2))

def match_ins_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, seq_sim_thresh):
    """
    Match one pair of insertion signals.
    """
    dist_ref = abs(sig2[1] - sig1[1])
    svlen1, svlen2 = abs(len(sig1[3]) - len(sig1[4])), abs(len(sig2[3]) - len(sig2[4]))
    size_sim = get_size_sim(svlen1, svlen2)

    if dist_ref <= dist_thresh and size_sim >= size_sim_thresh:
        seq_sim = edit_sim(sig1[4], sig2[4])
        return seq_sim >= seq_sim_thresh
    return False

def match_del_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, overlap_thresh):
    """
    Match one pair of deletion signals.
    """
    dist_ref = abs(sig2[1] - sig1[1])
    size_sim = get_size_sim(len(sig1[3]) - len(sig1[4]), len(sig2[3]) - len(sig2[4]))

    if dist_ref <= dist_thresh and size_sim >= size_sim_thresh:
        overlap = get_reciprocal_overlap(sig1, sig2)
        return overlap >= overlap_thresh
    return False

def get_reciprocal_overlap(sig1, sig2):
    """
    Calculate reciprocal overlap between two signals.
    """
    start1, end1 = sig1[1], sig1[1] + abs(len(sig1[3]) - len(sig1[4]))
    start2, end2 = sig2[1], sig2[1] + abs(len(sig2[3]) - len(sig2[4]))
    return (min(end1, end2) - max(start1, start2)) / max(abs(len(sig1[3]) - len(sig1[4])), abs(len(sig2[3]) - len(sig2[4])))

def match_del_chr(sig_list, dist_thresh, size_sim_thresh, overlap_thresh):
    """
    Match deletion signals within a single chromosome.
    """
    links = []
    for i, sig1 in enumerate(sig_list):
        pos1 = sig1[1]
        window = (pos1 - dist_thresh, pos1 + dist_thresh)
        comp_sig_list = [sig_list[j] for j in range(len(sig_list)) if i != j and window[0] <= sig_list[j][1] <= window[1]]

        for sig2 in comp_sig_list:
            if match_del_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, overlap_thresh):
                links.append((sig1[2], sig2[2]))  # Linking SV IDs
    return links

def match_del(sig_list, dist_thresh, size_sim_thresh, overlap_thresh):
    """
    Match deletions across all chromosomes.
    """
    links = []
    for i in range(1, 23):
        chr_name = f'chr{i}'
        sig_list_chr = [sig for sig in sig_list if sig[0] == chr_name]
        links.extend(match_del_chr(sig_list_chr, dist_thresh, size_sim_thresh, overlap_thresh))

    G = nx.Graph()
    G.add_edges_from(links)
    components = nx.connected_components(G)
    nodes_list = [sorted(nodes) for nodes in components]
    return nodes_list

def match_ins(sig_list, dist_thresh, size_sim_thresh, seq_sim_thresh):
    """
    Match insertions across all chromosomes.
    """
    links = []
    for i in range(1, 23):
        chr_name = f'chr{i}'
        sig_list_chr = [sig for sig in sig_list if sig[0] == chr_name]
        for i, sig1 in enumerate(sig_list_chr):
            pos1 = sig1[1]
            window = (pos1 - dist_thresh, pos1 + dist_thresh)
            comp_sig_list = [sig_list_chr[j] for j in range(len(sig_list_chr)) if i != j and window[0] <= sig_list_chr[j][1] <= window[1]]

            for sig2 in comp_sig_list:
                if match_ins_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, seq_sim_thresh):
                    links.append((sig1[2], sig2[2]))  # Linking SV IDs

    G = nx.Graph()
    G.add_edges_from(links)
    components = nx.connected_components(G)
    nodes_list = [sorted(nodes) for nodes in components]
    return nodes_list

def pick_best_sv_one_cluster(vcf_dc, index_list):
    """
    Select the best SV from a cluster based on proximity to the region center and length.
    """
    sig_list = [vcf_dc[idx] for idx in index_list]
    svtype = sig_list[0][7].split('SVTYPE=')[1].split(';')[0]
    ll = [abs(len(sig[3]) - len(sig[4])) for sig in sig_list]

    # Center positions and region centers
    sv_center_list = np.array([int(sig[1]) for sig in sig_list]) if svtype == 'INS' else \
                     np.array([int((int(sig[1]) + ll[i]) / 2) for i, sig in enumerate(sig_list)])
    region_start_list = np.array([int(sig[7].split("Region_start=")[1].split(';')[0]) for sig in sig_list])
    region_end_list = np.array([int(sig[7].split("Region_end=")[1].split(';')[0]) for sig in sig_list])
    region_ct_list = (region_start_list + region_end_list) / 2

    # Find the closest SV to the region center
    dist_to_rgcenter = abs(sv_center_list - region_ct_list)
    min_dist = dist_to_rgcenter.min()
    md_ids = np.where(dist_to_rgcenter == min_dist)[0]
    ll_md = [ll[i] for i in md_ids]
    max_len_under_md_constraint = max(ll_md)
    best_idx = ll.index(max_len_under_md_constraint)
    return index_list[best_idx]

def pick_best_sv(vcf_dc, nodes_list):
    """
    Pick the best SV from each cluster and return indices for retention or removal.
    """
    retain_index, remove_index = {}, {}

    for index_list in nodes_list:
        best_index = pick_best_sv_one_cluster(vcf_dc, index_list)
        retain_index[best_index] = index_list
        for idx in index_list:
            if idx != best_index:
                remove_index[idx] = index_list

    return retain_index, remove_index

def write_vcf(output_dir, prefix, header, retain_index_del, remove_index_del, retain_index_ins, remove_index_ins):
    """
    Write the final VCF file after collapsing SVs.
    """
    retain_sig, remove_sig = [], []

    for idx, sig in vcf_dc.items():
        if idx in retain_index_del:
            sig[7] += f";CollapseId=DEL{retain_index_del[idx]}"
            retain_sig.append(sig)
        elif idx in retain_index_ins:
            sig[7] += f";CollapseId=INS{retain_index_ins[idx]}"
            retain_sig.append(sig)
        elif idx in remove_index_del:
            sig[7] += f";CollapseId=DEL{remove_index_del[idx]}"
            remove_sig.append(sig)
        elif idx in remove_index_ins:
            sig[7] += f";CollapseId=INS{remove_index_ins[idx]}"
            remove_sig.append(sig)
        else:
            retain_sig.append(sig)

    retain_sig = sort_sig(retain_sig)
    remove_sig = sort_sig(remove_sig)

    # Write the retained and removed signals to separate VCF files
    rd_path = os.path.join(output_dir, f"{prefix}_redundancy.vcf")
    nrd_path = os.path.join(output_dir, f"{prefix}_no_redundancy.vcf")

    with open(rd_path, 'w') as f:
        f.writelines(header)
        for sig in remove_sig:
            sig[1] = str(sig[1])
            f.write('\t'.join(sig) + '\n')

    with open(nrd_path, 'w') as f:
        f.writelines(header)
        for sig in retain_sig:
            sig[1] = str(sig[1])
            f.write('\t'.join(sig) + '\n')

    print(f"Original: {len(retain_sig) + len(remove_sig)} lines")
    print(f"New VCF: {len(retain_sig)} lines")
    print(f"Redundancy: {len(remove_sig)} lines")

if __name__ == "__main__":
    # Extract command line arguments
    target_sv = args.target_sv
    regions = args.regions
    flanking = args.flanking
    output_dir = args.output_dir
    vcf_prefix = args.vcf_prefix
    dist_thresh = args.dist_thresh
    dist_thresh_del = args.dist_thresh_del
    overlap_thresh = args.overlap_thresh
    size_sim_thresh = args.size_sim_thresh
    size_sim_thresh_del = args.size_sim_thresh_del
    seq_sim_thresh = args.seq_sim_thresh

    # Initialize logger
    logger = setup_logging("remove_redundancy", output_dir)

    # Merge VCFs
    vcf_path = os.path.join(output_dir, 'variants.vcf')
    merge_vcf(target_sv, regions, flanking, vcf_path, vcf_prefix, logger)

    # Process VCF to get signatures
    prefix = 'AquilaSV_variant'
    del_sig, ins_sig, vcf_dc, header = vcf_to_sig(vcf_path)

    # Perform matching and redundancy removal
    links_del = match_del(del_sig, dist_thresh_del, size_sim_thresh_del, overlap_thresh)
    links_ins = match_ins(ins_sig, dist_thresh, size_sim_thresh, seq_sim_thresh)
    retain_index_del, remove_index_del = pick_best_sv(vcf_dc, links_del)
    retain_index_ins, remove_index_ins = pick_best_sv(vcf_dc, links_ins)

    # Write the final VCFs
    write_vcf(output_dir, prefix, header, retain_index_del, remove_index_del, retain_index_ins, remove_index_ins)
