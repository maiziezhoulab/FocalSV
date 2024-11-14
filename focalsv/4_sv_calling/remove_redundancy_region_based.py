import os
import numpy as np
import networkx as nx
import edlib
from argparse import ArgumentParser
from utils import setup_logging  # Assuming `setup_logging` is available in utils.py
from tqdm import tqdm

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

import os

def merge_vcf(target_sv, regions, flanking, outvcf, vcf_prefix="dippav_variant_no_redundancy.vcf"):
    rg_list =[]
    folder_list = []
    cnt = 0
    with open(target_sv,'r') as f:
        for line in f:
            if line[0]!='#':
                data = line.split()
                pos = int(data[1])
                rg_start = pos - flanking
                rg_end = pos + flanking 
                if 'SVTYPE=DEL' in data[7]:
                    svlen = abs(len(data[3])-len(data[4]))
                    rg_end += svlen 
                rg_list.append((rg_start,rg_end))
                folder = 'Out%d_%d'%(cnt+1,pos)
                folder_list.append(folder)
                cnt+=1
    vcf_list = [os.path.join(regions, fd,'results/final_vcf', vcf_prefix) for fd in os.listdir(regions) if fd.startswith("Region")]
    header = []
    with open(vcf_list[0],'r') as f:
        for line in f:
            if line[0]=='#':
                header.append(line)
            else:
                break
    add_line1 = "##INFO=<ID=Region_start,Number=1,Type=Integer,Description=\"Region starting postion\">\n"
    add_line2 = "##INFO=<ID=Region_end,Number=1,Type=Integer,Description=\"Region ending postion\">\n"
    header = header[:-2]+[add_line1,add_line2]+header[-2:]
    lines = []
    cnt = 0 
    for i in range(len(vcf_list)):
        with open(vcf_list[i],'r') as f:
            for line in f:
                if line[0]!='#':
                    data = line.split()
                    data[2] = 'FocalSV%d'%cnt
                    cnt+=1
                    data[7] = data[7]+';Region_start=%d;Region_end=%d'%(rg_list[i][0],rg_list[i][1])
                    line = '\t'.join(data)+'\n'
                    lines.append(line)
    with open(outvcf,'w') as f:
        f.writelines(header + lines)

    return 

def sort_sig(sig_list):
    sorted_sig = []
    for i in range(1,23):
        chr_name = 'chr'+str(i)
        sig_list_chr = []
        for sig in sig_list:
            if sig[0]==chr_name:
                sig_list_chr.append(sig)
        sorted_sig.extend(sort_sig_per_chr(sig_list_chr))
    return sorted_sig
        

def sort_sig_per_chr(sig_list):
    pos_list = []
    for sig in sig_list:
        pos_list.append(sig[1])
    idx_list = np.argsort(pos_list)
    
    sorted_sig_list = []
    for idx in idx_list:
        sorted_sig_list.append(sig_list[idx])
    return sorted_sig_list

def vcf_to_sig(vcf_path):
    with open(vcf_path,'r') as f:
        s = f.readlines()
    del_sig = []
    ins_sig = []
    dc = {}
    header = []
    for line in s:
        if line[0]!='#':
            data = line.split()
            data[1] = int(data[1])
            data[3] = data[3].upper()
            data[4] = data[4].upper()
            dc[data[2]] = data
            if 'SVTYPE=DEL' in line:
                del_sig.append(data)
            elif 'SVTYPE=INS' in line:
                ins_sig.append(data)
        else:
            header.append(line)
    add_line = "##INFO=<ID=CollapseId,Number=1,Type=Integer,Description=\"collapse match ID\">\n"
    header = header[:-2]+[add_line]+header[-2:]
    return sort_sig(del_sig),sort_sig(ins_sig),dc,header

def edit_sim(seq1,seq2):
	
	scr = edlib.align(seq1, seq2)
	totlen = len(seq1) + len(seq2)
	sim = (totlen - scr["editDistance"]) / totlen
	# print('edlib', (totlen - scr["editDistance"]) / totlen)
	return sim

def jaro_sim(seq1, seq2):
	import jellyfish
	sim = jellyfish.jaro_distance(seq1,seq2)
	return sim
def get_size_sim(svlen1,svlen2):
	svlen1,svlen2 = abs(svlen1),abs(svlen2)
	sim = min(svlen1,svlen2)/max(svlen1,svlen2)
	return sim

def match_ins_one_pair(sig1,sig2,dist_thresh, size_sim_thresh, seq_sim_thresh):
    dist_ref = abs(sig2[1]-sig1[1])
    svlen1 = abs(len(sig1[3])-len(sig1[4]))
    svlen2 = abs(len(sig2[3])-len(sig2[4]))
    size_sim = get_size_sim(svlen1, svlen2)
    match_result = 0
    if (dist_ref <= dist_thresh) and (size_sim >= size_sim_thresh):
        seq_sim = edit_sim(sig1[4],sig2[4])
        if seq_sim>= seq_sim_thresh:
            match_result = 1
    return match_result

def get_reciprocal_overlap(sig1,sig2):
    svlen1 = abs(len(sig1[3])-len(sig1[4]))
    svlen2 = abs(len(sig2[3])-len(sig2[4]))
    start1 = sig1[1]
    start2 = sig2[1]
    end1 = start1+svlen1
    end2 = start2+svlen2
    
    overlap = (min(end1,end2)-max(start1,start2))/(max(svlen1,svlen2))
    return overlap 
    
def match_del_one_pair(sig1,sig2,dist_thresh, size_sim_thresh, overlap_thresh):
    dist_ref = abs(sig2[1]-sig1[1])
    svlen1 = abs(len(sig1[3])-len(sig1[4]))
    svlen2 = abs(len(sig2[3])-len(sig2[4]))
    size_sim = get_size_sim(svlen1, svlen2)
    match_result = 0
    if (dist_ref <= dist_thresh) and (size_sim >= size_sim_thresh):
        overlap  = get_reciprocal_overlap(sig1,sig2)
        if overlap>= overlap_thresh:
            match_result = 1
    return match_result
            
def match_del_chr(sig_list,dist_thresh, size_sim_thresh, overlap_thresh):
    cluster_list = [-1]*len(sig_list)
    links = []
    for i in range(len(sig_list)):
        sig1 = sig_list[i]
        pos1 = sig1[1]
        window = (pos1-dist_thresh, pos1+dist_thresh)
        comp_sig_list = []
        for j in range(len(sig_list)):
            pos2 = sig_list[j][1]
            if pos2>window[1]:
                break
            if (i!=j) and (window[0]<=pos2<=window[1]):
                comp_sig_list.append(sig_list[j]) 
        for sig2 in comp_sig_list:
            match_result = match_del_one_pair(sig1,sig2,dist_thresh, size_sim_thresh, overlap_thresh)
            if match_result == 1:
                links.append((sig1[2],sig2[2]))
    return links


def match_del(sig_list,dist_thresh, size_sim_thresh, overlap_thresh):
    links = []
    for i in tqdm(range(1,23)):
        chr_name = 'chr'+str(i)
        sig_list_chr = []
        for sig in sig_list:
            if sig[0]==chr_name:
                sig_list_chr.append(sig)
        links.extend(match_del_chr(sig_list_chr,dist_thresh, size_sim_thresh, overlap_thresh))
    G = nx.Graph()
    G.add_edges_from(links)
    components = nx.connected_components(G)
    nodes_list = [sorted(nodes) for nodes in components]
    return nodes_list

def match_ins_chr(sig_list,dist_thresh, size_sim_thresh, seq_sim_thresh):
    cluster_list = [-1]*len(sig_list)
    links = []
    for i in range(len(sig_list)):
        sig1 = sig_list[i]
        pos1 = sig1[1]
        window = (pos1-dist_thresh, pos1+dist_thresh)
        comp_sig_list = []
        for j in range(len(sig_list)):
            pos2 = sig_list[j][1]
            if pos2>window[1]:
                break
            if (i!=j) and (window[0]<=pos2<=window[1]):
                comp_sig_list.append(sig_list[j])
                
        for sig2 in comp_sig_list:
            match_result = match_ins_one_pair(sig1,sig2,dist_thresh, size_sim_thresh, seq_sim_thresh)
            if match_result == 1:
                links.append((sig1[2],sig2[2]))
    return links
            
                
def match_ins(sig_list,dist_thresh, size_sim_thresh, seq_sim_thresh):
    links = []
    for i in tqdm(range(1,23)):
        chr_name = 'chr'+str(i)
        sig_list_chr = []
        for sig in sig_list:
            if sig[0]==chr_name:
                sig_list_chr.append(sig)
        links.extend(match_ins_chr(sig_list_chr,dist_thresh, size_sim_thresh, seq_sim_thresh))
    G = nx.Graph()
    G.add_edges_from(links)
    components = nx.connected_components(G)
    nodes_list = [sorted(nodes) for nodes in components]
    return nodes_list
        
                    
def pick_best_sv_one_cluster(vcf_dc,index_list):
    sig_list = [vcf_dc[idx] for idx in index_list]
    svtype = sig_list[0][7].split('SVTYPE=')[1].split(';')[0]
    ll = [abs(len(sig[3])-len(sig[4])) for sig in sig_list]

    if svtype == 'INS':
        sv_center_list = np.array([ int(sig[1]) for sig in sig_list])
    else:
        sv_center_list = np.array([ int((int(sig_list[i][1])+ll[i])/2) for i in range(len(sig_list))])

    region_start_list = np.array([int(sig[7].split("Region_start=")[1].split(';')[0]) for sig in sig_list])
    region_end_list = np.array([int(sig[7].split("Region_end=")[1].split(';')[0]) for sig in sig_list])
    region_ct_list = (region_start_list + region_end_list) / 2
    dist_to_rgcenter = abs(sv_center_list - region_ct_list )
    min_dist = dist_to_rgcenter.min()
    md_ids = np.where(dist_to_rgcenter==min_dist)[0]
    ll_md = [ll[i] for i in md_ids]
    max_len_under_md_constraint = max(ll_md)
    best_idx = ll.index(max_len_under_md_constraint)
    return index_list[best_idx]
        
        
def pick_best_sv(vcf_dc,nodes_list):
    retain_index = {}
    remove_index = {}
    
    for i in range(len(nodes_list)):
        index_list = list(nodes_list[i])
        best_index = pick_best_sv_one_cluster(vcf_dc,index_list)
        retain_index[best_index] = i
#         retain_index[0].append(best_index)
#         retain_index[1].append(i)
        for idx in index_list:
            if idx!=best_index:
                remove_index[idx] = i
#                 remove_index[0].append((idx,i))

    return retain_index,remove_index
                
        
def write_one_vcf(sig_list,header,path):
    with open(path,'w') as f:
        f.writelines(header)
        for sig in sig_list:
            sig[1]=str(sig[1])
            line = '\t'.join(sig)+'\n'
            f.write(line)
    return

def write_vcf(output_dir,prefix,header,retain_index_del,remove_index_del,retain_index_ins,remove_index_ins):
    retain_sig = []
    remove_sig = []
    for idx,sig in vcf_dc.items():
        if idx in retain_index_del:
            sig[7] = sig[7]+";CollapseId=DEL%d"%retain_index_del[idx]
            retain_sig.append(sig)
        elif idx in retain_index_ins:
            sig[7] = sig[7]+";CollapseId=INS%d"%retain_index_ins[idx]
            retain_sig.append(sig)
        elif idx in remove_index_del:
            sig[7] = sig[7]+";CollapseId=DEL%d"%remove_index_del[idx]
            remove_sig.append(sig)
        elif idx in remove_index_ins:
            sig[7] = sig[7]+";CollapseId=INS%d"%remove_index_ins[idx]
            remove_sig.append(sig)
        else:
            retain_sig.append(sig)
    retain_sig = sort_sig(retain_sig )
    remove_sig = sort_sig(remove_sig )
    rd_path = output_dir + "/"+ prefix + '_redundancy.vcf'
    nrd_path = output_dir + "/"+ prefix + '_no_redundancy.vcf'
    write_one_vcf(remove_sig,header,rd_path)
    write_one_vcf(retain_sig,header,nrd_path)
    print("original %d lines"%(len(retain_sig)+len(remove_sig)))
    print("new vcf %d lines"%len(retain_sig))
    print("redundancy %d lines"%len(remove_sig))
    return 

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
    logger = setup_logging("4_remove_redundancy", output_dir)

    # Merge VCFs
    vcf_path = os.path.join(output_dir, 'variants.vcf')
    merge_vcf(target_sv, regions, flanking, vcf_path, vcf_prefix)

    # Process VCF to get signatures
    prefix = 'FocalSV_variant'
    del_sig, ins_sig, vcf_dc, header = vcf_to_sig(vcf_path)

    # Perform matching and redundancy removal
    links_del = match_del(del_sig, dist_thresh_del, size_sim_thresh_del, overlap_thresh)
    links_ins = match_ins(ins_sig, dist_thresh, size_sim_thresh, seq_sim_thresh)
    retain_index_del, remove_index_del = pick_best_sv(vcf_dc, links_del)
    retain_index_ins, remove_index_ins = pick_best_sv(vcf_dc, links_ins)

    # Write the final VCFs
    write_vcf(output_dir, prefix, header, retain_index_del, remove_index_del, retain_index_ins, remove_index_ins)
