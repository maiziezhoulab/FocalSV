import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--bam_BL','-bl')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--bed_file','-bed')
parser.add_argument('--output_dir','-o')
parser.add_argument('--n_thread','-t', type = int, default = 10 )
# parser.add_argument('--prefix','-p', default="sample")
args = parser.parse_args()
# bam_BL = args.bam_BL
bam_TM = args.bamfile
bed_file = args.bed_file
output_dir = args.output_dir
n_thread = args.n_thread
# prefix = args.prefix

import os
import pysam
from collections import defaultdict
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd 
def reverse_tuple(forward_list):
    readlen = forward_list[-1][1]
    reverse_list = []
    for x in forward_list[::-1]:
        a,b,c,d = x
        reverse_list.append((readlen - b, readlen -a , c,d))
    return reverse_list
    
def collect_inv_qname(bam_file,  search_region1,search_region2 ):

    # Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")

    rn_strand_dc = defaultdict(list)
    for read in samfile.fetch( *search_region1):
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+" 
        rn_strand_dc[read.qname].append( strand)
    for read in samfile.fetch( *search_region2):
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+" 
        rn_strand_dc[read.qname].append( strand)

    inv_qnames = []
    for qn in rn_strand_dc:
        if len(set(rn_strand_dc[qn]))>1:
            inv_qnames.append(qn)
    
    inv_qnames = set(inv_qnames)
    samfile.close()
    # print(len(rn_strand_dc), len(inv_qnames))
    return inv_qnames
    

def collect_bnd(bam_file,  search_region, inv_qnames ):
    chrom = search_region[0]
    breakpoints = []
    # Define regions for searching, considering resolution
    # Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # Search for split reads in the first region
    qid = 0
    for read in samfile.fetch( *search_region):
        if read.is_unmapped:
            continue
        if read.qname not in inv_qnames:
            continue
        cigartuples = read.cigartuples

        if (cigartuples[0][0] in [4,5]) or (cigartuples[-1][0] in [4,5]):
            qid+=1
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"     
            readlen = read.infer_read_length()
            left_clip_len = 0
            right_clip_len = 0
            if (cigartuples[0][0] in [4,5]) :  # left clipping
                left_clip_len = cigartuples[0][1]
            if (cigartuples[-1][0] in [4,5]) :  # right clipping
                right_clip_len = cigartuples[-1][1]

            aln_tuples = []

            if left_clip_len==0:
                match_len = readlen - right_clip_len
                aln_tuples.append((0, match_len,'match', match_len))
            else:
                aln_tuples.append((0, left_clip_len, 'clip', left_clip_len))
                match_len = readlen - left_clip_len - right_clip_len
                aln_tuples.append((left_clip_len , left_clip_len + match_len,'match', match_len))

            if right_clip_len!=0:
                aln_tuples.append((readlen - right_clip_len, readlen,'clip', right_clip_len ))
        
            # aln_tuples_rev = reverse_tuple(aln_tuples)

            if strand == "+":
                aln_use = aln_tuples
                match_start_ref = read.reference_start
                match_end_ref = read.reference_end
            else:
                aln_use = reverse_tuple(aln_tuples)
                match_start_ref = read.reference_end
                match_end_ref = read.reference_start
            breakpoints.append((qid, read.qname, readlen, left_clip_len, right_clip_len, chrom, match_start_ref, match_end_ref,strand,aln_use))
    samfile.close()
    return breakpoints

        
def organize_bnds(bnds):
    forward_dc = defaultdict(list)
    reverse_dc = defaultdict(list)

    for bnd in bnds:        
        qname = bnd[1]
        strand = bnd[8]
        # print(qname,"strand",strand)
        if strand == '+':
            forward_dc[qname].append(bnd)
        else:
            reverse_dc[qname].append(bnd)
    return forward_dc, reverse_dc 

def get_match_seq(aln_list):
    for x in aln_list:
        if x[2]=='match':
            return x[0],x[1]

def pair_bnds(bnds_fwd, bnds_rev, max_merge_dist):
    match_rg_list_fwd = [get_match_seq(bnd[-1]) for bnd in bnds_fwd]
    match_rg_list_rev = [get_match_seq(bnd[-1]) for bnd in bnds_rev]
    ref_rg_list_fwd = [  (bnd[6], bnd[7]) for bnd in bnds_fwd]
    ref_rg_list_rev = [  (bnd[6], bnd[7]) for bnd in bnds_rev]
    idx_list = []
    dist_list = []
    inv_bnd_list = []
    for i in range(len(bnds_fwd)):
        for j in range(len(bnds_rev)):
            start_fwd, end_fwd = match_rg_list_fwd[i]
            start_rev, end_rev = match_rg_list_rev[j]
            ref_start_fwd, ref_end_fwd = ref_rg_list_fwd[i]
            ref_start_rev, ref_end_rev = ref_rg_list_rev[j]
            if (start_fwd < start_rev):
                dist = abs(end_fwd - start_rev)
                if dist <= max_merge_dist:
                    # inv_bnd = int((ref_end_fwd+ref_start_rev)/2)
                    inv_bnd = tuple(sorted([ref_end_fwd,ref_start_rev]))
                    inv_bnd_list.append(inv_bnd)
                    dist_list.append(dist)
            else:
                dist = abs( end_rev - start_fwd)
                if dist <= max_merge_dist:
                    # inv_bnd = int((ref_end_rev+ref_start_fwd)/2)
                    inv_bnd = tuple(sorted([ref_end_rev,ref_start_fwd]))
                    inv_bnd_list.append(inv_bnd)
                    dist_list.append(dist)
            # idx_list.append((i,j))
    # best_k = dist_list.index(min(dist_list))
    # best_i, best_j = idx_list[best_k]
    # best_inv_bnd = inv_bnd_list[best_k]
    # return best_inv_bnd,min(dist_list),bnds_fwd[best_i], bnds_rev[best_j]
    return inv_bnd_list, dist_list


def reduce_cluster(bnd_list):
    pos1_list = [ x[0] for x in bnd_list]
    pos2_list = [ x[1] for x in bnd_list]
    pos1 = int(sum(pos1_list)/len(pos1_list))
    pos2 = int(sum(pos2_list)/len(pos2_list))
    new_bnd = (pos1, pos2,len(pos1_list))
    return new_bnd

def cluster_bnd(bnd_list, cluster_dist):
    if len(bnd_list)==0:
        return []
    
    clusters = [[bnd_list[0]]]
    if len(bnd_list)==1:
        return [(bnd_list[0][0],bnd_list[0][1],1)]
    else:
        for new_bnd in bnd_list[1:]:
            # print(new_bnd)
            merge_flag = 0
            for k in range(len(clusters)):
                per_cluster = clusters[k]
                sample_bnd = per_cluster[-1]
                pos1_new,pos2_new= new_bnd
                pos1_samp,pos2_samp = sample_bnd
                if  ( abs( pos1_new - pos1_samp) <= cluster_dist) \
                    & ( abs( pos2_new - pos2_samp) <= cluster_dist) :

                    clusters[k].append(new_bnd)
                    merge_flag = 1
                    # print(new_bnd, "added to cluster",k)
                    break
            if merge_flag == 0:
                clusters.append([new_bnd])
                # print(new_bnd, "generate new cluster",k+1)
    clusters_center = [   reduce_cluster(bnd_list) for bnd_list in clusters]
    return clusters_center

def process_bnds(bnds, max_merge_dist):
    forward_dc, reverse_dc = organize_bnds(bnds)
    inv_bnd_list = []
    info_list = []
    # print(forward_dc)
    # print(reverse_dc)
    for qname in forward_dc:
        bnds_fwd = forward_dc[qname]
        bnds_rev = reverse_dc[qname]
        # print("\n\n***",qname)
        # print("forward:", bnds_fwd)
        # print("reverse:",bnds_rev)
        inv_bnds,dist_list,  = pair_bnds(bnds_fwd, bnds_rev, max_merge_dist)
        info_list.extend(dist_list)
        inv_bnd_list.extend(inv_bnds)
    # final_invs = cluster_inv_bnd(inv_bnd_list, max_inv_dist)
    
    return  inv_bnd_list, info_list

def find_inv(bam_file,chromosome, rough_breakpoint_1, rough_breakpoint_2, 
             resolution,  cluster_dist,
             max_merge_dist, min_svlen ):
    assert rough_breakpoint_1 < rough_breakpoint_2
    search_region1 = (chromosome, 
                    max(0, rough_breakpoint_1 - resolution), 
                    rough_breakpoint_1 + resolution)
    search_region2 = (chromosome, 
                    max(0, rough_breakpoint_2 - resolution), 
                    rough_breakpoint_2 + resolution)
    inv_qnames = collect_inv_qname(bam_file, search_region1, search_region2 )
    bnds1 = collect_bnd(bam_file,search_region1,inv_qnames)
    bnds2 = collect_bnd(bam_file,search_region2,inv_qnames)



    bnds = bnds1 + bnds2

    # for bnd in bnds:
    #     print(bnd)
    inv_bnd_list, info_list= process_bnds(bnds,max_merge_dist)

    # final_invs1, inv_bnd_list1, info_list1 = process_bnds(bnds1, max_inv_dist)
    # final_invs2, inv_bnd_list2, info_list2 = process_bnds(bnds2, max_inv_dist)
    cluster_centers = cluster_bnd(inv_bnd_list, cluster_dist)

    final_bnd = []
    ## filter that is in the range
    for bnd in cluster_centers:
        if  ((rough_breakpoint_1 - resolution) < bnd[0] < (rough_breakpoint_1 + resolution)) \
        & ((rough_breakpoint_2 - resolution) < bnd[1] < (rough_breakpoint_2 + resolution)) \
            & ((bnd[1]-bnd[0])>= min_svlen):
            final_bnd.append(bnd)
    # print("cluster centers:")
    # print(cluster_centers)
    return final_bnd

def call_inv(bam_file, n_thread, cluster_dist1, cluster_dist2, max_merge_dist,
             min_svlen):
   sequences = Parallel(n_jobs=n_thread)(delayed(find_inv)\
                                       (bam_file, chrom1_list[i], pos1_list[i],
                                          pos2_list[i], resolution,
                                             cluster_dist1,
                                             max_merge_dist, min_svlen)
                                             for i in tqdm(range(len(chrom1_list))))

   # flatten result
   r = 0
   dc_inv_raw = defaultdict(list)
   for i in range(len(sequences)):
      chrom = chrom1_list[i]
      inv_list = []
      for y in sequences[i]:
         inv_list.append((y[0], y[1]))
      dc_inv_raw[chrom].extend(inv_list)
   n = 0
   for chrom in dc_inv_raw:
      r+= len(dc_inv_raw[chrom])
      dc_inv_raw[chrom] = cluster_bnd(dc_inv_raw[chrom], cluster_dist2)
      n += len(dc_inv_raw[chrom])
   print("raw call:", r," reduced call:",n)
   return dc_inv_raw

def eval_inv(dc_inv_raw, max_eval_dist):
    dc_inv = {}
    for x in dc_inv_raw:
        dc_inv[x] = dc_inv_raw[x].copy()
    # max_eval_dist = 500
    tp = 0
    n = 0
    for x in dc_inv:
        n+= len(dc_inv[x])

    for i in range(len(pos1_list)):
        chrom = chrom1_list[i]
        bnds1 = dc_inv[chrom]
        for bnd in bnds1:
            left_bnd = bnd[0]
            right_bnd = bnd[1]
            pos1_gt = pos1_list[i]
            pos2_gt = pos2_list[i]
            if (abs(left_bnd - pos1_gt) <= max_eval_dist)\
            & (abs(right_bnd - pos2_gt) <= max_eval_dist):
                tp+=1
                # print(bnd)
                bnds1.remove(bnd)
                break
    print('TP:',tp)
    recall = round(tp/len(pos1_list),2)
    print("recall:", recall)

    print("total call:",n)


def deep_copy(dc_inv_raw):
    dc_inv = {}
    for x in dc_inv_raw:
        dc_inv[x] = dc_inv_raw[x].copy()
    return dc_inv 

def count_dc(dc):
    n = 0
    for x in dc:
        n+=len(dc[x])
    return n 

def get_somatic(dc_tm, dc_bl, max_dist):
    dc_tm1 = deep_copy(dc_tm)
    dc_bl1 = deep_copy(dc_bl)
    n = 0
    print("(before) num of tumor call:", count_dc(dc_tm))
    for chrom in dc_tm1:
        bnds_tm = dc_tm1[chrom]
        if chrom in dc_bl1:
            bnds_bl = dc_bl1[chrom]
            keep_bnd_tm = []
            for bnd_tm in bnds_tm:
                match_flag = 0
                for bnd_bl in bnds_bl:
                    if (abs(bnd_bl[0] - bnd_tm[0])<= max_dist)\
                    & (abs(bnd_bl[1] - bnd_tm[1])<= max_dist):
                        match_flag = 1
                        bnds_bl.remove(bnd_bl)
                        n+=1
                        break
                if match_flag==0:
                    keep_bnd_tm.append(bnd_tm)

            dc_tm1[chrom] = keep_bnd_tm
    print("num match:",n)
    print("(after) num of tumor call:", count_dc(dc_tm1))
    return dc_tm1

def write_vcf(dc, outfile):
    with open(f"{code_dir}/header_hg38",'r') as f:
        header = f.readlines()
    '''#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample'''
    # outfile = "test.vcf"
    with open(outfile,'w') as f:
        f.writelines(header)
        n = 0
        for chrom in dc:
            bnds = dc[chrom]
            for bnd in bnds:
                n+=1
                svid = f'INV.{n}'
                line = f"{chrom}\t{bnd[0]}\t{svid}\tN\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN={bnd[1]-bnd[0]};END={bnd[1]}\t.\t.\n"
                f.write(line)
    return 

if __name__  == '__main__':
 
    import os
    code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

    # ------------- Load Data
    # df = pd.read_excel(f"{code_dir}/High_confidence_callset.xlsx")
    # df_inv = df[(df['SV_type'] == 'INV')].reset_index(drop = True)
    df_inv = pd.read_csv(bed_file, sep = '\t', header = None)
    # chrom1_list = df_inv['Chrom1'].tolist()
    # chrom2_list = df_inv['Chrom2'].tolist()
    # pos1_list_raw = df_inv['Pos1'].tolist()
    # pos2_list_raw = df_inv['Pos2'].tolist()

    chrom1_list = df_inv.iloc[:,0]
    chrom2_list = df_inv.iloc[:,0]
    pos1_list_raw = df_inv.iloc[:,1]
    pos2_list_raw = df_inv.iloc[:,2]

    pos1_list = []
    pos2_list = []
    for i in range(len(pos1_list_raw)):
        pos1 = min(pos1_list_raw[i], pos2_list_raw[i] )
        pos2 = max(pos1_list_raw[i], pos2_list_raw[i] )
        pos1_list.append(pos1)
        pos2_list.append(pos2)
    print("Num of INVs:",len(pos1_list))
    # print("All chromosomes:", set(chrom1_list))

    # ------------Parameters
    # n_thread = 40
    # bam_file="/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_ONT/minimap2/HCC1395_ONT.bam"
    # bam_file_BL="/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395BL_ONT/minimap2/HCC1395BL_ONT.bam"
    cluster_dist1=100
    cluster_dist2=500
    # resolution = 50000  # 50kb
    resolution = 0 # target regions already have flanking
    max_merge_dist = 100 # higher value for better recall
    min_svlen = 30

    # ---------- call INV
    print("call  INV...")
    dc_inv_TM = call_inv(bam_TM,n_thread, cluster_dist1, cluster_dist2, max_merge_dist, min_svlen)
    # print("call Normal INV...")
    # dc_inv_BL = call_inv(bam_BL,n_thread, cluster_dist1, cluster_dist2, max_merge_dist, min_svlen)

    # ---------- get somatic INV
    # max_dist = 500
    # dc_inv_SM = get_somatic(dc_inv_TM, dc_inv_BL, max_dist)

    # ----------- Evaluation
    # max_eval_dist = 500
    # print("*** Tumor eval")
    # eval_inv(dc_inv_TM, max_eval_dist)
    # print("*** Normal eval")
    # eval_inv(dc_inv_BL, max_eval_dist)
    # print("*** Somatic eval")
    # eval_inv(dc_inv_SM, max_eval_dist)

    # ------------ write vcf
    if not os.path.exists(output_dir):
        os.system("mkdir -p " + output_dir)

    vcf_tm = output_dir+"/INV.vcf"
    # vcf_bl = output_dir+"/"+prefix+"_INV_normal.vcf"

    write_vcf(dc_inv_TM, vcf_tm )
    # write_vcf(dc_inv_BL, vcf_bl )