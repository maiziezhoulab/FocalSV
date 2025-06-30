import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--bam_BL','-bl')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--output_dir','-o')
parser.add_argument('--bed_file','-bed')
parser.add_argument('--n_thread','-t', type = int, default = 10 )
# parser.add_argument('--prefix','-p', default="sample")
args = parser.parse_args()
# bam_BL = args.bam_BL
bam_TM = args.bamfile
output_dir = args.output_dir
n_thread = args.n_thread
bed_file = args.bed_file
# prefix = args.prefix

import pysam
from collections import defaultdict
import pandas as pd  
from joblib import Parallel, delayed
from tqdm import tqdm
import os


def reverse_tuple(forward_list):
    readlen = forward_list[-1][1]
    reverse_list = []
    for x in forward_list[::-1]:
        a,b,c,d = x
        reverse_list.append((readlen - b, readlen -a , c,d))
    return reverse_list
    

def collect_bnd(bam_file, chrom, start,end, resolution):
    breakpoints = []
    # Define regions for searching, considering resolution
    search_region = (max(0, start - resolution), end + resolution)
    # Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")
    # Search for split reads in the first region
    qid = 0
    for read in samfile.fetch(chrom, *search_region):
        if read.is_unmapped:
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
    

def get_match_seq(aln_list):
    for x in aln_list:
        if x[2]=='match':
            return x 
        
def verify_2_aln(aln1, aln2, thresh):
    match_tuple_1 = get_match_seq(aln1[-1])
    match_tuple_2 = get_match_seq(aln2[-1])
    chrom1, ref_start1, ref_end1, strand1, = aln1[5],aln1[-4], aln1[-3],aln1[-2]
    chrom2, ref_start2, ref_end2, strand2 = aln2[5],aln2[-4], aln2[-3],aln2[-2]
    start1, end1, _,_ = match_tuple_1
    start2, end2, _,_ = match_tuple_2
    if strand1 == strand2 :
        direction = "++"
    else:
        direction = "+-"

    if (abs(start1 - end2) <= thresh):
        return (chrom1,ref_start1, chrom2,ref_end2,direction)
    elif (abs(start2 - end1) <= thresh):
        return (chrom1,ref_end1, chrom2,ref_start2,direction)
    else:
        return -1

def reduce_cluster(bnd_list):
    pos1_list = [ x[1] for x in bnd_list]
    pos2_list = [ x[3] for x in bnd_list]
    pos1 = int(sum(pos1_list)/len(pos1_list))
    pos2 = int(sum(pos2_list)/len(pos2_list))
    new_bnd = list(bnd_list[0]).copy()
    new_bnd[1] = pos1 
    new_bnd[3] = pos2
    new_bnd.append(len(bnd_list))
    return new_bnd

def cluster_bnd(bnd_list, cluster_dist):
    if len(bnd_list)==0:
        return []
    clusters = [[bnd_list[0]]]
    if len(bnd_list)==1:
        return clusters
    else:
        for new_bnd in bnd_list[1:]:
            # print(new_bnd)
            merge_flag = 0
            for k in range(len(clusters)):
                per_cluster = clusters[k]
                sample_bnd = per_cluster[-1]
                chrom1_new,pos1_new,chrom2_new,pos2_new,direction_new = new_bnd
                chrom1_samp,pos1_samp,chrom2_samp,pos2_samp,direction_samp = sample_bnd

                
                if ((chrom1_samp == chrom1_new) \
                    & (chrom2_samp == chrom2_new) \
                    & (direction_new == direction_samp) \
                    & ( abs( pos1_new - pos1_samp) <= cluster_dist) \
                    & ( abs( pos2_new - pos2_samp) <= cluster_dist)) :
                    clusters[k].append(new_bnd)
                    merge_flag = 1
                    # print(new_bnd, "added to cluster",k)
                    break
            if merge_flag == 0:
                clusters.append([new_bnd])
                # print(new_bnd, "generate new cluster",k+1)
    
    clusters_center = [   reduce_cluster(bnd_list) for bnd_list in clusters]
    return clusters_center
    
def get_bnd(bam_file, n_thread):
    sequences = Parallel(n_jobs=n_thread)(delayed(infer_breakpoints)\
                                        (bam_file, chrom1_list[i], start1_list[i], end1_list[i],
                                            chrom2_list[i], start2_list[i], end2_list[i], resolution,
                                                thresh, cluster_dist)
                                                for i in tqdm(range(len(chrom1_list))))
    bnd_list = []
    for x in sequences:
        if x:
            bnd_list.extend(x)

    bnd_tuple_list = [ tuple(x[:-1]) for x in bnd_list if len(x)>1]
    new_bnd_list = cluster_bnd(bnd_tuple_list,cluster_dist = 2000)
    print(len(bnd_tuple_list),len(new_bnd_list))
    return new_bnd_list

                    



def infer_breakpoints(bam_file, chrom1, start1, end1,chrom2, start2,end2, resolution, thresh, cluster_dist):

    breakpoints1 = collect_bnd(bam_file, chrom1, start1,end1, resolution)
    breakpoints2 = collect_bnd(bam_file, chrom2, start2,end2, resolution)
    
    name1 = [x[1] for x in breakpoints1]
    name2 = [x[1] for x in breakpoints2]
    olp_names = set(name1) & set(name2)
    # print(olp_names)

    dc1 = defaultdict(list)
    dc2 = defaultdict(list)
    for x in breakpoints1:
        if x[1] in olp_names:
            dc1[x[1]].append(x)
    for x in breakpoints2:
        if x[1] in olp_names:
            dc2[x[1]].append(x)

    bnd_list = []
    support_list = []
    for rn in olp_names:

        for a in dc1[rn]:
            for b in dc2[rn]:
                infer_bnd = verify_2_aln(a,b,thresh)
                if infer_bnd!=-1:
                    bnd_list.append(infer_bnd)
                    support_list.append((a,b))
    
    # for x in bnd_list:
    #     print(x)
    # print(len(bnd_list))
    final_bnd_list = cluster_bnd(bnd_list, cluster_dist)
    # print(final_bnd_list)
    return final_bnd_list

def reformat_GT(df_tra):
    dc = defaultdict(list)
    chrom1_list = df_tra['Chrom1'].tolist()
    chrom2_list = df_tra['Chrom2'].tolist()
    pos1_list = df_tra['Pos1'].tolist()
    pos2_list = df_tra['Pos2'].tolist()
    for i in range(len(chrom1_list)):
        chrom1 = chrom1_list[i]
        chrom2 = chrom2_list[i]
        pos1 = pos1_list[i]
        pos2 = pos2_list[i]
        key = tuple(sorted([chrom1, chrom2]))
        if key == (chrom1, chrom2):
            dc[key].append((pos1,pos2))
        else:
            dc[key].append((pos2,pos1))
    return dc 

def reformat_final_bnd(new_bnd_list, cluster_dist):
    dc = defaultdict(list)
    for bnd in new_bnd_list:
        chrom1,pos1, chrom2, pos2,_,_ = bnd
        key = tuple(sorted([chrom1, chrom2]))
        if key == (chrom1, chrom2):
            dc[key].append((pos1,pos2))
        else:
            dc[key].append((pos2,pos1))
    a = 0
    b = 0
    for key in dc :
        bnd_list = dc[key]
        a += len(bnd_list)
        dc[key] = cluster_bnd2(bnd_list, cluster_dist)
        b+= len(dc[key])
    print("original:",a,"after clustering:",b)
    return dc 


def eval_bnd(bench_dc, call_dc, max_dist):
    tp = 0
    total_gt = 0
    total_call = sum([ len(call_dc[x]) for x in call_dc])
    tp_set = []
    all_gt = []

    for key in bench_dc:
        total_gt += len(bench_dc[key])
        if key in call_dc:
            bench_list = bench_dc[key]
            call_list = call_dc[key]
            for pos_bench in bench_list:
                for pos_call in call_list:
                    # print(pos_bench[0],pos_call[0])
                    # print(abs(pos_bench[0] - pos_call[0]))
                    if (abs(pos_bench[0] - pos_call[0]) <= max_dist)\
                    & (abs(pos_bench[1] - pos_call[1]) <= max_dist):
                        call_list.remove(pos_call)
                        tp+=1
                        tp_set.append((key[0],key[1],pos_bench[0],pos_bench[1]))
                        break
                all_gt.append((key[0],key[1],pos_bench[0],pos_bench[1]))
    
    fp = total_call - tp 
    fn = total_gt -  tp 
    recall = round(tp/total_gt,2)
    precision = round(tp/total_call,2)
    if tp == 0:
        f1 = 0
    else:
        f1 = round(2 * recall * precision / (recall + precision),2)
    print("total gt:", total_gt)
    print("total call:", total_call)
    print("tp:", tp)
    print("fp:",fp)
    print("fn:",fn)
    print("recall:",recall)
    print("precision:", precision)
    print("f1:",f1)
    return set(tp_set),set(all_gt)


def reduce_cluster2(bnd_list):
    pos1_list = [ x[0] for x in bnd_list]
    pos2_list = [ x[1] for x in bnd_list]
    pos1 = int(sum(pos1_list)/len(pos1_list))
    pos2 = int(sum(pos2_list)/len(pos2_list))
    new_bnd = (pos1, pos2)
    return new_bnd

def cluster_bnd2(bnd_list, cluster_dist):
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
    clusters_center = [   reduce_cluster2(bnd_list) for bnd_list in clusters]
    return clusters_center


def write_vcf(dc, outfile):
    with open(f"{code_dir}/header_hg38",'r') as f:
        header = f.readlines()
    '''#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample'''
    # outfile = "test.vcf"
    # print(outfile)
    with open(outfile,'w') as f:
        f.writelines(header)
        n = 0
        for chrom1,chrom2 in dc:
            # print(chrom1,chrom2)
            bnds = dc[(chrom1,chrom2)]
            for bnd in bnds:
                # print(bnd)
                n+=1
                svid = f'TRA.{n}'
                line = f"{chrom1}\t{bnd[0]}\t{svid}\tN\tN[{chrom2}:{bnd[1]}[\t.\tPASS\tSVTYPE=BND\t.\t.\n"
                f.write(line)
    return 

import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

# ------------- Load Data


# df = pd.read_excel(f"{code_dir}/High_confidence_callset.xlsx")
df_tra = pd.rad_csv(bed_file,sep='\t', header = None)
df_tra.columns = ['Chrom1','Start1','End1','Chrom2','Start2','End2','svtype']
# df_tra = df[df['SV_type'] == 'TRA'].reset_index(drop = True)
chrom1_list = df_tra['Chrom1'].tolist()
start1_list = df_tra['Start1'].tolist()
end1_list = df_tra['End1'].tolist()
chrom2_list = df_tra['Chrom2'].tolist()
start2_list = df_tra['Start2'].tolist()
end2_list = df_tra['End2'].tolist()

# --------------- Call TRA
# n_thread = 40
# resolution = 50000  # 10kb
resolution = 0  # 10kb
thresh = 1000
cluster_dist = 100

# bam_file="/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_ONT/minimap2/HCC1395_ONT.bam"
# bam_file_BL="/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395BL_ONT/minimap2/HCC1395BL_ONT.bam"
bnd_list = get_bnd(bam_TM, n_thread)
# bnd_list_BL = get_bnd(bam_BL, n_thread)

# ---------------- Evaluation
# dc_gt = reformat_GT(df_tra)
dc_tumor = reformat_final_bnd(bnd_list,cluster_dist)
# dc_normal = reformat_final_bnd(bnd_list_BL,cluster_dist)

# ----------------- Write VCF
if not os.path.exists(output_dir):
    os.system("mkdir -p " + output_dir)

vcf_tm = output_dir+"/TRA.vcf"
# vcf_bl = output_dir+"/"+prefix+"_TRA_normal.vcf"
write_vcf(dc_tumor, vcf_tm)
# write_vcf(dc_normal, vcf_bl)

# print("\n****eval tumor:")
# tp_tumor, gt = eval_bnd(dc_gt,dc_tumor, max_dist= 500)
# print("\n****eval normal:")
# tp_normal, gt = eval_bnd(dc_gt,dc_tumor, max_dist= 500)

