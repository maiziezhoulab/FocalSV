# %%

import argparse
from argparse import ArgumentParser


parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_dir','-i')
parser.add_argument('--bam_file','-bam')
parser.add_argument('--bed_file','-bed')
parser.add_argument('--reference','-ref')
# parser.add_argument('--excel_file','-excel')
# parser.add_argument('--bed_file','-bed')
# parser.add_argument('--vcffile','-vcf')
parser.add_argument('--outdir','-o')
parser.add_argument('--datatype','-d', choices=['HIFI','CLR','ONT'])
parser.add_argument('--n_thread','-t', type = int, default = 80 )
parser.add_argument('--max_d_cluster','-dc', type = int, default = 100 )
parser.add_argument('--max_d_eval','-de', type = int, default = 100 )

# need heavy editing this file
args = parser.parse_args()
input_dir = args.input_dir
bam_file = args.bam_file 
outdir = args.outdir 
bed_file = args.bed_file
# vcffile = args.vcffile 
n_thread = args.n_thread
max_d_cluster = args.max_d_cluster
max_d_eval = args.max_d_eval
reference = args.reference
datatype = args.datatype

# bam_file="/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_ONT/minimap2/HCC1395_ONT.bam"

# vcffile = "/data/maiziezhou_lab/CanLuo/FocalSV/Result/Cancer_ONT/complex_sv/DUP/DUP_final.vcf"
# outdir = "./test"

# n_thread = 80
# max_d_cluster = 100
# max_d_eval = 100


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
import pysam
from collections import defaultdict,Counter
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import os 
from subprocess import Popen


def call_dup_from_contig(in_dir, bamfile, reference, datatype, out_dir, t):

    cmd = f'''python3 {code_dir}/call_DUP_from_contigs.py \
    --input_dir {in_dir} \
    --bamfile {bamfile} \
    --reference {reference} \
    --datatype {datatype} \
    --out_dir {out_dir} \
    --n_thread {t} '''
    Popen(cmd, shell = True).wait()



def retrive_read(qname, ort, chrom, pos, samfile):
    for read in samfile.fetch(chrom,pos, pos + 1):
        if read.is_reverse:
            new_ort = -1
        else:
            new_ort = 1
        if (read.qname == qname) & (read.pos == pos) & (new_ort == ort):
            return read 
    return -1
def get_read_info(read):
    ref_start,ref_end,start_cigar, end_cigar, readlen =\
        read.pos,read.reference_end,read.cigar[0],read.cigar[-1], read.infer_read_length()
    
    if start_cigar[0] in [4,5]:
        read_start = 0 + start_cigar[1]
    else:
        read_start = 0 
    
    if end_cigar[0] in [4,5]:
        read_end = readlen - end_cigar[1]
    else:
        read_end = readlen 

    assert read_end > read_start 

    return (read_start, read_end, ref_start, ref_end)

def pair_split_aln(aln1, aln2, max_seg_dist):
    if  aln1[0] < aln2[0]:
        left_seg = aln1
        right_seg = aln2
    else:
        left_seg = aln2
        right_seg = aln1 

    left_ref_start = left_seg[2] 
    right_ref_start = right_seg[2] 

    if (abs(left_seg[1] - right_seg[0]) <= max_seg_dist) & (left_ref_start > right_ref_start):
        # pass QA
        return tuple(sorted([left_seg[3],right_seg[2]]))
    else:
        # fail QA
        return -1


def cluster_dup(cand_dup_list,max_cluster_dist):
    if len(cand_dup_list)==0:
        return [],[]
    elif len(cand_dup_list)==1:
        return cand_dup_list, [1]
    # Sort the list of tuples by the first element
    cand_dup_list = sorted(cand_dup_list, key=lambda x: x[0])
    cluster_list = [[cand_dup_list[0]]]
    for i in range(1, len(cand_dup_list)):
        new_dup = cand_dup_list[i]
        old_dup = cluster_list[-1][-1]
        if (abs(new_dup[0] - old_dup[0]) <= max_cluster_dist) &  (abs(new_dup[1] - old_dup[1]) <= max_cluster_dist):
            
            cluster_list[-1].append(new_dup)
        else:
            cluster_list.append([new_dup])
    dup_list = []
    n_list = []
    for cluster in cluster_list:
        n = len(cluster)
        center = int(sum([dup[0] for dup in cluster])/n),int(sum([dup[1] for dup in cluster])/n)
        dup_list.append(center)
        n_list.append(n)
    return dup_list, n_list

    
def reduce_call(dc_call, max_cluster_dist):
    '''
    remove redundancy in DUPs
    '''
    dc_call_new = {}
    a = 0
    b = 0
    for chrom, dups in dc_call.items():
        a += len(dups)
        dc_call_new[chrom] = cluster_dup(dups, max_cluster_dist)[0]
        b+= len(dc_call_new[chrom])
    return dc_call_new 

def extract_reads(bam_file, chrom,start, end, min_clip):
    samfile = pysam.AlignmentFile(bam_file)

    samiter = samfile.fetch(chrom, start, end , min_clip)
    qnames_raw = []
    qn_pos = defaultdict(list)
    for read in samiter:
        if read.is_reverse:
            ort = -1
        else:
            ort = 1
        qnames_raw.append((read.qname, ort))
        qn_pos[(read.qname, ort)].append(read.pos)
    return qnames_raw, qn_pos

def call_dup(df_dup, bam_file,i,flank, min_clip):
    chrom = df_dup['chrom'][i]
    start = df_dup['start'][i]
    end = df_dup['end'][i]

    # flank = 50000
    max_seg_dist = 500
    max_cluster_dist = 1000


    qnames_raw, qn_pos = extract_reads(bam_file, chrom, start - flank, end + flank, min_clip)

    qnames_split = []
    for (qn,ort), cnt in Counter(qnames_raw).items():
        if cnt > 1:
            qnames_split.append((qn,ort))
    qnames_split = set(qnames_split)
    samfile = pysam.AlignmentFile(bam_file)
    qn_record = defaultdict(list)
    for read in samfile.fetch(chrom, start - flank, end+ flank):
        if read.is_reverse:
            ort = -1
        else:
            ort = 1
        qn = read.qname
        if (qn,ort) in qnames_split:
            qn_record[(qn,ort)].append(read)
    cand_dup_list = []
    for (qn,ort), reads in qn_record.items():
        # print("===========",qn,ort)
        aln_list = []
        for read in reads:
            aln_list.append(get_read_info(read))
        for i in range(len(aln_list)-1):
            for j in range (i+1, len(aln_list)):
                pair_result = pair_split_aln(aln_list[i], aln_list[j], max_seg_dist)
                if pair_result!= -1:
                    cand_dup_list.append(pair_result)

    dups, support_list = cluster_dup(cand_dup_list, max_cluster_dist)
    return dups, support_list, (chrom, start,end)


def benchmark(gt, call_list, max_eval_dist):
    for call in call_list:
        if (abs(call[0]-gt[0]) <=max_eval_dist) & (abs(call[1]-gt[1]) <=max_eval_dist):
            return 1
    return 0

def benchmark_all(dc_gt, dc_call, max_d):
    dc_call1 = {}
    for a,b in dc_call.items():
        dc_call1[a] = b.copy()

    tp = 0
    N = 0
    n = 0

    for chrom in dc_gt:
        N += len(dc_gt[chrom])
        if chrom in dc_call1:
            n += len(dc_call1[chrom])
            a = dc_gt[chrom]
            b = dc_call1[chrom]
            for gt in a:
                for call in b:
                    if ( abs(gt[0] - call[0]) <= max_d) & ( abs(gt[1] - call[1]) <= max_d):
                        tp +=1
                        b.remove(call)
    rec = round(tp/N,4)
    prec = round(tp/n,4)
    if (rec+prec) == 0:
        f1 = 0
    else:
        f1 = round ( 2 * rec * prec / (rec + prec),4)
    print("total bench:", N)
    print("total call:",n)
    print("tp:",tp)
    print("fp:", n - tp)
    print("fn:", N - tp)
    print("recall:", rec)
    print("precision:",prec)
    print("f1:",f1)





def load_gt(df_dup):
    dc_gt = defaultdict(list)
    for i in range(df_dup.shape[0]):
        chrom = df_dup['Chrom1'][i]
        start = df_dup['Pos1'][i]
        end = df_dup['Pos2'][i]
        dc_gt[chrom].append([start,end])
    return dc_gt

def write_vcf(vcffile, outdir, dc_final):
    if not os.path.exists(outdir):
        os.system("mkdir -p " + outdir)
    asm_vcf = outdir+"/asm.vcf"
    out_vcf = outdir + "/DUP.vcf"

    cmd = f'''sed "s/volcanosv/focalsv/g" {vcffile} > {asm_vcf}'''
    Popen(cmd, shell = True).wait()
    cnt = 0
    with open(asm_vcf, 'a') as f:
        for chrom in dc_final:
            dup_list = dc_final[chrom]
            for dup in dup_list:
                cnt +=1
                start, end = dup 
                line = f'{chrom}\t{start}\tfocalsv.DUP.aln.{cnt}\t.\t<DUP>\t20\tPASS\tSVTYPE=DUP;SVEND={end};SVLEN={end-start}\tGT\t0/1\n'
                f.write(line)

    cmd = f'''vcf-sort {asm_vcf} > {out_vcf};rm {asm_vcf}'''
    Popen(cmd, shell = True).wait()

def get_target_region_excel(excel_file):

    # excel_file = code_dir + "High_confidence_callset.xlsx"
    df = pd.read_excel(excel_file)

    # For faster demonstration, only used < 5M DUP
    # It is less reliable to when a read has split-alignment further than 5MB, so we only focus on DUPs <5M
    df_dup_target = df[(df['SV_type']=='DUP') & (df['SV_Size or breakpoints distance']<5000000)].reset_index(drop = True)
    df_dup = df[(df['SV_type']=='DUP') ].reset_index(drop = True)
    dc_gt = load_gt(df_dup)

    return df_dup_target, df_dup, dc_gt
    
def load_vcf(vcffile):
    # load FocalSV  DUP vcf
    dc_asm = defaultdict(list)
    with open(vcffile,'r') as f:
        for line in f:
            if line[0]!='#':
                data = line.split()
                chrom,pos = data[0], data[1]
                pos = int(pos)
                size = int(data[7].split("SVLEN=")[1].split(';')[0])
                end = pos + size
                dc_asm[chrom].append((pos, end))
    return dc_asm

def load_bed(bed_file):
    df = pd.read_csv(bed_file, sep = '\t', header = None)
    df.columns = ['chrom','start','end','svtype']
    return df 

contig_call_dir = outdir + "/dup_call_from_contig/"
vcffile = contig_call_dir+"/DUP/DUP_final.vcf"

## load bed 

if not os.path.exists(outdir):
    os.system("mkdir -p " + outdir)
df_dup_target = load_bed(bed_file)

print("-------------------call_dup_from_contig")
call_dup_from_contig(input_dir, bam_file, reference, datatype, contig_call_dir, n_thread)
# exit()
dc_asm = load_vcf(vcffile)

# %%
# evaluate somatic DUP from asm result
# print("\n*****Only asm result:")
# benchmark_all(dc_gt, dc_asm, max_d_eval)


print("Now start adding DUPs from alignment file...")

# print("Number of bench DUPs: ",df_dup.shape[0])
# print("Number of target DUPs(<5M): ",df_dup_target.shape[0])






# call dup from bam file
sequences = Parallel(n_jobs=n_thread)(delayed(call_dup)(df_dup_target,bam_file, i,flank=0, min_clip = 0) for i in tqdm(range(df_dup_target.shape[0])))

# collect dup into dict
dc_call = defaultdict(list)
for seq in sequences:
    chrom = seq[-1][0]
    dc_call[chrom].extend( seq[0])
        
# reduce dup call        
dc_call_new = reduce_call(dc_call, max_d_cluster)
# print("\n*****pure aln result(before reduce)")
# benchmark_all(dc_gt, dc_call, max_d_eval)
# print("\n*****pure aln result(after reduce)")
# benchmark_all(dc_gt, dc_call_new, max_d_eval)


# merge DUPs from alignment and asm
dc_new = defaultdict(list)

for chrom in dc_call_new:
    dc_new[chrom].extend(dc_call_new[chrom])

for chrom in dc_asm:
    dc_new[chrom].extend(dc_asm[chrom])

# reduce redundancy for merged calls
dc_new1 = reduce_call(dc_new, max_d_cluster)


# # evaluate merged result
# print("\n*****merged result:")
# benchmark_all(dc_gt, dc_new, max_d_eval)


# print("\n*****merged no redun result:")
# # evaluate merged call without redundancy
# benchmark_all(dc_gt, dc_new1, max_d_eval)


print("We use raw merged result as final output")

write_vcf(vcffile, outdir, dc_new)