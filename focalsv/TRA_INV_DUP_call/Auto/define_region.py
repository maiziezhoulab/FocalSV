
import pysam
from joblib import Parallel, delayed
from tqdm import tqdm
from collections import defaultdict
import pickle
import numpy as np
import pandas as pd
import os

def obj2pickle(my_dict, outfile):
    # Save dictionary to a pickle file
    with open(outfile, 'wb') as f:
        pickle.dump(my_dict, f)


def extract_relevant_reads_one_region(bamfile, region, min_mapq):
    bam = pysam.AlignmentFile(bamfile, "rb")
    inv_reads = []
    dup_reads = []
    tra_reads = []
    chrom, start, end =  region
    for read in bam.fetch(chrom, start, end, until_eof=True):

        if read.is_reverse:
            cur_strand = '-'
        else:
            cur_strand = '+'

        if read.has_tag("SA"):
            sa_tag = read.get_tag("SA")


            for sa in sa_tag.strip(';').split(';'):
                rname, pos, strand, cigar, mapq, nm = sa.split(',')
                mapq = int(mapq)
                read_dc = {
                'qname': read.query_name,
                'flag': read.flag,
                'rname': read.reference_name,
                'pos': read.reference_start,
                # 'cigar': read.cigarstring,
                'strand': '-' if read.is_reverse else '+',
                'is_supplementary': read.is_supplementary,
                'query_alignment_start': read.query_alignment_start,
                'query_alignment_end': read.query_alignment_end,
                'reference_end': read.reference_end,
                'read_len':read.infer_read_length(),
                'mapq':read.mapq
                # 'seq':read.seq
                    }

                if ( rname == read.reference_name ) & (read.mapq >= min_mapq) & (mapq>= min_mapq) :
                    inv_reads.append(read_dc)
                    if (strand == cur_strand) & ( rname == read.reference_name ):
                        dup_reads.append(read_dc)
                        break
                elif ( rname != read.reference_name ) & (read.mapq >= min_mapq) & (mapq>= min_mapq):
                    tra_reads.append(read_dc)
                    break
    bam.close()

    
    return inv_reads,dup_reads, tra_reads




def process_inv_reads(results):
    dc = defaultdict(list)

    for inv_reads in results:
        for read in inv_reads:
            dc[read['qname']].append(read)

    qnames = list(dc.keys())

    for qname in qnames:
        strands = [ read['strand'] for read in dc[qname]]
        if len(set(strands)) == 1:
            del dc[qname]

    return dc

def process_dup_reads(results):
    dc = defaultdict(list)

    for dup_reads in results:
        for read in dup_reads:
            dc[read['qname']].append(read)

    qnames = list(dc.keys())

    for qname in qnames:
        if len(dc[qname]) == 1:
            del dc[qname]
    return dc


def process_tra_reads(results):
    dc = defaultdict(list)

    for dup_reads in results:
        for read in dup_reads:
            dc[read['qname']].append(read)

    return dc

def extract_relevant_reads(bamfile, region_list, chrom, outdir, n_thread, min_mapq):

    results = Parallel(n_jobs=n_thread)(delayed(extract_relevant_reads_one_region)(bamfile, region,min_mapq) for region in tqdm(region_list, desc = f"extract relevant reads({chrom})"))
    
    inv_results = [ out[0] for out in results]
    dup_results = [ out[1] for out in results]
    tra_results = [ out[2] for out in results]
    dc_inv = process_inv_reads(inv_results)
    dc_dup = process_dup_reads(dup_results)
    dc_tra = process_tra_reads(tra_results)

    # inv_reads = []
    # dup_reads = []
    # for invs,dups in results:
    #     inv_reads.extend(invs)
    #     dup_reads.extend(dups)

    inv_outfile = outdir+f"/inv_relevant_reads_{chrom}.pkl"
    dup_outfile = outdir+f"/dup_relevant_reads_{chrom}.pkl"
    tra_outfile = outdir+f"/tra_relevant_reads_{chrom}.pkl"
    print('num inv reads:',len(dc_inv))
    print('num dup reads:',len(dc_dup))
    print('num tra reads:',len(dc_tra))
    
    obj2pickle(dc_inv, inv_outfile)
    obj2pickle(dc_dup, dup_outfile)
    obj2pickle(dc_tra, tra_outfile)

    return 

def get_chrom_blocks(bam_path, block_size, autusome_only):
    bam = pysam.AlignmentFile(bam_path, "rb")
    blocks = []
    autusome_chroms = [ 'chr'+str(i) for i in range(1,23)]

    dc = defaultdict(list)
    for chrom, length in zip(bam.references, bam.lengths):
        if autusome_only:
            if chrom not in autusome_chroms:
                continue

        for start in range(0, length, block_size):
            end = min(start + block_size, length)
            blocks.append((chrom, start, end))
            dc[chrom].append((chrom, start, end))
    print("block size:", block_size)
    print("num blocks:", len(blocks))

    return dc

def extract_coordinates(read):
    if read['strand'] == '+':
        start_cor_read_forward = read['query_alignment_start']
        end_cor_read_forward = read['query_alignment_end']
        start_cor_ref_forward = read['pos']
        end_cor_ref_forward = read['reference_end']
    else:
        read_len = read['read_len']
        start_cor_read_forward = read_len - read['query_alignment_end']
        end_cor_read_forward = read_len - read['query_alignment_start']
        start_cor_ref_forward = read['reference_end']
        end_cor_ref_forward = read['pos']
    return start_cor_read_forward,end_cor_read_forward,start_cor_ref_forward,end_cor_ref_forward


def process_a_pair_dup(read1, read2, dist_thresh_read, min_read_cov):
    read_len_R1 = read1['read_len']
    read_len_R2 = read2['read_len']
    assert read_len_R1==read_len_R2
    read_len = read_len_R1

    start_cor_read_forward_R1,end_cor_read_forward_R1,start_cor_ref_forward_R1,end_cor_ref_forward_R1 = extract_coordinates(read1)
    start_cor_read_forward_R2,end_cor_read_forward_R2,start_cor_ref_forward_R2,end_cor_ref_forward_R2 = extract_coordinates(read2)
    starts_list = [(start_cor_read_forward_R1,start_cor_ref_forward_R1),(start_cor_read_forward_R2,start_cor_ref_forward_R2)]
    ends_list = [(end_cor_read_forward_R1,end_cor_ref_forward_R1),(end_cor_read_forward_R2,end_cor_ref_forward_R2)]
    assert end_cor_read_forward_R1 > start_cor_read_forward_R1
    assert end_cor_read_forward_R2 > start_cor_read_forward_R2
    # get max start
    max_start_read, max_start_ref = sorted(starts_list, key = lambda x: x[0])[-1]
    # get min end
    min_end_read, min_end_ref = sorted(ends_list, key = lambda x: x[0])[0]
    # print(starts_list)
    # print(ends_list)
    # print(max_start_read, max_start_ref)
    # print(min_end_read, min_end_ref)
    # check distance on read
    cov_R1 = end_cor_read_forward_R1 -  start_cor_read_forward_R1
    cov_R2 = end_cor_read_forward_R2 -  start_cor_read_forward_R2
    cov = cov_R1 + cov_R2
    if (abs(max_start_read - min_end_read) <= dist_thresh_read) & ( cov/read_len >= min_read_cov ) :
        # is INV, return INV start end
        # print(sorted([max_start_ref, min_end_ref ]))
        return tuple(sorted([max_start_ref, min_end_ref,  ]) + [(read1['mapq'] + read2['mapq'])/2])
    else:
        # is not INV, return None
        return None


def process_a_pair_inv(read1, read2, dist_thresh_read, min_read_cov):
    read_len_R1 = read1['read_len']
    read_len_R2 = read2['read_len']
    assert read_len_R1==read_len_R2
    read_len = read_len_R1

    start_cor_read_forward_R1,end_cor_read_forward_R1,start_cor_ref_forward_R1,end_cor_ref_forward_R1 = extract_coordinates(read1)
    start_cor_read_forward_R2,end_cor_read_forward_R2,start_cor_ref_forward_R2,end_cor_ref_forward_R2 = extract_coordinates(read2)
    starts_list = [(start_cor_read_forward_R1,start_cor_ref_forward_R1),(start_cor_read_forward_R2,start_cor_ref_forward_R2)]
    ends_list = [(end_cor_read_forward_R1,end_cor_ref_forward_R1),(end_cor_read_forward_R2,end_cor_ref_forward_R2)]
    assert end_cor_read_forward_R1 > start_cor_read_forward_R1
    assert end_cor_read_forward_R2 > start_cor_read_forward_R2
    # get max start
    max_start_read, max_start_ref = sorted(starts_list, key = lambda x: x[0])[-1]
    # get min end
    min_end_read, min_end_ref = sorted(ends_list, key = lambda x: x[0])[0]
    # print(starts_list)
    # print(ends_list)
    # print(max_start_read, max_start_ref)
    # print(min_end_read, min_end_ref)
    # check distance on read
    cov_R1 = end_cor_read_forward_R1 -  start_cor_read_forward_R1
    cov_R2 = end_cor_read_forward_R2 -  start_cor_read_forward_R2
    cov = cov_R1 + cov_R2
    est_size = abs(max_start_ref - min_end_ref)
    if (abs(max_start_read - min_end_read) <= max(dist_thresh_read, est_size*0.15)) & ( cov/read_len >= min_read_cov ) :
        # is INV, return INV start end
        # print(sorted([max_start_ref, min_end_ref ]))
        return tuple(sorted([max_start_ref, min_end_ref,  ]) + [(read1['mapq'] + read2['mapq'])/2])
    else:
        # is not INV, return None
        return None
    

def process_a_pair_tra(read1, read2, dist_thresh_read, min_read_cov):
    chrom_R1 = read1['rname']
    chrom_R2 = read2['rname']
    assert chrom_R1!=chrom_R2
    strand_R1 = read1['strand']
    strand_R2 = read2['strand']
    read_len_R1 = read1['read_len']
    read_len_R2 = read2['read_len']
    assert read_len_R1==read_len_R2
    read_len = read_len_R1
    start_cor_read_forward_R1,end_cor_read_forward_R1,start_cor_ref_forward_R1,end_cor_ref_forward_R1 = extract_coordinates(read1)
    start_cor_read_forward_R2,end_cor_read_forward_R2,start_cor_ref_forward_R2,end_cor_ref_forward_R2 = extract_coordinates(read2)
    starts_list = [(start_cor_read_forward_R1,start_cor_ref_forward_R1, chrom_R1, strand_R1),(start_cor_read_forward_R2,start_cor_ref_forward_R2, chrom_R2, strand_R2)]
    ends_list = [(end_cor_read_forward_R1,end_cor_ref_forward_R1, chrom_R1, strand_R1),(end_cor_read_forward_R2,end_cor_ref_forward_R2, chrom_R2, strand_R2)]
    assert end_cor_read_forward_R1 > start_cor_read_forward_R1
    assert end_cor_read_forward_R2 > start_cor_read_forward_R2
    #--------- get the BND part, they are attached to each other
    # get max start
    max_start_read, max_start_ref, max_start_chrom, max_start_strand = sorted(starts_list, key = lambda x: x[0])[-1]
    # get min end
    min_end_read, min_end_ref, min_end_chrom, min_end_strand = sorted(ends_list, key = lambda x: x[0])[0]
    #--------- get the non-BND part, they are distant to each other
    # get min start
    min_start_read, min_start_ref, min_start_chrom, min_start_strand = sorted(starts_list, key = lambda x: x[0])[0]
    # get max end
    max_end_read, max_end_ref, max_end_chrom, max_end_strand = sorted(ends_list, key = lambda x: x[0])[-1]
    cov_R1 = end_cor_read_forward_R1 -  start_cor_read_forward_R1
    cov_R2 = end_cor_read_forward_R2 -  start_cor_read_forward_R2
    cov = cov_R1 + cov_R2


    if (max_start_chrom!=min_end_chrom)  & (abs(max_start_read - min_end_read) <= dist_thresh_read) & ( cov/read_len >= min_read_cov ):
        # one segment can't totally contain the other, 
        # bnd distance on read should be small
        # so it can be a TRA
        assert min_start_chrom == min_end_chrom
        assert max_start_chrom == max_end_chrom
        assert min_start_strand == min_end_strand
        assert max_start_strand == max_end_strand
        
        min_chrom, max_chrom  = min_start_chrom, max_start_chrom
        min_strand, max_strand = min_start_strand, max_start_strand

        #---------- get the segments coordinates
        # on read, the order is
        min_start_read, min_end_read, max_start_read, max_end_read
        # on reference, the order is
        min_start_ref, min_end_ref, max_start_ref, max_end_ref
        if (min_strand == '+') & (max_strand == '+'):
            assert min_start_ref<min_end_ref
            assert max_start_ref<max_end_ref
            if int(min_chrom[3:]) < int(max_chrom[3:]):
                bnd = [min_chrom, min_end_ref, f"N[{max_chrom}:{max_start_ref}["]
            else:
                bnd = [max_chrom, max_start_ref, f"]{min_chrom}:{min_end_ref}]N"]
            
        elif (min_strand == '-') & (max_strand == '-'):
            assert min_start_ref > min_end_ref
            assert max_start_ref > max_end_ref
            if int(min_chrom[3:]) < int(max_chrom[3:]):
                bnd = [min_chrom, min_end_ref, f"]{max_chrom}:{max_start_ref}]N"]
            else:
                bnd = [max_chrom, max_start_ref, f"N[{min_chrom}:{min_end_ref}["]

        elif (min_strand == '+') & (max_strand == '-'):
            assert min_start_ref < min_end_ref
            assert max_start_ref > max_end_ref
            if int(min_chrom[3:]) < int(max_chrom[3:]):
                bnd = [min_chrom, min_end_ref, f"N]{max_chrom}:{max_start_ref}]"]
            else:
                bnd = [max_chrom, max_start_ref, f"N]{min_chrom}:{min_end_ref}]"]
        else :
            assert min_start_ref > min_end_ref
            assert max_start_ref < max_end_ref
            
            if int(min_chrom[3:]) < int(max_chrom[3:]):
                bnd = [min_chrom, min_end_ref, f"[{max_chrom}:{max_start_ref}[N"]
            else:
                bnd = [max_chrom, max_start_ref, f"[{min_chrom}:{min_end_ref}[N"]

        return bnd + [(read1['mapq'] + read2['mapq'])/2 ]

    else:
        return None


def process_pair_one_chunk_inv(dc, qname_list, dist_thresh_read, min_read_cov):
    # ------------- need to be done by chromosome
    inv_list = []
    
    for qname in qname_list:
        #---------------INV
        forward_reads  = []
        reverse_reads = []
        for read in dc[qname]:
            if read['strand'] == '+':
                forward_reads.append(read)
            else:
                reverse_reads.append(read)

        for for_read in forward_reads:
            for rev_read in reverse_reads:
                result = process_a_pair_inv(for_read, rev_read, dist_thresh_read, min_read_cov)
                if result is not None:
                    inv_list.append(result)


    return inv_list


def process_pair_one_chunk_dup(dc, qname_list, dist_thresh_read,min_read_cov):
    # ------------- need to be done by chromosome
    dup_list = []
    for qname in qname_list:
        #----------------DUP
        reads = dc[qname]
        for i in range(len(reads)):
            for j in range(i+1,len(reads)):
                result = process_a_pair_dup(reads[i], reads[j], dist_thresh_read,min_read_cov)
                if result is not None:
                    dup_list.append(result)
                    # print(result, reads[i], reads[j])

    return dup_list


def process_pair_one_chunk_tra(dc, qname_list, dist_thresh_read, min_read_cov):
    # ------------- need to be done by whole genome
    tra_list = []
    for qname in qname_list:
        #---------------- TRA
        reads = dc[qname]
        for i in range(len(reads)):
            for j in range(i+1,len(reads)):
                if reads[i]['rname']!=reads[j]['rname']:
                    result = process_a_pair_tra(reads[i], reads[j], dist_thresh_read, min_read_cov)
                    if result is not None:
                        tra_list.append(result)
                        # print(result, reads[i], reads[j])

    return tra_list



def chunk_list(lst, chunk_size):
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def flatten_list(deep_list):
    # print(len(deep_list[0]))

    shallow_list = []

    for sub_list in deep_list:
        shallow_list.extend(sub_list)
    return shallow_list


def infer_complex_sv_one_chromosome(outdir, chrom, dist_thresh_read, chunk_size , n_thread, min_read_cov_inv, min_read_cov_dup):
    inv_file = outdir+f"/inv_relevant_reads_{chrom}.pkl"
    dup_file = outdir+f"/dup_relevant_reads_{chrom}.pkl"
    # tra_file = outdir+f"/tra_relevant_reads_{chrom}.pkl"

    with open(inv_file , 'rb') as f:
        dc_inv = pickle.load(f)
    with open(dup_file , 'rb') as f:
        dc_dup = pickle.load(f)
    # with open(tra_file , 'rb') as f:
    #     dc_tra = pickle.load(f)

    #---------INV
    qnames = list(dc_inv.keys())
    qname_chunks = chunk_list(qnames, chunk_size)
    results_inv = Parallel(n_jobs=n_thread)(delayed(process_pair_one_chunk_inv)
                                        (dc_inv, qname_chunk, dist_thresh_read, min_read_cov_inv, ) 
                                        for qname_chunk in tqdm(qname_chunks, desc = "infer INV from "+ chrom))
    #---------DUP
    qnames = list(dc_dup.keys())
    qname_chunks = chunk_list(qnames, chunk_size)
    results_dup = Parallel(n_jobs=n_thread)(delayed(process_pair_one_chunk_dup)
                                        (dc_dup, qname_chunk, dist_thresh_read, min_read_cov_dup) 
                                        for qname_chunk in tqdm(qname_chunks, desc = "infer DUP from "+ chrom))
    # print(len(results_dup))

    
    cand_inv_list_one_chrom = flatten_list(results_inv)
    cand_dup_list_one_chrom = flatten_list(results_dup)
    # print(len(cand_inv_list_one_chrom))
    # print(len(cand_dup_list_one_chrom))

    return cand_inv_list_one_chrom,cand_dup_list_one_chrom 



def infer_tra_wgs(outdir, dist_thresh_read, chunk_size , n_thread, min_read_cov):
    tra_file = outdir+f"/By_Chrom_Out/tra_relevant_reads.pkl"


    with open(tra_file , 'rb') as f:
        dc_tra = pickle.load(f)

    #---------TRA
    qnames = list(dc_tra.keys())
    qname_chunks = chunk_list(qnames, chunk_size)
    results_tra = Parallel(n_jobs=n_thread)(delayed(process_pair_one_chunk_tra)
                                        (dc_tra, qname_chunk, dist_thresh_read, min_read_cov) 
                                        for qname_chunk in tqdm(qname_chunks, desc = "infer TRA from whole genome "))
    
    cand_tra_list = flatten_list(results_tra)

    print(len(cand_tra_list))

    return cand_tra_list




def cluster_by_end(cand_inv_list, dist_thresh_clustering):
    # print("input",cand_inv_list )
    # sort by end
    cand_inv_list = sorted(cand_inv_list, key = lambda x : x[1])

    # cluster by end
    inv_cluster_list = [ [ cand_inv_list[0]]]
    for cand_inv in cand_inv_list[1:]:
        if abs(inv_cluster_list[-1][-1][1] - cand_inv[1]) <= dist_thresh_clustering:
            inv_cluster_list[-1].append(cand_inv)
        else:
            inv_cluster_list.append([cand_inv])
    
    return inv_cluster_list

def reduce_cluster(inv_list):
    # for inv in inv_list:
    #     print(inv)
    # exit()
    starts = [inv[0] for inv in inv_list]
    ends = [inv[1] for inv in inv_list]
    avg_start = int(np.mean(starts))
    avg_end = int(np.mean(ends))
    avg_mapq = round(np.mean([inv[2] for inv in inv_list]),1)
    return avg_start, avg_end, len(inv_list), avg_mapq, round(np.std(starts),4), round(np.std(ends),4)

def cluster_sig(cand_inv_list, dist_thresh_clustering, min_sig_support):
    # sort by start
    cand_inv_list = sorted(cand_inv_list, key = lambda x : x[0])

    # cluster by start
    inv_cluster_list = [ [ cand_inv_list[0]]]
    for cand_inv in cand_inv_list[1:]:
        if abs(inv_cluster_list[-1][-1][0] - cand_inv[0]) <= dist_thresh_clustering:
            inv_cluster_list[-1].append(cand_inv)
        else:
            inv_cluster_list.append([cand_inv])
    
    # print(inv_cluster_list[:10])
    # exit()


    # cluster by end
    new_inv_cluster_list = []
    for inv_cluster in inv_cluster_list:
        new_inv_cluster_list.extend(cluster_by_end(inv_cluster,dist_thresh_clustering))
    # print(new_inv_cluster_list[:10])
    # reduce cluster
    final_invs =  [reduce_cluster(invs) for invs in new_inv_cluster_list if len(invs) >= min_sig_support]
    return final_invs

def write_bed(sv_list,chrom, outfile):
    all_svs = []
    for sv in sv_list:
        # all_svs.append((chrom,sv[0],sv[1], sv[2], sv[3]))
        all_svs.append([chrom] + list(sv) )
    df = pd.DataFrame(all_svs)
    df.to_csv(outfile, sep = '\t', index = False, header= False)

def write_bed_wgs(sv_list, outfile):
    df = pd.DataFrame(sv_list)
    df.to_csv(outfile, sep = '\t', index = False, header= False)

def process_one_chromosome(bamfile, outdir, blocks_one_chrom, chrom, n_thread, dist_thresh_read, dist_thresh_clustering, min_sig_support,chunk_size, min_mapq, min_read_cov_inv, min_read_cov_dup):
    
    #---------- extrat reads relevant to INV DUP TRA
    extract_relevant_reads(bamfile, blocks_one_chrom, chrom, outdir, n_thread, min_mapq )

    cand_inv_list, cand_dup_list = infer_complex_sv_one_chromosome(outdir, chrom, dist_thresh_read, chunk_size , n_thread, min_read_cov_inv, min_read_cov_dup)
    # print(cand_dup_list[:10])

    final_inv_list = cluster_sig(cand_inv_list, dist_thresh_clustering= 3000, min_sig_support= 1)
    final_dup_list = cluster_sig(cand_dup_list, dist_thresh_clustering , min_sig_support)

    bed_dup = outdir + f"/DUP_{chrom}.bed"
    bed_inv = outdir + f"/INV_{chrom}.bed"

    write_bed(final_inv_list, chrom, bed_inv)
    write_bed(final_dup_list, chrom,bed_dup)

    print("num INVs:", len(final_inv_list))
    print("num DUPs:", len(final_dup_list))
    return final_inv_list ,final_dup_list 

def format_tra(tra, key):

    bnd1, bnd2, cnt, mapq, std1, std2 = tra

    chrom1, chrom2, direction = key 
    bnd2_str = f"{chrom2}:{bnd2}"
    if direction == 'N[':
        alt = direction+bnd2_str+'['
    elif direction == 'N]':
        alt = direction+bnd2_str+']'
    elif direction == '[N':
        alt = '['+bnd2_str+direction
    else:
        alt = ']'+bnd2_str+direction

    line = f"{chrom1}\t{bnd1}\t{chrom2}\t{bnd2}\t{alt}\t{cnt}\t{mapq}\t{std1}\t{std2}\n"
    return line

def write_tra(dc, outfile):
    with open(outfile,'w') as f:
        for key in dc:
            tra_list = dc[key]
            for tra in tra_list:
                line = format_tra(tra, key)
                f.write(line)

def process_tra_wgs( outdir,  n_thread, dist_thresh_read, dist_thresh_clustering, min_sig_support,chunk_size, min_read_cov):
    


    cand_tra_list = infer_tra_wgs(outdir, dist_thresh_read, chunk_size , n_thread, min_read_cov)
    # print(cand_inv_list[:10])

    # organze tra list 
    dc = defaultdict(list)

    for tra in cand_tra_list:
        chrom1, bnd1, bnd2_str, mapq = tra 

        if 'N[' in bnd2_str:
            direction = 'N['
        elif 'N]' in bnd2_str:
            direction = 'N]'
        elif ']N' in bnd2_str:
            direction = ']N'
        else:
            direction = '[N'

        if '[' in bnd2_str:
            sep = '['
        else:
            sep = ']'
        
        chrom2, bnd2 = bnd2_str.split(sep)[1].split(':')
        bnd2 = int(bnd2)

        dc[(chrom1, chrom2, direction)].append((bnd1, bnd2, mapq))

    # cluster bnd
    cnt = 0
    dc_new = {}
    for key in dc:
        cand_tra_list = dc[key]
        final_tra_list = cluster_sig(cand_tra_list, dist_thresh_clustering, min_sig_support)
        dc_new[key] = final_tra_list
        cnt+=len(final_tra_list)

    tra_file = outdir + f"/TRA.tsv"
    write_tra(dc_new, tra_file)
    print("num TRAs:", cnt)
    return 


def merge_result_tra(outdir):
    tra_pickle_list = [outdir + f'/By_Chrom_Out/tra_relevant_reads_chr{i}.pkl' for i in range(1,23)]
    dc = defaultdict(list)
    for tra_pickle in tra_pickle_list:
        with open(tra_pickle, 'rb') as f:
            dc_cur = pickle.load(f) 
        for qname,reads in dc_cur.items():
            dc[qname].extend(reads)
    with open(outdir + '/By_Chrom_Out/tra_relevant_reads.pkl','wb' ) as f:
        pickle.dump(dc, f)
    print("Total TRA relevant reads: ", len(dc))

        


def infer_complex_sv(bamfile,outdir, dist_thresh_read, dist_thresh_clustering,min_sig_support, n_thread,min_mapq, min_read_cov_inv, min_read_cov_dup, min_read_cov_bnd):
    by_chrom_outdir = outdir+"/By_Chrom_Out"
    os.system("mkdir -p " + by_chrom_outdir )
    block_size = 1000000  # 1 Mb blocks
    chunk_size_paring = 1000
    dc_blocks = get_chrom_blocks(bamfile, block_size, autusome_only=True)
    all_invs = []
    all_dups = []
    #---------------------INV DUP
    print("\n****************** INV and DUP extraction \n")
    for chrom, blocks in dc_blocks.items():
        # if chrom != 'chr22':
        #     continue
        print("\n************** Process "+chrom+'\n')
        inv_chrom, dup_chrom = process_one_chromosome(bamfile, by_chrom_outdir, blocks, chrom, n_thread, dist_thresh_read, dist_thresh_clustering, min_sig_support,chunk_size_paring, min_mapq, min_read_cov_inv, min_read_cov_dup)

        for bnds in inv_chrom:
            # all_invs.append((chrom, bnds[0], bnds[1], bnds[2], bnds[3]))
            all_invs.append([chrom] + list(bnds))

        for bnds in dup_chrom:
            all_dups.append([chrom] + list(bnds))
            # all_dups.append((chrom, bnds[0], bnds[1], bnds[2], bnds[3]))

    print("\n\n************** INV and DUP Summary \n")
    print("Total INVs: ", len(all_invs))
    print("Total DUPs: ", len(all_dups))
    inv_bed = outdir + "/INVs.bed"
    dup_bed = outdir + "/DUPs.bed"
    write_bed_wgs(all_invs, inv_bed)
    write_bed_wgs(all_dups, dup_bed)
    #--------------------- TRA
    print("\n****************** TRA extraction \n")
    merge_result_tra(outdir)
    process_tra_wgs( outdir,  n_thread, dist_thresh_read, dist_thresh_clustering, min_sig_support,chunk_size_paring , min_read_cov_bnd)







import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--dtype','-d', choices = ['HIFI','CLR','ONT'])
parser.add_argument('--n_thread','-t', type = int, default = 22 )

args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
dtype = args.dtype
n_thread = args.n_thread

# Example usage
# bamfile = "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_Pacbio/minimap2/HCC1395_Pacbio.bam"
# outfile = 'clr_tm_inv.bed'
dist_thresh_read = 1000 
if dtype == 'HIFI':
    
    dist_thresh_clustering = 100
    min_read_cov_inv = 0.2
elif dtype =='CLR':
    dist_thresh_clustering = 200
    min_read_cov_inv = 0.2
else:
    dist_thresh_clustering = 300
    min_read_cov_inv = 0.7
min_sig_support = 1
min_mapq = 10

min_read_cov_dup = 0.9
min_read_cov_bnd = 0.9
# n_thread = 50
infer_complex_sv(input_path ,output_dir,  dist_thresh_read, dist_thresh_clustering,min_sig_support, n_thread, min_mapq, min_read_cov_inv, min_read_cov_dup, min_read_cov_bnd)




