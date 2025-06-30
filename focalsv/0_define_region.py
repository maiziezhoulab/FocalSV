
from tqdm import tqdm
import numpy as np
from collections import defaultdict
import json
import sys
import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
sys.path.append(f'{code_dir}/TRA_INV_DUP_call/Auto')
from estimate_coverage import estimate_bam_cov
from subprocess import Popen

import gzip

def smart_open(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def load_vcf_pg(vcffile):
	dc = defaultdict(list)
	with smart_open(vcffile) as f:
		for line in f:
			if line[0]!='#':
				sv = line.split()
				chrom =sv[0]
				pos = int(sv[1])
				dc[chrom].append(pos)
	return dc 

def call_sig(dtype,bamfile, sigdir, reference, chromosome,t):

	os.system("mkdir -p " + sigdir)

	real_sig = f"{sig_dir}/signatures"
	if os.path.exists(real_sig):
		os.system("rm -r "+ real_sig)

	if  dtype == "Hifi":
		para ='''--max_cluster_bias_INS      1000 \
		--diff_ratio_merging_INS    0.9 \
		--max_cluster_bias_DEL  1000 \
		--diff_ratio_merging_DEL    0.5 '''
	elif dtype == 'CLR':
		para = '''--max_cluster_bias_INS      100 \
		--diff_ratio_merging_INS    0.3 \
		--max_cluster_bias_DEL  200 \
		--diff_ratio_merging_DEL    0.5 '''
	else:
		para = '''--max_cluster_bias_INS      100 \
		--diff_ratio_merging_INS    0.3 \
		--max_cluster_bias_DEL  100 \
		--diff_ratio_merging_DEL    0.3 '''
	cmd = f'''python3 {code_dir}/5_post_processing/Reads_Based_Scan/Reads_Based_Scan.py \
	{bamfile} \
	{reference} \
	{sigdir}/reads_draft_variants.vcf \
	{sigdir} \
	{para} \
	-chr {chromosome} -t {t} --genotype --retain_work_dir '''
	print(cmd)
	Popen(cmd, shell = True).wait()




def extract_gt30(sig_dir):
	cmd = '''awk -F'\t' '$4 > 30 {print $2 "\t" $3 "\t" $4}' %s/DEL.sigs > %s/DEL_gt30.sigs'''%(sig_dir,sig_dir)
	Popen(cmd, shell = True).wait()

	cmd = '''awk -F'\t' '$4 > 30 {print $2 "\t" $3 "\t" $4}' %s/INS.sigs > %s/INS_gt30.sigs'''%(sig_dir,sig_dir)
	Popen(cmd, shell = True).wait()


def load_signature_file(file_path):
    """
    Load a signature file into a list of tuples (chrom, pos, svlen).
    
    Args:
        file_path (str): Path to the signature file (tab-delimited, no header).
    
    Returns:
        List[Tuple[str, int, int]]: A list of tuples representing (chrom, pos, svlen).
    """
    signature_list = []
    with open(file_path, 'r') as f:
        for line in f:
            chrom, pos, svlen = line.strip().split()
            signature_list.append((chrom, int(pos), int(svlen)))
    return signature_list



# Define clustering function
def cluster_svs(signature_list, distance_threshold):
	clustered_results = []
	cluster_id = 0
	prev_chrom = None
	prev_pos = None

	for chrom, pos, svlen in tqdm(signature_list):
		# Check if we should start a new cluster
		if prev_chrom != chrom or (prev_pos is not None and pos - prev_pos > distance_threshold):
			clustered_results.append([(chrom, pos, svlen)])
		else:
			clustered_results[-1].append((chrom, pos, svlen))
			
		prev_chrom, prev_pos = chrom, pos

	return clustered_results

def cluster_stats(clustered_list):
	cl_sizes = np.array([len(sigs) for sigs in clustered_list])
	cl_spans = np.array([sigs[-1][1] - sigs[0][1] for sigs in clustered_list])
	for q in np.arange(0,1,0.1):
		print(q, np.quantile(cl_sizes,q))
	for q in np.arange(0,1,0.1):
		print(q, np.quantile(cl_spans,q))

def cluster_pos(poss, dt, fl):
	assert dt >= 2*fl
	clusters = [[poss[0]]]
	# fl = 10e3
	# dt = 2 * fl
	for i in range(1, len(poss)):
		pos = poss[i]
		if (pos - clusters[-1][-1]) < dt:
			clusters[-1].append(pos)
		else:
			clusters.append([pos])
		
	cl_size = np.array([  len(cl) for cl in clusters])
	print("mean cl size:",cl_size.mean().round(), "max cl size:",cl_size.max())

	centroids = [  (cl[0]-fl, cl[-1]+fl) for cl in clusters]
	# centroids = [  (cl[0], cl[-1]) for cl in clusters]
	return centroids


def cluster_wgs(dc,dt,fl):
	b = 0
	a = 0
	total_size = 0
	# fl = 25e3
	# dt = 50e3
	for chrom, poss in dc.items():
		# print(chrom)
		# print(len(poss))
		a+= len(poss)
		dc[chrom] = cluster_pos(sorted(poss), dt, fl)
		size = sum([reg[1]- reg[0] for reg in dc[chrom]])
		total_size += size
		# print(len(dc[chrom]))
		b+= len(dc[chrom])
	print(f"clustering distance threshold {int(dt/1e3)}k")
	print(f"flanking size {int(fl/1e3)}k")
	print(f"before clustering {a} regions")
	print(f"after clustering {b} regions")
	print(f"after clustering region size:{round(total_size/1e6)}M out of 3200M human genome ( {round(total_size/3.2e7)}%)")
	return dc 

def reduce_cluster(clustered_list, min_sig, svtype):
	if svtype == 'INS':
		final_sv = [(sigs[0][0],sigs[0][1] , sigs[-1][1]) for sigs in clustered_list if len(sigs) >= min_sig]
	else:
		final_sv =   [(sigs[0][0],sigs[0][1] , sigs[-1][1]+max([sig[2] for sig in sigs ])) for sigs in clustered_list if len(sigs) >= min_sig]
	print(f"before reduce, {len(clustered_list)} SVs")
	print(f"after reduce, {len(final_sv)} SVs")

	dc = defaultdict(list)
	for chrom, start, end in final_sv:
		dc[chrom].append((start, end))

	return final_sv,dc

def recluster_regions(regions, dt):
	sorted_regions = sorted(regions, key=lambda x: (x[0], x[1]))
	clusters = [sorted_regions[0]]
	cluster_components = [[sorted_regions[0]]]

	for i in range(1, len(regions)):
		prev_start, prev_end = clusters[-1]
		cur_region = sorted_regions[i]
		cur_start, cur_end = cur_region
		olp = min(prev_end, cur_end) - max(prev_start, cur_start)
		if olp > -dt:
			cluster_components[-1].append(cur_region)
			clusters[-1] = (min(cur_start, prev_start), max(prev_end, cur_end))
		else:
			clusters.append(cur_region)
			cluster_components.append([cur_region])
	assert len(cluster_components) == len(clusters)

	for i in range(len(cluster_components)):
		regions = cluster_components[i]
		start = min([reg[0] for reg in regions])
		end = max([reg[1] for reg in regions])
		assert start == clusters[i][0]
		assert end == clusters[i][1]
	
	return clusters 

def recluster_wgs(dc,dt,fl):

	assert re_dt>= 2*fl
	prev_cnt = 0
	after_cnt = 0
	total_size=0
	for chrom, regions in dc.items():
		prev_cnt+=len(regions)
		new_regions = [ (reg[0]-fl, reg[1]+fl) for reg in recluster_regions(regions,dt)]
		total_size+=sum([ reg[1]-reg[0] for reg in new_regions])
		after_cnt +=len(new_regions)
		dc[chrom] = new_regions 
	print(f"before recluster {prev_cnt}")
	print(f"after recluster {after_cnt}")
	print(f"after clustering region size:{round(total_size/1e6)}M out of 3200M human genome ( {round(total_size/3.2e7)}%)")
	return dc

def load_one_kind(sig_dir, svtype,dt_fine, min_sig,re_dt,fl ):
	
	sig_file = f"{sig_dir}/{svtype}_gt30.sigs"
	signature_list = load_signature_file(sig_file)
	# Run the clustering algorithm
	clustered_list = cluster_svs(signature_list, dt_fine)
	# Print results
	print(len(clustered_list))
	# cluster_stats(clustered_list)
	final_sv,dc = reduce_cluster(clustered_list, min_sig, svtype)
	dc = recluster_wgs(dc,re_dt,fl)

	return dc

def merge_dict(dca,dcb,dt, fl):
	new_dc ={}
	cnt1 = 0
	cnt2 = 0
	total_size = 0
	valid_chroms = set(['chr'+str(i) for i in range(1,23)])
	for chrom in valid_chroms:
		regions = dca[chrom] + dcb[chrom]
		cnt1+=len(regions)
		new_dc[chrom] = [(reg[0]-fl, reg[1]+fl) for reg in recluster_regions(regions,dt)]
		total_size += sum([ reg[1] - reg[0] for reg in new_dc[chrom]])
		cnt2+=len(new_dc[chrom])
	print(f"before merging {cnt1}")
	print(f"after merging {cnt2}")
	print(f"after clustering region size:{round(total_size/1e6)}M out of 3200M human genome ( {round(total_size/3.2e7)}%)")
	return new_dc,total_size,cnt2
	

def eval_ins(dc_pos, dc_bed,dt_edge, dc_id, suffix):
	cnt_tp = 0
	dc = defaultdict(list)
	pos_cnt = 0
	cnt_bench = 0
	fn_ids = []
	for chrom, poss in dc_pos.items():
		cnt_bench+=len(poss)
		region_list = dc_bed[chrom]
		high_conf_poss = []
		start_k = 0
		pos_cnt+=len(poss)
		for pos in sorted(poss):
			high_conf = 0
			hit_cnt = 0
			for k in range(start_k, len(region_list)):
				region = region_list[k]
				if (region[0]+dt_edge) < pos < (region[1]-dt_edge):
					dc[(chrom,region[0], region[1])].append(pos)
					start_k = k
					high_conf = 1
					hit_cnt+=1
					break
			if high_conf:
				high_conf_poss.append(pos)
			else:
				id = dc_id[(chrom, pos)]
				fn_ids.append(id)

		dc_pos[chrom] = high_conf_poss
		cnt_tp+=len(high_conf_poss)
	print(f"covered {cnt_tp} out of {cnt_bench} bench SVs(recall {round(cnt_tp/cnt_bench*100,2)}%)")

	cnt = 0
	missed_dc_bed = defaultdict(list)
	tp_dc_bed = defaultdict(list)
	for chrom, regions in dc_bed.items():
		for region in regions:
			full_region = (chrom, region[0], region[1])  
			olp = len(dc[full_region])
			if olp==0:
				cnt+=1
				missed_dc_bed[chrom].append(region)
			else:
				tp_dc_bed[chrom].append(region)
	print(f"num region not overlapping with any gold SV: {cnt}")

	def write_bed(dc, ofile):
		with open(ofile,'w') as f:
			for chrom,k in dc.items():
				for reg in k:
					f.write(f"{chrom}\t{reg[0]}\t{reg[1]}\n")
	write_bed(missed_dc_bed, "fp_poss")
	write_bed(tp_dc_bed, "tp_poss")

	with open(f"FN_INS_{suffix}_merge",'w') as f:
		f.write('\n'.join(fn_ids)+'\n')
	return  cnt_tp, cnt_bench


def eval_del(dc_pos, dc_bed, dt_edge,dc_id, suffix):
	cnt_tp = 0
	dc = defaultdict(list)
	pos_cnt = 0
	cnt_bench = 0
	fn_ids = []
	for chrom, poss in dc_pos.items():
		cnt_bench+=len(poss)
		region_list = dc_bed[chrom]
		high_conf_poss = []
		start_k = 0
		pos_cnt+=len(poss)
		for pos in sorted(poss):
			high_conf = 0
			hit_cnt = 0
			for k in range(start_k, len(region_list)):
				region = region_list[k]
				if ((region[0]+dt_edge) < pos[0]) & (pos[1] < (region[1] - dt_edge)):
					dc[(chrom,region[0], region[1])].append(pos)
					start_k = k
					high_conf = 1
					hit_cnt+=1
					break
			if high_conf:
				high_conf_poss.append(pos)
			else:
				id = dc_id[(chrom, pos[0])]
				fn_ids.append(id)
		dc_pos[chrom] = high_conf_poss
		cnt_tp+=len(high_conf_poss)
	print(f"covered {cnt_tp} out of {cnt_bench} bench SVs(recall {round(cnt_tp/cnt_bench*100,2)}%)")

	cnt = 0
	missed_dc_bed = defaultdict(list)
	tp_dc_bed = defaultdict(list)
	for chrom, regions in dc_bed.items():
		for region in regions:
			full_region = (chrom, region[0], region[1])  
			olp = len(dc[full_region])
			if olp==0:
				cnt+=1
				missed_dc_bed[chrom].append(region)
			else:
				tp_dc_bed[chrom].append(region)
	print(f"num region not overlapping with any gold SV: {cnt}")

	def write_bed(dc, ofile):
		with open(ofile,'w') as f:
			for chrom,k in dc.items():
				for reg in k:
					f.write(f"{chrom}\t{reg[0]}\t{reg[1]}\n")
	write_bed(missed_dc_bed, "fp_poss")
	write_bed(tp_dc_bed, "tp_poss")
	with open(f"FN_DEL_{suffix}_merge",'w') as f:
		f.write('\n'.join(fn_ids)+'\n')
	return cnt_tp, cnt_bench

def load_vcf(vcffile):
	dc = defaultdict(list)
	dc_ins = defaultdict(list)
	dc_del = defaultdict(list)
	dc_id = {}

	with open(vcffile,'r') as f:
		for line in f:
			if line[0]!='#':
				sv = line.split()
				id = sv[2]
				chrom =sv[0]
				pos = int(sv[1])
				dc[chrom].append(pos)
				dc_id[(chrom,pos)] = id
				if 'SVTYPE=INS' in line:
					dc_ins[chrom].append(pos)
				else:
					svsize = int(sv[7].split('SVLEN=')[1].split(';')[0])
					dc_del[chrom].append((pos, pos+abs(svsize)))
	
	return dc,dc_ins,dc_del,dc_id

def write_summary(n_region, total_size,cnt_tp_ins, cnt_bench_ins,cnt_tp_del, cnt_bench_del):
	dc_summary = {
		"num_region": n_region,
		"total_size(MB)": round(total_size/1e6),
		"perct_of_wg": round(total_size/(3.2*1e9)*100,2),
		"INS_bench": cnt_bench_ins,
		"INS_TP": cnt_tp_ins,
		"INS_miss": cnt_bench_ins - cnt_tp_ins,
		"INS_recall": round(cnt_tp_ins/cnt_bench_ins*100,2),
		"DEL_bench": cnt_bench_del,
		"DEL_TP": cnt_tp_del,
		"DEL_miss": cnt_bench_del - cnt_tp_del,
		"DEL_recall": round(cnt_tp_del/cnt_bench_del*100,2),
		
	}
	print(dc_summary)

	# Save the dictionary as a JSON file with indentation
	with open(outfile, "w") as json_file:
		json.dump(dc_summary, json_file, indent=4)


def dict2bed(dc, bedfile):
	with open(bedfile,'w') as f:
		for i in range(1,23):
			chrom = 'chr'+str(i)
			regions = dc[chrom]
			sorted_regions = sorted(regions, key=lambda x: x[1])
			for reg in sorted_regions:
				f.write(f"{chrom}\t{int(reg[0])}\t{int(reg[1])}\n")




import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--sig_dir','-s')
parser.add_argument('--bam_file','-bam')
parser.add_argument('--ref_file','-r')
parser.add_argument('--prior_file','-p')
parser.add_argument('--out_dir','-o')
parser.add_argument('--data_type','-d', choices=['HIFI','CLR','ONT'])
parser.add_argument('--lib','-l')
parser.add_argument('--num_threads','-thread', type = int, default = 8)

args = parser.parse_args()
# sig_dir = args.sig_dir
bamfile = args.bam_file
reference = args.ref_file
out_dir = args.out_dir
prior_file = args.prior_file
dtype = args.data_type
suffix = args.lib
t = args.num_threads

if dtype == 'HIFI':
	dtype = 'Hifi'

if dtype == 'Hifi':
	#----------hifi
	# Set your distance threshold (e.g., 1000 bp)
	dt_fine = 500

	min_sig = 1
	re_dt = 15e3
	fl = 7e3
	dt_edge = 5e3 # for ins eval
	# suffix = f"HiFi_L{lib}"
	# sig_dir = f"/lio/lfs/maiziezhou_lab/maiziezhou_lab/CanLuo/long_reads_project/Variant_Caller_Result/cuteSV/hifi/Hifi_L{lib}/cuteSV-1.0.11"

elif dtype == 'CLR':
	#---------CLR
	# Set your distance threshold (e.g., 1000 bp)
	dt_fine = 200
	re_dt = 15e3
	fl = 7e3
	dt_edge = 5e3 

	# lib = 3
	# suffix = f"{dtype}_L{lib}"
	# cov_list = [65,89,29]

	r = 0.17
	estimate_bam_cov(bamfile, out_dir, t)
	cov_file = out_dir+"/mean_cov"
	with open(cov_file,'r') as f:
		mean_cov = eval(f.read())
	min_sig = int(r * mean_cov)
	# print(f"CLR L{lib} cov {cov_list[lib-1]}, min sig {min_sig}")
	# sig_dir = f"/lio/lfs/maiziezhou_lab/maiziezhou_lab/CanLuo/long_reads_project/Variant_Caller_Result/cuteSV/CLR/CLR_L{lib}/cuteSV-1.0.11"
elif dtype == 'ONT':
	#---------ONT
	# Set your distance threshold (e.g., 1000 bp)
	dt_fine = 500
	re_dt = 15e3
	fl = 7e3
	dt_edge = 5e3 
	# lib = 3
	# suffix = f"{dtype}_L{lib}"
	# cov_list = [57,47,51]
	r = 0.17
	estimate_bam_cov(bamfile, out_dir, t)
	cov_file = out_dir+"/mean_cov"
	with open(cov_file,'r') as f:
		mean_cov = eval(f.read())
	min_sig = int(r * mean_cov)
	# print(f"ONT L{lib} cov {cov_list[lib-1]}, min sig {min_sig}")
	# if lib==1:
	# 	use_lib = 2
	# elif lib==2:
	# 	use_lib=5
	# else:
	# 	use_lib = 6

	
	# sig_dir = f"/lio/lfs/maiziezhou_lab/maiziezhou_lab/CanLuo/long_reads_project/Variant_Caller_Result/cuteSV/ONT_L{use_lib}/"


else:
	print("dtype only support Hifi, CLR ONT")
	exit()


######### call individual signature
sig_dir = out_dir + "/read_signature"
call_sig(dtype, bamfile, sig_dir, reference, 'wgs',t)
extract_gt30(sig_dir)

outfile = f"{out_dir}/summary_{dtype}_{suffix}.json"
print(outfile)
# benchfile = "bench.vcf"
# pg_file = f"/data/maiziezhou_lab/Yichen/Projects/MARS_long_reads/Pangenome_paper_reproduction/NA24385_MARS_HifiTrio9_assemblies_hg19/result_genotyping_noSmallBub.vcf"
# dc_bench, dc_ins_bench, dc_del_bench,dc_id_bench = load_vcf(benchfile)
dc_pg = load_vcf_pg(prior_file)
dc_pg = cluster_wgs(dc_pg,re_dt,0)
dc_del = load_one_kind(sig_dir, 'DEL',dt_fine, min_sig,re_dt , fl = 0)
dc_ins = load_one_kind(sig_dir, 'INS',dt_fine, min_sig,re_dt, fl = 0 )
dc_del,_,_ = merge_dict(dc_del, dc_pg, re_dt, fl)
dc_ins,_,_ = merge_dict(dc_ins, dc_pg, re_dt, fl)
dc, total_size, n_region = merge_dict(dc_ins, dc_del, re_dt, 0)
# cnt_tp_ins, cnt_bench_ins = eval_ins(dc_ins_bench, dc_ins, dt_edge, dc_id_bench, suffix)
# cnt_tp_del, cnt_bench_del = eval_del(dc_del_bench, dc_del, dt_edge, dc_id_bench, suffix)
# eval_ins(dc_bench, dc, dt_edge)

# write_summary(n_region, total_size,cnt_tp_ins, cnt_bench_ins,cnt_tp_del, cnt_bench_del)
bedfile = f"{out_dir}/SV_Regions_{dtype}_{suffix}.bed"
dict2bed(dc, bedfile)
