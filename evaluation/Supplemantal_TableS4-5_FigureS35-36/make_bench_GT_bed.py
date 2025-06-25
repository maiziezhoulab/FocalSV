import numpy as np
from collections import defaultdict
import pandas as pd
import pickle
import os
import subprocess
import argparse
import re
from argparse import ArgumentParser
import collections
from highlight_excel import highlight_excel_rows_by_extremes 
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--config','-config')
parser.add_argument('--input_dir','-i')
# parser.add_argument('--normal','-normal')
# parser.add_argument('--dtype','-d', choices= ['hifi', 'ont','clr'])
parser.add_argument('--out_dir','-o')
parser.add_argument('--sample','-s', default=['sample'])
parser.add_argument('--mode','-mode', choices=['strict','relax','prior','GT'])
parser.add_argument('--fp_thresh','-fpt', type = int, default = 1)

args = parser.parse_args()
# config_file = args.config
input_dir = args.input_dir
out_dir = args.out_dir
# dtype = args.dtype
sample = args.sample
global mode, fp_thresh
mode = args.mode
fp_thresh = args.fp_thresh
global use_chrs 
use_chrs = [ 'chr'+str(i) for i in range(1,23)]


def get_vcf(truvari_dir, out_dir, out_vcf,tool):
	cmd = f"zcat {truvari_dir}/INS_50_/tp-comp.vcf.gz |grep '#' > {out_dir}/header ;\
	zcat {truvari_dir}/INS_50_/tp-comp.vcf.gz {truvari_dir}/INS_50_/fp.vcf.gz {truvari_dir}/DEL_50_/tp-comp.vcf.gz {truvari_dir}/DEL_50_/fp.vcf.gz > {out_dir}/body;\
	cat {out_dir}/header {out_dir}/body | vcf-sort > {out_dir}/{tool}_raw.vcf"
	os.system(cmd)




def filter_vcf_by_bed(input_vcf, bed_file, output_vcf_gz):
	"""
	Compress a VCF with bgzip, index with tabix, and filter by a BED file using bcftools.
	Uses full shell commands instead of argument lists.
	"""

	# Step 1: bgzip if not already compressed
	if not input_vcf.endswith(".gz"):
		# compressed_vcf = input_vcf + ".gz"
		compressed_vcf = os.path.dirname(output_vcf_gz)+"/temp.vcf.gz"
		cmd_bgzip = f"bgzip -c {input_vcf} > {compressed_vcf}"
		subprocess.run(cmd_bgzip, shell=True, check=True)
	else:
		compressed_vcf = input_vcf

	# Step 2: tabix index
	cmd_tabix = f"tabix -p vcf {compressed_vcf}"
	subprocess.run(cmd_tabix, shell=True, check=True)

	# Step 3: bcftools view using BED file
	cmd_bcftools = f"ml  GCC/11.3.0 BCFtools/1.18;bcftools view -R {bed_file} {compressed_vcf} -Oz -o {output_vcf_gz}"
	subprocess.run(cmd_bcftools, shell=True, check=True)

	# Step 4: index the filtered output
	cmd_tabix_out = f"tabix -p vcf {output_vcf_gz}"
	subprocess.run(cmd_tabix_out, shell=True, check=True)
	


def preprocess(vcf,
               prefix,
               out_dir,
               chroms = None,
               minsize=50,
               svtypes=["ALL"],
               passonly=False,
			   bedfile = None):
	if ".vcf" not in vcf:
		get_vcf(vcf, out_dir, out_vcf,prefix)
		vcf = out_dir+f"/{prefix}_raw.vcf"



	out_vcf = out_dir+'/'+prefix+'.vcf'

	if chroms is None:
		chroms = ["chr"+str(i) for i in range(1,23)] # default for all autosomes

	fin = gzip.open(vcf, 'rt') if vcf.lower().endswith(".gz") else open(vcf,'r')

	# with open(vcf,'r') as fin:
	with open(out_vcf,'w') as fout:
		i = 0
		for line in fin:
			if line[0] == '#':
				fout.write(line)
			else:
				line = line.replace("SVLEN=>", "SVLEN=") #for NanoVar
				svinfo = line.split("\t")

				if svinfo[0] not in chroms:
					continue

				FILTER = svinfo[6]
				if passonly and FILTER != "." and FILTER != "PASS":
					continue

				try:
					svlen = float(re.findall("SVLEN=-?(\d+)",line)[0])
				except:
					svlen = abs(len(svinfo[3])-len(svinfo[4]))
				if svlen < minsize:
					continue

				if "ALL" in svtypes:
					i += 1
					svinfo[2] = prefix+"."+str(i)
					line = "\t".join(svinfo)
					fout.write(line)
				else:
					_svtype = "UNKNOWN" # TODO: make it compatitable with SV types other than DEL and INS 
					if "SVTYPE=" not in svinfo[7]:
						if "<" in svinfo[4] and ">" in svinfo[4]: # <DEL>
							_svtype = svinfo[4].strip("<>")
						elif (len(svinfo[3])-len(svinfo[4]))>0:
							_svtype = "DEL"
						elif (len(svinfo[3])-len(svinfo[4]))<0:
							_svtype = "INS"
					else:
						if 'SVTYPE=DEL' in svinfo[7]:
							_svtype = "DEL"
						elif 'SVTYPE=INS' in svinfo[7]:
							_svtype = "INS"
					if _svtype in svtypes:
						i += 1
						svinfo[2] = prefix +"."+str(i)
						line = "\t".join(svinfo)
						fout.write(line)
					else:
						continue

	fin.close()

	if bedfile is not None:
		out_vcf_gz = out_dir+'/'+prefix+'_bed.vcf.gz'
		filter_vcf_by_bed(out_vcf, bedfile, out_vcf_gz)

		return out_vcf_gz
	else:
		return out_vcf



def load_bed(bedfile,min_sup,min_mapq,min_size):
	dc = defaultdict(list)
	df = pd.read_csv(bedfile, sep = '\t', header = None)
	df['size'] = df.iloc[:,2] - df.iloc[:,1]
	
	df = df[(df.iloc[:,3]>=min_sup)      & (df.iloc[:,4]>=min_mapq)  & (df.iloc[:,5]>=min_size)].reset_index(drop=True)
	
	cnt = 0
	for i in range(df.shape[0]):
		chrom = df.iloc[i,0]
		start = df.iloc[i,1]
		end = df.iloc[i,2]
		size = end - start 
		sup = df.iloc[i,3]
		mapq = df.iloc[i,4]
		# print(mapq>=min_mapq)
		# print(size> 20e3)
		# if (mapq >= 60 ) & (( (size> 7e3) & (size< 10e3)& (sup >= 17) ) |  ( (size> 10e3) & (size< 200e3)& (sup >= 35) ) ):
		dc[chrom].append((start,end,))
		cnt+=1

	print(df.shape[0],cnt)
	return dc

def load_tra(bedfile,min_sup,min_mapq):
	dc = defaultdict(list)
	df = pd.read_csv(bedfile, sep = '\t', header = None)
	df = df[(df.iloc[:,5]>=min_sup) & (df.iloc[:,6]>=min_mapq)].reset_index(drop=True)
	
	cnt = 0
	for i in range(df.shape[0]):
		chrom1 = df.iloc[i,0]
		bnd1 = df.iloc[i,1]
		chrom2 = df.iloc[i,2]
		bnd2 = df.iloc[i,3]
		cnt+=1
		if int(chrom1[3:]) > int(chrom2[3:]):
			dc[(chrom2,chrom1)].append((bnd2,bnd1))
		else:
			dc[(chrom1,chrom2)].append((bnd1,bnd2))
			
	print(df.shape[0],cnt)
	return dc


def load_focalsv(indir,min_sup_bnd, min_mapq_bnd,min_sup_inv,min_mapq_inv,min_size_inv,min_sup_dup,min_mapq_dup,min_size_dup):

	bndfile = indir+"/TRA.tsv"
	invfile = indir+"/INVs.bed"
	dupfile = indir+"/DUPs.bed"

	dc_bnd = load_tra(bndfile, min_sup_bnd, min_mapq_bnd)
	dc_inv = load_bed(invfile, min_sup_inv, min_mapq_inv, min_size_inv)
	dc_dup = load_bed(dupfile, min_sup_dup, min_mapq_dup, min_size_dup)
	dc= {}
	dc['BND'] = dc_bnd 
	dc['INV'] = dc_inv 
	dc['DUP'] = dc_dup
	return dc 

import gzip

def open_vcf(path):
	"""
	Open a VCF (plain or gzipped) transparently.

	:param path: filepath ending in .vcf, .vcf.gz, or .vcf.bgz
	:param mode: file mode, e.g. 'rt' for read‐text or 'wt' for write‐text
	:return: a file‐like object
	"""
	opener = gzip.open if path.endswith(('.gz', '.bgz')) else open
	if path.endswith('gz'):
		mode = 'rt'
	else:
		mode = 'r'
	return opener(path, mode)


def load_vcf(vcffile, add_chr = False, out_prefix = None, tool = None):
	put_del = 0
	put_inv = 0
	put_dup = 0

	with open_vcf(vcffile) as f:
		dc = { v:defaultdict(list) for v in ['INS','DEL']}
		# dc_bnd = defaultdict(list)
		mate_ids = []
		for line in f:
			if line[0]=='#':
				continue
			data = line.split()
			filter = data[6]
			chrom1 = data[0]
			if add_chr:
				chrom1 = 'chr' + chrom1
			pos1 = int(data[1])
			try:
				svtype = data[7].split('SVTYPE=')[1].split(';')[0]
			except:
				if len(data[4]) > len(data[3]):
					svtype = 'INS'
				else:
					svtype = 'DEL'
			if 'MATEID' in data[7]:
				mateid = data[7].split('MATEID=')[1].split(';')[0]
				mate_ids.append(mateid)
			
			cur_id = data[2]
			if cur_id in set(mate_ids):
				# its mate has already been processed
				continue

			if  (chrom1 in use_chrs) & (svtype in ['DEL','INS']):
				data[2] = tool +'.'+data[2]
				line = '\t'.join(data)
				if svtype == 'INS':
					try:
						svlen = abs(int(data[7].split('SVLEN=')[1].split(';')[0]))
					except:
						svlen = abs(len(data[3])-len(data[4]))
					if svlen >= min_size :
						sv = (pos1, svlen, line)
						dc[svtype][chrom1].append(sv)
				else:
					if 'SVLEN' in data[7]:
						svlen = abs(int(data[7].split('SVLEN=')[1].split(';')[0]))
						sv = (pos1, pos1+svlen, line)
					elif 'END' in data[7]:
						end = abs(int(data[7].split('END=')[1].split(';')[0]))
						sv = (pos1, end, line)
					else:
						svlen = abs(len(data[3])-len(data[4]))
						sv = (pos1, pos1+svlen, line)


					if svlen >= min_size:
						dc[svtype][chrom1].append(sv)
					

	#----sanity check
	
	for key,dcv in dc.items():		
		cnt = 0
		bad_cnt = 0
		for chrom, vals in dcv.items():
			valid_vals = []
			if key not in ['INS', 'BND']:
				for val in vals:
					try:
						assert val[0]<val[1]
					except:
						print(val)
						bad_cnt+=1
						# exit()
					if val[0]<val[1]:
						valid_vals.append(val)
			else:
				valid_vals = vals
			dcv[chrom] = sorted(valid_vals, key = lambda x: x[0])
			cnt+=len(valid_vals)
		print(key, ': ', cnt )
		if bad_cnt:
			print(key,'(bad): ', bad_cnt)

	if out_prefix is not None:
		for key,dcv in dc.items():
			outfile = out_prefix+"_"+key 
			write_var(dcv, outfile)

	

	return dc 


def write_var_all_type(dc_all, out_prefix):

	for key,dcv in dc_all.items():
		outfile = out_prefix+"_"+key 
		write_var(dcv, outfile)

def match_sv(v_tm, v_nm, d,s, vtype):
	if vtype!='INS':
		return max( np.abs(np.array(v_tm[:2]) - np.array(v_nm[:2]) )) <= d
	else:
		sim = min(v_tm[1],v_nm[1])/max(v_tm[1],v_nm[1])
		return (abs(v_tm[0] - v_nm[0]) <= d ) & (sim>=s)


def derive_met(tp,fp,fn):
	bc_cnt = tp + fn 
	cp_cnt = tp+fp 
	rec = round(tp/bc_cnt,4) 
	prec = round(tp/cp_cnt,4)
	if tp:	
		f1 = round(2 * rec * prec / (rec + prec),4)
	else:
		f1 = 0
	return (bc_cnt, cp_cnt, tp, fp, fn, rec, prec, f1)



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
from collections import Counter

def most_voted_gt(gt_list, tie_value='./.'):
    count = Counter(gt_list)
    most_common = count.most_common()
    
    if not most_common:
        return tie_value  # Empty list case
    
    top_vote = most_common[0][1]
    top_candidates = [gt for gt, cnt in most_common if cnt == top_vote]
    
    if len(top_candidates) > 1:
        return tie_value
    else:
        return top_candidates[0]
	
def reduce_cluster(inv_list):
	# for inv in inv_list:
	#     print(inv)
	# exit()

	avg_start = int(np.mean([inv[0] for inv in inv_list]))
	avg_end = int(np.mean([inv[1] for inv in inv_list]))
	gt_list = []
	for inv in inv_list:
		gt = '/'.join(sorted(inv[3].split()[-1].split(':')[0].replace("|","/").split("/")))
		gt_list.append(gt)

	svid_list = []
	for inv in inv_list:
		svid = inv[3].split()[2]
		svid_list.append(svid)
	

	# avg_mapq = round(np.mean([inv[2] for inv in inv_list]),1)
	tools = [inv[2] for inv in inv_list]
	return avg_start, avg_end, len(inv_list), tools, ":".join(gt_list)+'\t'+':'.join(svid_list)

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




def compare_dict(dc_nm, dc_tm, d, s, eval_mode = False, target_vs = ['INS','DEL','DUP','INV','BND'], outdir= None):
	# dc_nm    -------> bench set
	# dc_tm    -------> comp set
	dc_somatic = { v:defaultdict(list) for v in target_vs}
	dc_eval = {}
	dc_eval_var = defaultdict(dict)
	for v_type in target_vs:
		dc_tm_v = dc_tm[v_type]
		dc_nm_v = dc_nm[v_type]
		tp = 0
		fp = 0
		fn = 0	

		tp_base = defaultdict(list)
		tp_call = defaultdict(list)
		fp_call = defaultdict(list)
		fn_base = defaultdict(list)
		for key in dc_tm_v:
			# print(len(dc_tm_v))
			v_list_tm = dc_tm_v[key]
			

			if key in dc_nm_v:
				# print('k')
				
				v_list_nm = dc_nm_v[key]
				start_i = 0
				used_is = []
				# print(len(dc_tm_v))
				
				for v_tm in v_list_tm:
					germline_flag = 0
					in_window_cnt = 0
					for i in range(start_i, len(v_list_nm)):
						v_nm = v_list_nm[i]
						# print(v_tm,v_nm)
						if (abs(v_tm[0] - v_nm[0]) <= d):
							# ----- in window
							if ( in_window_cnt < 1):
								start_i = i 
								in_window_cnt+=1
							if (match_sv(v_tm,v_nm,d,s ,v_type)) & (i not in used_is):
								used_is.append(i)
								germline_flag=1
								tp_call[key].append(v_tm)
								tp_base[key].append(v_nm)
								break
						if (v_nm[0] - v_tm[0]) > d:
							break

					if germline_flag==0:
						# somatic
						# print(v_tm)
						# exit()
						dc_somatic[v_type][key].append(v_tm)
						fp_call[key].append(v_tm)
						fp+=1

				
				tp += len(set(used_is))
				fn_is = set(list(range(len(v_list_nm)))) - set(used_is)
				for fn_i in fn_is:
					fn_base[key].append(v_list_nm[fn_i])
				fn += len(v_list_nm) - len(set(used_is))

			else:
				dc_somatic[v_type][key].extend(v_list_tm)
				fp_call[key].extend(v_list_tm)

				tp += 0
				fn += 0
				fp += len(v_list_tm) 
			# print(dc_tm_v.keys())
		
		for key in dc_nm_v:
			# print(len(dc_tm_v))
			if key not in dc_tm_v:
				fn += len(dc_nm_v[key])
				fn_base[key].extend(dc_nm_v[key])

		if eval_mode:
			dc_eval[v_type] = derive_met(tp,fp,fn)
			dc_eval_var[v_type]['TP-call'] = tp_call
			dc_eval_var[v_type]['TP-base'] = tp_base
			dc_eval_var[v_type]['FP-call'] = fp_call
			dc_eval_var[v_type]['FN-base'] = fn_base
		

	if not eval_mode:
		print("\n\nSomatic SV\n")
		for key,dc in dc_somatic.items():		
			cnt = 0
			for chrom, vals in dc.items():
				cnt+=len(vals)
			print(key, ': ', cnt )

		return dc_somatic
	else:
		print("\n\n Evaluation \n")
		print(dc_eval)
		os.system("mkdir -p " + outdir)
		df = pd.DataFrame(dc_eval)
		df.insert(0, 'Metric', ["bench_cnt", "comp_cnt", "TP", "FP", "FN", "Recall", "Precision", "F1"])
		df.to_csv(outdir+"/summary.csv", index = False)

		for vtype in dc_eval_var:
			for eval_type, dc_var in dc_eval_var[vtype].items():
				# print("\n***********",vtype, eval_type)
				eval_dir = outdir+"/"+vtype
				os.system("mkdir -p " + eval_dir)
				outfile = eval_dir+"/"+eval_type 
				write_var(dc_var, outfile)
		return dc_eval


def write_var(dc, outfile):
	with open(outfile,'w') as f:
		for key, vars in dc.items():
			for var in vars:
				# try:
				# 	print(var[0],var[1])
				# except:
				# 	print(var)
				# 	exit()
				f.write(f"{str(key)}\t{var[0]}\t{var[1]}\n")
	return 


def get_cnt(dc, v_list):
	cnt_list = []
	# ['INS','DEL','DUP','INV','BND']
	for v in v_list :		
		dcv = dc[v]
		cnt = 0
		for chrom, vals in dcv.items():
			cnt+=len(vals)
		cnt_list.append(cnt)
	return cnt_list

def write_vcf(dc, outfile):
	with open(outfile,'w') as f:
		for vtype in dc:
			dcv = dc[vtype]
			for chrom, vals in dcv.items():
				for val in vals:
					line = val[2]
					if ('SVTYPE=BND' in line) & (vtype!='BND'):
						data = line.split()
						info = data[7].split(';')

						for i in range(len(info)):
							x = info[i]
							if 'SVTYPE' in x:
								info[i] = 'SVTYPE='+vtype 
						data[7] = ';'.join(info)
						line = '\t'.join(data)+'\n'
					f.write(line)
	return


def get_summary_one_tool(dc_tm, dc_nm, dc_somatic, v_list):
	cnt_tm = get_cnt(dc_tm,v_list)	
	cnt_nm = get_cnt(dc_nm,v_list)	
	cnt_sm = get_cnt(dc_somatic ,v_list)	
	return cnt_tm+cnt_nm+cnt_sm				

# def load_focalsv(dtype):
# 	dc = {}
# 	indir = "/data/maiziezhou_lab/CanLuo/FocalSV/GR_Revision/ComplexSV/COLO829/eval"
# 	inv_dup_file = indir+"/summary_focalsv_inv_dup.csv"
# 	tra_file = indir+"/summary_focalsv_tra.csv"
# 	df1 = pd.read_csv(inv_dup_file)
# 	df2 = pd.read_csv(tra_file)
# 	dc['DUP'] = df1[f"{dtype}_DUP_focalsv"]
# 	dc['INV'] = df1[f"{dtype}_INV_focalsv"]
# 	dc['BND'] = df2[f"{dtype}_TRA_focalsv"]
# 	return dc

def write_eval_summary(dc_eval, tool_list, v_list, out_dir,dtype):
	df_list = []
	dc_focalsv = load_focalsv(dtype)
	for v in v_list:
		dc = {v: ['bench','call','TP','FP','FN','Recall','Precion','F1']}
		for tool in tool_list:
			vals = dc_eval[tool][v]
			dc[tool] = vals 
		df = pd.DataFrame(dc)
		if v in ['BND','INV','DUP']:
			df['FocalSV'] = dc_focalsv[v]
		df.to_csv(out_dir+f"/{tool}_eval.csv", index = False)
		df_list.append(df)

	# Write to Excel
	excel_file = out_dir+ '/All_tool_eval.xlsx'
	with pd.ExcelWriter(excel_file, engine='xlsxwriter') as writer:
		for df, name in zip(df_list, v_list):
			df.to_excel(writer, sheet_name=name, index=False)

	max_row_ids = [4,7,8,9]
	min_row_ids = [5,6]
	highlight_excel_rows_by_extremes(
		excel_file,
		max_row_ids,
		min_row_ids,)
		

def refomat_list(svs, tool):

	return [(sv[0],sv[1],tool,sv[2]) for sv in svs]

def reformat_dc_list(dc_list,v_list, tool_list):

	dc_new = {v: defaultdict(list) for v in v_list} # vtype -> chrom -> svs (start,end/size, tool)
	print(len(dc_list))
	print(len(tool_list))

	for i,dc_one_tool in enumerate(dc_list):
		tool = tool_list[i]
		for v in v_list:
			if v in dc_one_tool:
				dc_one_tool_var = dc_one_tool[v]
				for chrom, svs in dc_one_tool_var.items():
					dc_new[v][chrom].extend(refomat_list(svs, tool))
	return dc_new

def cluster_one_var(dc,d,n):
	new_dc = {}

	for chrom, svs in dc.items():
		clustered_svs = cluster_sig(svs,d,n)
		new_dc[chrom] = clustered_svs
	return new_dc


def cluster_all_var(dc,d,n):
	new_dc = {}
	for v,dcv in dc.items():
		dcv_clustered = cluster_one_var(dcv, d, n)
		new_dc[v] = dcv_clustered 

	return new_dc 



def get_stats(id_tool, id_bench):
	b = len(id_bench)
	c = len(id_tool)

	tp = len(id_tool & id_bench)
	fp = c - tp 
	fn = b - tp
	if b:
		r = tp/b
	else:
		r = 0
	p = tp/c
	if tp:
		f1 = 2 * r * p/(r+p)
	else:
		f1 = 0
	return b,c,tp,fp,fn,r,p,f1


def get_stats_GT(id_tool, id_tool_cor_gt,id_bench):
	b = len(id_bench)
	c = len(id_tool)

	tp = len(id_tool & id_bench & id_tool_cor_gt)
	fp = c - tp 
	fn = b - tp
	if b:
		r = tp/b
	else:
		r = 0
	p = tp/c
	if tp:
		f1 = 2 * r * p/(r+p)
	else:
		f1 = 0
	return b,c,tp,fp,fn,r,p,f1

def get_stats_prior(id_tool, id_bench,tooln ):
	b = len(id_bench)
	c = len(id_tool)

	tp = len(id_tool & id_bench)
	fp = c - tp 
	fn = b - tp
	if b:
		r = tp/b
	else:
		r = 0
	p_dict = {"FocalSV-auto":0.88,
	"FocalSV":0.88,
		   "PAV":0.86,
		   "SVIM-asm":0.86,
		   "Dipcall":0.73,
		   "sawfish":0.84,
		   "cuteSV":0.90,
		   "SVIM":0.64,
		   "PBSV":0.89,
		   "pbsv":0.89,
		   "Sniffles2":0.85,
		   "SKSV":0.8
	}
	p = tp/c
	p = p_dict[tooln]
	if tp:
		f1 = 2 * r * p/(r+p)
	else:
		f1 = 0
	return b,c,tp,fp,fn,r,p,f1


def get_stats_relax(id_tool, id_bench, id_bench_lowc):

	b = len(id_bench)
	c = len(id_tool)
	tp = len(id_tool & id_bench)
	fp = len(id_tool - id_bench_lowc)
	fn = b - tp
	if b:
		r = tp/b
	else:
		r = 0
	p = (c - fp)/c
	if tp:
		f1 = 2 * r * p/(r+p)
	else:
		f1 = 0
	return b,c,tp,fp,fn,r,p,f1


def check_GT(cols,):

	tools = cols[3]
	gts = cols[4].split('\t')[0].split(':')
	# svids = cols[5].split(':')
	if len(tools) != len(gts):
		print("len(tools) != len(gts)")
		exit()

	if len(tools)==1:
		return {tools[0]:'unknown'}

	tool_gt_map = dict(zip(tools, [gt.split(':')[0] for gt in gts]))
	# tool_svid_map = dict(zip(tools, [gt.split(':')[0] for gt in svids]))

	# majority vote
	all_gts = [gt for gt in tool_gt_map.values() if "." not in gt]
	try:
		majority_gt = collections.Counter(all_gts).most_common(1)[0][0]
	except:
		return {tool:'unknown' for tool in tools}

	dc = {}
	for tool, gt in tool_gt_map.items():
		# tool_total[tool] += 1
		# svid = tool_svid_map[tool]
		if gt == majority_gt:
			dc[tool] = 'correct'
		else:
			dc[tool] = 'wrong'
	return dc 


def get_one_var_stats(dc, min_tool_eval):
	idx = 0
	all_id = []
	bench_id = []
	bench_id_lowc = []
	dc_id_by_tool = defaultdict(list)
	dc_id_by_tool_cor_gt = defaultdict(list)
	id_sup_cnt_dc = {}
	for chrom, svs in dc.items():
		for sv in svs:
			idx+=1
			all_id.append(idx)
			cur_tools = sv[3]
			id_sup_cnt_dc[idx] = len(set(cur_tools))
			# print(cur_tools)
			if len(set(cur_tools))>= min_tool_eval:
				bench_id.append(idx)
			if len(set(cur_tools))>=  fp_thresh:
				bench_id_lowc.append(idx)

			for tool in cur_tools:
				dc_id_by_tool[tool].append(idx)
			
			dc_gt = check_GT(sv)
			for tool,gt_check in dc_gt.items():
				if gt_check == 'correct':
					dc_id_by_tool_cor_gt[tool].append(idx)



	bench_id = set(bench_id)
	bench_id_lowc = set(bench_id_lowc)
	print(len(all_id))
	print(len(bench_id))
	dc_eval = {'metric':['bench','comp','TP','FP','FN','Recall','Precision','F1','avg_support']}
	for tool in dc_id_by_tool:
		if mode  == "relax":
			stats = get_stats_relax(set(dc_id_by_tool[tool]), bench_id, bench_id_lowc)
		elif mode == 'strict':
			stats = get_stats(set(dc_id_by_tool[tool]), bench_id)
		elif mode == 'GT':
			stats = get_stats_GT(set(dc_id_by_tool[tool]), set(dc_id_by_tool_cor_gt[tool]),bench_id)
		else:
			stats = get_stats_prior(set(dc_id_by_tool[tool]), bench_id, tool)
		# print(tool,stats)
		mean_sup_cnt = np.mean([id_sup_cnt_dc[tool_svid] for tool_svid in set(dc_id_by_tool[tool])])
		dc_eval[tool] = list(stats)+[mean_sup_cnt]
	df = pd.DataFrame(dc_eval)
	return df

	
def get_all_var_stats(dc,suffix, min_tool_eval):
	dfs = []
	vs = []

	for v,dcv in dc.items():
		print(f"\n\n************{v}\n")
		df_one_var = get_one_var_stats(dcv, min_tool_eval)
		dfs.append(df_one_var)
		vs.append(v)

	# Write to Excel
	excel_file = out_dir+ f'/All_tool_eval_union_{sample}_{suffix}_{min_tool_eval}.xlsx'
	with pd.ExcelWriter(excel_file, engine='xlsxwriter') as writer:
		for df, name in zip(dfs, vs):
			df.to_excel(writer, sheet_name=name, index=False)

	max_row_ids = [4,7,8,9]
	min_row_ids = [5,6]
	highlight_excel_rows_by_extremes(
		excel_file,
		max_row_ids,
		min_row_ids,)


def get_all_var_stats_diff_support(dc, suffix, n_tool):

	for min_tool_eval in range(2, n_tool+1):
		print("min tools support:", min_tool_eval)
		get_all_var_stats(dc, suffix, min_tool_eval)

def get_high_conf_sv_one_var(dc):
	idx = 0
	all_id = []
	bench_id = []
	dc_id_by_tool = defaultdict(list)
	high_conf_sv_list = []
	dc_new = defaultdict(list)
	for chrom, svs in dc.items():
		for sv in svs:
			idx+=1
			all_id.append(idx)
			cur_tools = sv[3]
			gt = sv[4]
			# print(cur_tools)
			n_cur_tools = len(set(cur_tools))
			if n_cur_tools >= min_tool:
				bench_id.append(idx)
				high_conf_sv_list.append((chrom, sv, n_cur_tools,gt))
				dc_new[chrom].append(sv)
	return high_conf_sv_list, dc_new



def write_bnd(bnd_list, outfile, vtype):
	with open(outfile,'w') as f:
		for bnd in bnd_list:
			# print(bnd)
			if vtype == 'BND':
				chrom1, chrom2 = bnd[0]
				
				pos1,pos2, _ , tools = bnd[1]
				line = f"{chrom1}\t{chrom2}\t{pos1}\t{pos2}\t{';'.join(tools)}\n"
			else:
				chrom1  = bnd[0]
				pos1,pos2, _,tools, gt = bnd[1]
				line = f"{chrom1}\t{pos1}\t{pos2}\t{';'.join(tools)}\t{gt}\n"
			f.write(line)

	
def get_high_conf_sv_all_var(dc,suffix):
	dc_all_v = {}
	
	dc_cnt = {"data":[f"{sample}_{suffix}"]}
	for v in ['INS','DEL']:
		dcv = dc[v]
		print(f"\n\n************{v}\n")
		svs,dcv_bench = get_high_conf_sv_one_var(dcv)
		outfile = out_dir+ f'/Bench_{sample}_{suffix}_{v}.txt'
		write_bnd(svs, outfile, v)
		dc_all_v[v] = dcv_bench
		dc_cnt[v] = [len(svs)]
	df = pd.DataFrame(dc_cnt)
	df.to_csv(out_dir+ f'/Bench_{sample}_{suffix}_summary.csv', index = False)
	return dc_all_v


def write_somatic(dc,suffix):	
	dc_cnt = {"data":[f"{sample}_{suffix}"]}
	for v in ['INS','DEL','DUP','INV','BND']:
		dcv = dc[v]
		print(f"\n\n************{v}\n")
		svs = []
		for chrom,bnds in dcv.items():
			for bnd in bnds:
				svs.append((chrom, bnd))
		outfile = out_dir+ f'/Bench_{sample}_{suffix}_{v}.txt'
		write_bnd(svs, outfile, v)
		dc_cnt[v] = [len(svs)]
	df = pd.DataFrame(dc_cnt)
	df.to_csv(out_dir+ f'/Bench_{sample}_{suffix}_summary.csv', index = False)
	return 




def process_dict(dcs,d,n, v_list, tool_list, suffix):
	dc = reformat_dc_list(dcs, v_list, tool_list)

	dc_clustered =  cluster_all_var(dc,d,n)
	dc_bench = get_high_conf_sv_all_var(dc_clustered, suffix)
	get_all_var_stats_diff_support(dc_bench, suffix, len(tool_list))
	return dc_bench

	
def write_dict(dc, outfile):
	with open(outfile,'wb') as f:
		pickle.dump(dc,f, )

def load_pickle(infile):
	with open(infile,'rb') as f:
		obj = pickle.load(f)
	return obj

d=1000
s=0.3
global min_size 
min_size = 30
v_list = ['INS','DEL']

os.system("mkdir -p "+ out_dir)
# tools_tm = []
# with open(config_file,'r') as f:
# 	for line in f:
# 		tools_tm.append(line.strip().split(','))

tools_tm = [(vcf.split('_')[-1].split('.')[0], vcf) for vcf in glob.glob(f"{input_dir}/*.vcf")]


n = len(v_list)
dc = {'Source': ['Tumor'] * n + ['Normal'] * n + ['Somatic'] * n,
	  'SVtype': v_list + v_list + v_list }
# print('\n\n********** bench set \n')
# dc_bench = load_vcf(benchfile, add_chr=True)
dc_eval = {}
tool_list = []

dcs_tm = []
dcs_nm = []
dcs_sm = []
global min_tool
min_tool = 1

import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

bedfile = code_dir +"/High_Confidence_Region.hg38.bed"
# bedfile = None
vcf_folder = out_dir+"/VCFs"
os.system("mkdir -p " + vcf_folder)
processed_vcfs = []
for i, ( caller, vcf_tm) in enumerate(tools_tm):
	print(caller)
	# out_vcf = preprocess(vcf_tm,
    #            caller,
    #            vcf_folder,
    #            chroms = None,
    #            minsize=50,
    #            svtypes=["INS","DEL"],
	# 		   passonly=True,
	# 		   bedfile = bedfile)
	out_vcf = vcf_folder+f"/{caller}_bed.vcf.gz"
	filter_vcf_by_bed(vcf_tm,bedfile, out_vcf)
	processed_vcfs.append(out_vcf)

for i, ( caller, vcf_tm) in enumerate(tools_tm):
	tool_list.append(caller)
	print('\n\n**********',caller,'\n')
	tool_folder = out_dir+"/"+caller
	os.system("mkdir -p " + tool_folder)
	dc_tm = load_vcf(f"{vcf_folder}/{caller}_bed.vcf.gz",out_prefix= tool_folder+"/"+caller, tool = caller)
	dcs_tm.append(dc_tm)


write_dict(dcs_tm, out_dir + "/dcs_tm.pckl")

dcs_tm = load_pickle(out_dir + "/dcs_tm.pckl")

d=500
n=1

dc_tm = process_dict(dcs_tm,d,n, v_list, tool_list, sample)
write_dict(dc_tm, out_dir+f"/{sample}.pkl")



