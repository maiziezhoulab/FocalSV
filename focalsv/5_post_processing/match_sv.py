
from collections import defaultdict
import numpy as np
import gzip

def open_vcf(vcffile):
    if vcffile.endswith(".gz"):
        return gzip.open(vcffile, "rt")  # Open as text file for reading
    else:
        return open(vcffile, "r")  # Open as regular text file



def load_vcf(vcffile):
	min_svlen = 30
	max_svlen = 50e3
	dc = defaultdict(list)
	header = []
	with open_vcf(vcffile) as f:
		for line in f:
			if line[0]!='#':
				if ("SVTYPE=INS" in line) :
					data = line.split()
					chrom, pos , svinfo= data[0], int(data[1]), data[7]
					pos = int(pos)
					svtype = svinfo.split("SVTYPE=")[1].split(";")[0]
					# print(svinfo)
					svlen = abs(int(svinfo.split("SVLEN=")[1].split(";")[0]))
					if (max_svlen >= svlen >= min_svlen) &  (data[6]=='PASS'):
						dc[(chrom,svtype)].append([pos,svlen, line])
			else:
				header.append(line)
	return dc, header

def merge_dc(dc1, dc2):

	union_keys = set(list(dc1.keys()) + list(dc2.keys()) )
	dc3 = {}
	cnt = 0
	for key in union_keys:
		if (key in dc1) & (key in dc2):
			vals = sorted(dc1[key] + dc2[key] , key= lambda x: x[0])
		elif key in dc1:
			vals = dc1[key]
		else:
			vals = dc2[key]
		dc3[key] = vals 
		cnt+=len(vals)
	print("num of SV:", cnt)

	return dc3 


def match_2_list(sv_ref_list, sv_comp_list):
	dist_thresh = 200
	start_i = 0
	match_ref_sv_all = []
	for sv_comp in sv_comp_list:
		pos_comp = sv_comp[0]
		match_ref_sv = []
		updata_cnt = 0
		for i in range(start_i, len(sv_ref_list)):
			sv_ref = sv_ref_list[i]
			pos_ref = sv_ref[0]
			dist = abs(pos_comp-pos_ref)
			if dist <= dist_thresh:
				match_ref_sv.append((dist,sv_ref))
				if updata_cnt<1:
					start_i = i 
			elif (pos_comp - pos_ref) > dist_thresh:
				break 
		

			

def match_2_list(sv_ref_list, sv_comp_list):
	dist_thresh = 200
	start_i = 0
	match_comp_sv_all = []
	match_cnt_list = []
	no_match_cnt = 0
	
	for sv_ref in sv_ref_list:
		pos_ref = sv_ref[0]
		match_comp_sv = []
		updata_cnt = 0
		for i in range(start_i, len(sv_comp_list)):
			
			sv_comp = sv_comp_list[i]

			
			pos_comp = sv_comp[0]
			dist = abs(pos_comp-pos_ref)
			# print("dist:",dist)
			if dist <= dist_thresh:
				match_comp_sv.append((dist,sv_comp))
				# print(dist)
				if updata_cnt<1:
					start_i = i 
					# print(start_i)
				updata_cnt+=1
				
			elif (pos_comp - pos_ref) > dist_thresh:
				break 
		# print(start_i)
		match_comp_sv_all.append(match_comp_sv)
		if len(match_comp_sv)==0:
			no_match_cnt+=1
		else:
			match_cnt_list.append(len(match_comp_sv))
		# break
	# print(no_match_cnt, len(match_cnt_list),np.mean(match_cnt_list))
	# exit()
	return match_comp_sv_all

def pick_sv(sv_ref, sv_cluster_match):
	new_sv_list =[]
	no_match_cnt = 0
	for i in range(len(sv_cluster_match)):
		sv_cluster = sv_cluster_match[i]
		ref_gt = sv_ref[i][2].split('\t')[-1].split(':')[0]
		
		if len(sv_cluster)>1:
			opt_sv = sorted(sv_cluster, key=lambda x : x[1][1])[-1][1]
		elif len(sv_cluster)==1:
			opt_sv = sv_cluster[0][1]
		else:
			opt_sv = sv_ref[i]
			no_match_cnt+=1
		# update GT
		data = opt_sv[2].split('\t')
		data[-1] = ref_gt
		opt_sv[2] = '\t'.join(data)+'\n'
		new_sv_list.append(opt_sv)
	match_cnt = len(new_sv_list) - no_match_cnt
	new_sv_list = sorted(new_sv_list, key= lambda x:x[0])
	# print(new_sv_list[:2])
	return new_sv_list, match_cnt, no_match_cnt
	


def match_2_dc(dc_ref, dc_comp):

	lines = []
	match_cnt_all = 0
	no_match_cnt_all = 0
	for key in dc_ref:
		if key in dc_comp:
			# print(key)
			match_comp_sv = match_2_list(dc_ref[key], dc_comp[key])
			new_sv_list, match_cnt, no_match_cnt = pick_sv(dc_ref[key], match_comp_sv)
			match_cnt_all+= match_cnt
			no_match_cnt_all+= no_match_cnt
			for sv in new_sv_list:
				lines.append(sv[-1])
				# print(sv[-1])
				# exit()
	print("final variants:",len(lines))
	print("match_cnt_all:",match_cnt_all)
	print("no_match_cnt_all:",no_match_cnt_all)
	return lines

def write_vcf(header_ref, header_comp,lines, outfile):
	# print(lines[:2])
	with open(outfile,'w') as f:
		f.writelines(header_ref[:-1] + header_comp+lines)


def match_union_ins(comp_vcf, ref_vcf, outfile):
	dc_ref, header_ref = load_vcf(ref_vcf)
	dc_comp, header_comp = load_vcf(comp_vcf)
	lines = match_2_dc(dc_ref,dc_comp)
	write_vcf(header_ref, header_comp,lines, outfile)


# if __name__ == '__main__':
# 	match_union_ins(comp_vcf, ref_vcf, outfile)