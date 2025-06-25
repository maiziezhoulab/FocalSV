from collections import defaultdict
from tqdm import tqdm
def load_vcf(vcf_file):
	header = []
	dc = defaultdict(list)
	with open(vcf_file,'r') as f:
		for line in f:
			if line[0]=='#':
				header.append(line)
			elif '0/0' not in line:
				data = line.split()
				chrom, pos = data[0],data[1]
				svlen = abs(int(data[7].split('SVLEN=')[1].split(';')[0]))
				svtype = data[7].split('SVTYPE=')[1].split(';')[0]
				pos = int(pos)
				gt = data[-1].split(':')[0]
				
				dc[chrom].append([pos,svlen, gt,svtype, line])
	for chrom in dc:
		my_list = dc[chrom]

		sorted_list = sorted(my_list, key=lambda x: x[0])
		dc[chrom] = sorted_list
	return  dc,header
				
def size_similar(len1, len2, ):
    return min(len1, len2) / max(len1, len2) 


def pick_best_match_gt(cand_match_list):

	sorted_list = sorted(cand_match_list, key=lambda x: (-x[0], x[1]))
	# if len(sorted_list)>1:
	# 	print(sorted_list)
	# 	exit()

	return sorted_list[0][2]




def gt_impute_one_chromosome(vars_cand, vars_gt, dist_thresh, sim_thresh):
	start_i = 0
	# new_gt_list = []
	match_cnt = 0
	no_match_cnt = 0
	update = 0
	for var_cand in vars_cand:
		cand_match_list = []
		update_start = 0
		for i in range(start_i, len(vars_gt)):
			var_gt = vars_gt[i]
			d = var_cand[0] - var_gt[0]
			sim = size_similar(var_cand[1] , var_gt[1])
			if (abs(d) <= dist_thresh) & (update_start == 0):
				start_i = i
				update_start = 1
			if (abs(d) <= dist_thresh) & (sim >= sim_thresh ) & (var_cand[3] == var_gt[3]):					
				cand_match_list.append([sim,d,var_gt[2]])
			if var_gt[0] - var_cand[0] > dist_thresh:
				break

		if len(cand_match_list):
			new_gt = pick_best_match_gt(cand_match_list)
			if new_gt!=var_cand[2]:
				update+=1
			var_cand[2] = new_gt

			match_cnt+=1
		else:
			# new_gt = var_cand[2]
			no_match_cnt+=1
		
		# new_gt_list.append(new_gt)
	
	return vars_cand,match_cnt, no_match_cnt, update

def write_vcf(dc_cand, outfile, header):
	with open(outfile,'w') as f:
		f.writelines(header)
		for chrom in dc_cand:
			vars_cand = dc_cand[chrom]
			
			for var in vars_cand:
				data = var[4].split()
				data[-1] = var[2]
				line = '\t'.join(data)+'\n'
				f.write(line)


def gt_impute(vcf_cand, vcf_gt,outfile, dist_thresh, sim_thresh):
	dc_cand, header_cand = load_vcf(vcf_cand)
	dc_gt, header_gt = load_vcf(vcf_gt)
	total_match = 0
	total_unmatch = 0
	total_update = 0

	for chrom in tqdm(dc_cand):
		vars_cand = dc_cand[chrom]
		vars_gt = dc_gt[chrom]
		dc_cand[chrom],match_cnt, no_match_cnt, update = gt_impute_one_chromosome(vars_cand, vars_gt,dist_thresh, sim_thresh)
		# print(match_cnt, no_match_cnt)
		total_match+= match_cnt
		total_unmatch+= no_match_cnt
		total_update += update
	print("total match:",total_match)
	print("total unmatch:",total_unmatch)
	print("total update:", total_update)

	write_vcf(dc_cand, outfile, header_cand)




	



