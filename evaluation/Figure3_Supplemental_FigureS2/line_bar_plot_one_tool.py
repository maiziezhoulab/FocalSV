from argparse import ArgumentParser
import os
import numpy as np
import pandas as pd
import gzip
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--del_dir','-del')
parser.add_argument('--ins_dir','-ins')
parser.add_argument('--output_dir','-o')
parser.add_argument('--prefix',default = 'line_bar_data')
parser.add_argument('--delete_intermediate_file','-d', action='store_true')
args = parser.parse_args()
del_dir = args.del_dir + '/'
ins_dir = args.ins_dir + '/'
output_dir = args.output_dir
prefix = args.prefix
os.system("mkdir -p "+output_dir)

fn_path = 'fn.vcf.gz'
fp_path = 'fp.vcf.gz'
tp_path = 'tp-base.vcf.gz'


# positive_cuts = [100,150,200,250,300,350,400,500,750,1000,1500,2000,2500,3000,4000,5000,10000,50000]

positive_cuts = [50,100,200,400,800,1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,10e3,50e3]
# positive_cuts = [50,500,1000,1500,2000,2500,5000,50e3]

negative_cuts = list(np.array(positive_cuts)*-1)[::-1]
cuts = negative_cuts+positive_cuts

# cuts = []
lenth_ranges = []
for i in range(len(cuts)-1):
	if cuts[i]!=-50:
		lenth_ranges.append((int(cuts[i]),int(cuts[i+1])))
# print(lenth_ranges)

def parse_vcf(vcf_path,lenth_ranges):
	with gzip.open(vcf_path,'rt')as f:
		s = f.readlines()
	zero_list = [0]*len(lenth_ranges)
	dc = dict(zip(lenth_ranges,zero_list))

	for line in s:
		if line[0]!='#':
			data = line.split('\t') 
			svlen = len(data[4])-len(data[3])
			for len_rg in lenth_ranges:
				if svlen>0:
					if svlen>=len_rg[0] and svlen<len_rg[1]:
						dc[len_rg]+=1
						break
				else:
					if svlen>len_rg[0] and svlen<=len_rg[1]:
						dc[len_rg]+=1
						break
	return dc

def parse_pair_vcf(ins_vcf,del_vcf,lenth_ranges):
	dc_ins = parse_vcf(ins_vcf,lenth_ranges)
	dc_del = parse_vcf(del_vcf,lenth_ranges)
	dc = {}
	for rg in dc_del:
		dc[rg] = dc_ins[rg]+dc_del[rg]

	return dc



dc_fn = parse_pair_vcf(ins_dir+fn_path,del_dir+fn_path,lenth_ranges)
dc_fp = parse_pair_vcf(ins_dir+fp_path,del_dir+fp_path,lenth_ranges)
dc_tp = parse_pair_vcf(ins_dir+tp_path,del_dir+tp_path,lenth_ranges)


sv_num_dc = {}
f1_dc = {}
sv_num_list = []
f1_list = []
tp_list = []
fp_list = []
fn_list = []
for len_rg in dc_fn:
	tp = dc_tp[len_rg]
	fp = dc_fp[len_rg]
	fn = dc_fn[len_rg]
	tp_list.append(tp)
	fp_list.append(fp)
	fn_list.append(fn)
	sv_num = tp+fn
	total_call = tp+fp
	if total_call==0:
		precision = 0
	else:
		precision = tp/total_call
	if sv_num==0:
		recall = 0
	else:
		recall = tp/sv_num
	if (recall+precision)==0:
		f1 = 0
	else:
		f1 = 2*recall*precision/(recall+precision)
	sv_num_dc[len_rg]=sv_num
	f1_dc[len_rg]=f1
	sv_num_list.append(sv_num)
	f1_list.append(f1)


df = pd.DataFrame({"lenth_range":lenth_ranges,'SV_number':sv_num_list,'TP':tp_list,'FP':fp_list,'FN':fn_list,'F1':f1_list})
df['Recall'] = df['TP']/df['SV_number']
df['Precision'] = df['TP']/(df['TP']+df['FP'])
df.to_csv(output_dir+'/%s.csv'%prefix,index = False)















