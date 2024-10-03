from argparse import ArgumentParser
import os
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd
from subprocess import Popen
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--asm_csv','-asm')
parser.add_argument('--read_signature_dir','-rdsig')
parser.add_argument('--output_dir','-o')
parser.add_argument('--data_type','-dtype',help='CCS;CLR;ONT')
## optional
parser.add_argument('--n_thread','-t',type = int, default = 11)
parser.add_argument('--n_thread_align','-ta',type = int, default = 10)
parser.add_argument('--mem_per_thread','-mempt', default = '768M',help = "Set maximum memory per thread; suffix K/M/G recognized; default = 768M")

args = parser.parse_args()
asm_csv = args.asm_csv
read_signature_dir = args.read_signature_dir
output_dir = args.output_dir
data_type = args.data_type
## optional
n_thread = args.n_thread
n_thread_align = args.n_thread_align
mem_per_thread = args.mem_per_thread


code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
os.system("mkdir -p "+output_dir)

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


df = pd.read_csv(asm_csv)
fasta_temp = df['fasta_pattern'][0]
chr_place_holder = df['chrom_place_holder'][0]
fasta_list =  [fasta_temp.replace(chr_place_holder,str(i+1)) for i in range(22)]


def run_one_chrom(i):
	logger.info("process chromosome %d..."%(i+1))
	cmd = "python3 "+code_dir+"/DipPAV_variant_call.py \
	-contig %s -ref /data/maiziezhou_lab/CanLuo/long_reads_project/DipPAV2/hg19_ref_by_chr/hg19_chr%d.fa \
	-sigd %s -o %s -dtype %s -t %d -chr %d"%(
		fasta_list[i],
		i+1,
		read_signature_dir,
		output_dir+'/chr%d'%(i+1),
		data_type,
		n_thread_align,
		i+1
		)
	logger.info(cmd)
	Popen(cmd,shell = True).wait()
	logger.info("finish chromosome %d"%(i+1))
	return 


results = Parallel(n_jobs=n_thread)(delayed(run_one_chrom)(i) for i in tqdm(range(22)))



### collect all result to be one 

def read_vcf(vcf_path):
	with open(vcf_path,'r') as f:
		header = []
		body = []
		for line in f:
			if line[0]=='#':
				header.append(line)
			else:
				body.append(line)

	return header,body 

logger.info("collect all chromosome result")
body_wgs = []
for i in range(22):
	vcf_path = output_dir+'/chr%d/final_vcf/dippav_variant_no_redundancy.vcf'%(i+1)
	header,body = read_vcf(vcf_path)
	body_wgs.extend(body)

with open(output_dir+"/variants.vcf",'w') as f:
	f.writelines(header + body_wgs)












