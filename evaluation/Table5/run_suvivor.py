from argparse import ArgumentParser
from joblib import Parallel, delayed
from tqdm import tqdm
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--normalvcf','-nvcf', help = "dtype,tool,vcffile in each row")
parser.add_argument('--cancervcf','-cvcf', help = "dtype,tool,vcffile in each row")
parser.add_argument('--output_dir','-o')
parser.add_argument('--survivor_dir','-s')
parser.add_argument('--n_thread','-t',type = int)
parser.add_argument('--pass_only','-p', action='store_true')

args = parser.parse_args()
normalvcf = args.normalvcf
cancervcf = args.cancervcf
output_dir = args.output_dir
n_thread = args.n_thread
survivor_dir = args.survivor_dir
global pass_only
pass_only = args.pass_only
import logging
from subprocess import Popen
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

with open(normalvcf,'r') as f:
	nvcf_list = [line.strip().split(',') for line in f if  line[0]!='#']

with open(cancervcf,'r') as f:
	cvcf_list =  [line.strip().split(',') for line in f if  line[0]!='#']


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

os.system("mkdir -p "+output_dir )

def one_run(n_path,c_path,idx, survivor_dir):
	if pass_only:
		flag = ' 1 '
	else:
		flag = ' 0 '
	cmd = code_dir+"/get_somatic.sh "+\
	n_path+' '+c_path+' '+output_dir+'/%s/'%idx + flag + f" {survivor_dir}"
	Popen(cmd, shell=True).wait()
	return 

# for i in range(len(nvcf_list)):
# 	one_run(nvcf_list[i][2],cvcf_list[i][2],'_'.join(nvcf_list[i][:2]))

sequences = Parallel(n_jobs=n_thread)\
(delayed(one_run)(nvcf_list[i][2],cvcf_list[i][2],'_'.join(nvcf_list[i][:2]), survivor_dir)\
	for i in tqdm(range(len(nvcf_list))))










