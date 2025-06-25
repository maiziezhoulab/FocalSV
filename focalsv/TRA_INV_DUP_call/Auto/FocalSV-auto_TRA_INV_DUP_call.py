from subprocess import Popen 

import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_file','-bam')
parser.add_argument('--out_dir','-o')
parser.add_argument('--data_type','-d', choices = ['Hifi','CLR','ONT'])
parser.add_argument('--patient', '-p',required=True, )
parser.add_argument('--state','-s', required=True, choices=['Tumor', 'Normal'])
parser.add_argument('--num_threads','-thread', type = int, default = 22 )

args = parser.parse_args()
input_path = args.bam_file
output_dir = args.out_dir
dtype = args.data_type
patient = args.patient
state = args.state
n_thread = args.num_threads

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")




import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'



cmd =f'''{code_dir}/define_region.py  \
    --input_path {input_path} \
    -o {output_dir} \
    -d {dtype} \
    -t {n_thread}
    '''
Popen(cmd, shell = True).wait()


cmd =f'''{code_dir}/process_dup.py  \
    -bed {output_dir}/DUPs.bed \
    -bam {input_path} \
    -p {patient} \
    -o {output_dir} \
    -d {dtype} \
    -s {state} \
    -t {n_thread}
    '''
Popen(cmd, shell = True).wait()


cmd =f'''{code_dir}/process_tra_inv.py  \
    -p {patient} \
    -o {output_dir} \
    -d {dtype} \
    -s {state} 
    '''
Popen(cmd, shell = True).wait()


#---------need to add merging script