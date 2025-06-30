import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
parser.add_argument('--hp','-hp', type = int, choices=[1,2] )

args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
hp = args.hp

cnt=0
with open(input_path,'r') as f:
	with open(output_path,'w') as fw:
		for line in f:
			if line[0]=='>':
				cnt+=1
				fw.write(f">a_hp{hp}_{cnt}\n")
			elif ".fa" not in line:
				fw.write(line)




