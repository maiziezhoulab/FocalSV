from argparse import ArgumentParser
from joblib import Parallel, delayed
from tqdm import tqdm
import glob
import os
import pandas as pd
from subprocess import Popen
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
# parser.add_argument('--input_path','-i',help = "config file, every row is a vcf path")
parser.add_argument('--input_dir','-i',help = "somatic call dir")
parser.add_argument('--n_thread','-t', type = int, default = 10, help = "default = 10")
args = parser.parse_args()
# input_path = args.input_path
input_dir = args.input_dir
n_thread = args.n_thread

# with open(input_path,'r') as f:
# 	vcf_list = [line.strip() for line in f if  line[0]!='#']

import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


def run_one_vcf(vcf_path):
	cmd = f"python3 {code_dir}/eval_real_data.py -i "+vcf_path
	Popen(cmd, shell = True).wait()
	return 



def merge_eval_summaries_to_excel(eval_dir, output_excel):
	svtypes_order = ['INV', 'DUP', 'TRA']
	datatypes_order = ['Hifi', 'CLR', 'ONT']
	sheet_order = [(dt, sv)  for sv in svtypes_order for dt in datatypes_order]

	# Structure: {(datatype, svtype): {tool: column_values}}
	merged_dict = {(dt, sv): {} for dt, sv in sheet_order}

	for subdir in os.listdir(eval_dir):
		subpath = os.path.join(eval_dir, subdir)
		if not os.path.isdir(subpath):
			continue

		try:
			datatype, tool = subdir.split('_', 1)
		except ValueError:
			continue  # skip folders not matching pattern

		if datatype not in datatypes_order:
			continue

		summary_path = os.path.join(subpath, 'summary.tsv')
		if not os.path.isfile(summary_path):
			continue
		
		try:
			df = pd.read_csv(summary_path, sep='\t')
		except:
			print(summary_path)
			exit()
		if 'metric' not in df.columns:
			continue

		for sv in svtypes_order:
			if sv not in df.columns:
				continue
			temp = df[['metric', sv]].rename(columns={sv: tool})
			if tool not in merged_dict[(datatype, sv)]:
				merged_dict[(datatype, sv)][tool] = temp
			else:
				merged_dict[(datatype, sv)][tool] = pd.merge(
					merged_dict[(datatype, sv)][tool],
					temp, on='metric', how='outer'
				)

	tool_order = ['FocalSV-target','FocalSV-auto','SVIM-asm','sawfish','cuteSV','SVIM','pbsv','Sniffles2']
	# Write to Excel
	with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
		for dt, sv in sheet_order:
			tool_dfs = merged_dict[(dt, sv)]
			if not tool_dfs:
				continue
			merged_df = None
			for tool in tool_order:
				if tool in tool_dfs:
					df_tool = tool_dfs[tool]
					if merged_df is None:
						merged_df = df_tool
					else:
						merged_df = pd.merge(merged_df, df_tool, on='metric', how='outer')
			merged_df.iloc[-3:, 1:] = merged_df.iloc[-3:, 1:].applymap(lambda x: round(x,4) if isinstance(x, (int, float)) else x)
			sheet_name = f"{sv}_{dt}"
			merged_df.to_excel(writer, sheet_name=sheet_name, index=False)

# Example usage:

vcf_list = glob.glob(input_dir+"/*/variants_all_somatic.vcf")
# vcf_list = [ vcf for vcf in vcf_list if 'sawfish' in vcf]
sequences = Parallel(n_jobs=n_thread)(delayed(run_one_vcf)(vcf) for vcf in tqdm(vcf_list))
merge_eval_summaries_to_excel(input_dir, input_dir+"/summary.xlsx")
