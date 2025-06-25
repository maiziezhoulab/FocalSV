from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
# parser.add_argument('--tools_list','-tl')
parser.add_argument('--input_dir','-i', help = "The truvari eval dir for one library across tools. Under this folder, each subfolder should be each tool name. Under each subfolder, there is INS_50_ and DEL_50_ folders containing the truvari result.")
parser.add_argument('--output_dir','-o')
args = parser.parse_args()
# tools_list = args.tools_list
input_dir = args.input_dir
output_dir = args.output_dir

import os
import numpy as np
import pandas as pd
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

# with open(tools_list,'r') as f:
# 	s = f.readlines()
tool_list = os.listdir(input_dir)

for tool in tool_list:
	tool_dir = f"{input_dir}/{tool}"
	cmd = '''python3 %s/line_bar_plot_one_tool.py -ins %s -del %s -o %s --prefix "%s"'''%(code_dir,tool_dir+'/INS_50_',
            tool_dir+'/DEL_50_',output_dir,tool)
	print(cmd)
	os.system(cmd)

cnt = 0 
for i in range(len(s)):
	prefix = s[i].split()[0]
	csv = output_dir+'/%s.csv'%prefix
	if os.path.exists(csv):
		df_new = pd.read_csv(csv)
		cnt+=1
	else:
		print( f"{csv} not exits")
		continue

	if cnt ==1:
		df = df_new[["lenth_range","SV_number"]]
	df.loc[:,'F1_%s'%prefix] = df_new['F1']
	df.loc[:,'Recall_%s'%prefix] = df_new['Recall']
	df.loc[:,'Precision_%s'%prefix] = df_new['Precision']
	df.loc[:,'TP_%s'%prefix] = df_new['TP']
	df.loc[:,'FP_%s'%prefix] = df_new['FP']
	df.loc[:,'FN_%s'%prefix] = df_new['FN']

start_list = []
end_list = []
for rg in df.lenth_range.values:
    start_list.append(str(eval(rg)[0]))
    end_list.append(str(eval(rg)[1]))
df['start']=start_list
df['end']=end_list 

df.to_csv(output_dir+"/tools.csv",index=False)

		

