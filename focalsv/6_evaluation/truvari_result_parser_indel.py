from argparse import ArgumentParser
import pandas as pd
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dir','-i')
# parser.add_argument('--output_dir','-o')

args = parser.parse_args()
input_dir = args.input_dir
# output_dir = args.output_dir

key_list = ["TP-call",
    "FP",
    "FN",
    "recall",
    "precision",
    "f1",
    "TP-call_TP-gt",
    "TP-call_FP-gt",
    "FN",
    "gt_recall",
    "gt_precision",
    "gt_f1"]

def load_summary(summary_path):
    with open(summary_path,'r')as f:
        s = f.read()
    dc = eval(s)

    metric_list = [ dc[key] for key in key_list]

    return metric_list

summary_path_ins = input_dir+'/INS_50_/summary.txt'
summary_path_del = input_dir+'/DEL_50_/summary.txt'
metric_ins = load_summary(summary_path_ins)
metric_del = load_summary(summary_path_del)


df = pd.DataFrame({"metric": key_list,
    "INS_50_": metric_ins,
    "DEL_50_": metric_del})

df.to_excel(input_dir+"/Truvari_results.xls")  

