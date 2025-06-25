



import argparse 
parser = argparse.ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument('--out_dir','-o')
parser.add_argument('--patient','-p', default='sample')
parser.add_argument('--state','-s', choices = ['Normal','Tumor'])
parser.add_argument('--dtype','-d', choices= ['hifi', 'ont','clr'])

args = parser.parse_args()
dtype = args.dtype

out_dir = args.out_dir
dtype = args.dtype
patient = args.patient
state = args.state




import os
import pandas as pd
import re

code_dir = os.path.dirname(os.path.realpath(__file__))+'/'




def chrom_sort_key(chrom):
    match = re.match(r'chr(\d+)$', chrom)
    return int(match.group(1)) if match else float('inf')

def sort_sv_df(df):
    df = df.copy()
    df['chrom_sort'] = df['chrom'].apply(chrom_sort_key)
    df = df.sort_values(by=['chrom_sort', 'start', 'end'])
    return df.drop(columns=['chrom_sort'])


def load_bed(bedfile, outfile, min_sup, min_mapq, min_size):
    # bedfile = f"/data/maiziezhou_lab/CanLuo/FocalSV/GR_Revision/ComplexSV/Reproduce_DUP/{patient}/out_{dtype.lower()}_{state}/{vtype}s.bed"
    

    df = pd.read_csv(bedfile, sep = '\t', header = None)
    df['size'] = df.iloc[:,2] - df.iloc[:,1]


    filter = (df.iloc[:,3]>=min_sup) & (df.iloc[:,4]>=min_mapq)  & (df.iloc[:,7]<= 160000000  ) & (df.iloc[:,7]>=min_size) 


    df = df[filter].reset_index(drop=True)
    df.columns = ['chrom','start','end','n_sup','mapq','std_left','std_right','size']
    df = sort_sv_df(df)
    df.to_csv(outfile, sep = '\t', index = False)
    
    cnt = df.shape[0]

    print(f"{patient} {dtype} {state} {vtype} FocalSV, num {vtype} = " ,cnt)

    return df

def reorder_bnd_columns(df):
    first_cols = ['chrom1', 'chrom2', 'pos1', 'pos2']
    other_cols = [col for col in df.columns if col not in first_cols]
    return df[first_cols + other_cols]




def chrom_sort_key(chrom):
    match = re.match(r'chr(\d+)$', chrom)
    return int(match.group(1)) if match else float('inf')

def sort_bnd_df(df):
    # Add numeric sort keys
    df = df.copy()
    df['chrom1_sort'] = df['chrom1'].apply(chrom_sort_key)
    df['chrom2_sort'] = df['chrom2'].apply(chrom_sort_key)

    # Sort by chrom1, chrom2, pos1, pos2
    df = df.sort_values(by=['chrom1_sort', 'chrom2_sort', 'pos1', 'pos2'])

    # Drop helper columns
    df = df.drop(columns=['chrom1_sort', 'chrom2_sort'])

    return df


def load_tsv(infile, outfile, min_sup, min_mapq, ):

    df = pd.read_csv(infile, sep = '\t', header = None)
    df.columns = ['chrom1','pos1','chrom2','pos2','bnd','n_sup','mapq','std_left','std_right']


    filter = (df.iloc[:,5]>=min_sup) & (df.iloc[:,6]>=min_mapq)  

    df = df[filter].reset_index(drop=True)

    cnt = df.shape[0]
    df = reorder_bnd_columns(df)
    df = sort_bnd_df(df)
    df.to_csv(outfile, sep = '\t', index = False)

    print(f"{patient} {dtype} {state} {vtype} FocalSV, num {vtype} = " ,cnt)
    return df




#-- hifi
if dtype == 'hifi':

    r_inv, min_mapq_inv, min_size_inv = 0.25, 60, 200
    r_bnd, min_mapq_bnd,  = 0.1, 55, 

#-- clr
elif dtype == 'clr':

    r_inv, min_mapq_inv, min_size_inv = 0.3, 58, 1000
    r_bnd, min_mapq_bnd,  = 0.15, 57, 

#-- ont
elif dtype == 'ont':

    r_inv, min_mapq_inv, min_size_inv = 0.35, 58, 1000
    r_bnd, min_mapq_bnd, = 0.2, 55, 


with open(out_dir+"/estimate_bam_cov/mean_cov",'r') as f:
    cov = eval(f.read())



for vtype in ['BND','INV']:
        
    print(f"\n-------------Processing {patient} {state} {vtype}")

    outfile = os.path.join(out_dir,f"{patient}_{state}_{dtype}_final_{vtype}.tsv")
    if  vtype == 'INV':
        min_sup_inv = cov * r_inv
        # print("min_sup_inv:",min_sup_inv)
        infile = out_dir+"/INVs.bed"
        dc_comp = load_bed(infile, outfile, min_sup_inv, min_mapq_inv, min_size_inv)
        # dc_comp = cluster_dict(dc_comp, std_t= 100)
    else:
        min_sup_bnd = cov * r_bnd
        # print("min_sup_bnd:",min_sup_bnd)
        infile = out_dir+"/TRA.tsv"
        dc_comp = load_tsv(infile, outfile, min_sup_bnd, min_mapq_bnd,)
    
    print(f"Save final {vtype} call to {outfile}")




