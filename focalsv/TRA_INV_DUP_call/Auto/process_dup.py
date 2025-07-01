import os
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse
import subprocess
from estimate_coverage import estimate_bam_cov
# Set up base directories
global code_dir
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'



parser = argparse.ArgumentParser()
parser.add_argument('--bed_file','-bed')
parser.add_argument('--bam_file','-bam')
parser.add_argument('--patient', '-p',required=True, choices=['HCC1395', 'COLO829'])
parser.add_argument('--out_dir','-o')
parser.add_argument('--dtype', '-d',required=True, choices=['HIFI', 'ONT','CLR'])
parser.add_argument('--state','-s', required=True, choices=['Tumor', 'Normal'])
parser.add_argument('--n_thread','-t', type = int, default = 22 )
args = parser.parse_args()
bed_file = args.bed_file
bam_file = args.bam_file
patient = args.patient
dtype = args.dtype
state = args.state
out_dir = args.out_dir
t = args.n_thread
os.makedirs(out_dir, exist_ok=True)


#---------------------- STEP 1: Estimate AVG BAM Coverage ----------------------

def estimate_avg_cov(bam_file, out_dir,t):
    print("---------------------- STEP 1: Estimate AVG BAM Coverage ----------------------")
    estimate_bam_cov(bam_file, out_dir,t,)



# ---------------------- STEP 2: Filter Raw BED ----------------------



def first_round_filter(bedfile,out_dir,out_bed, patient, dtype, state):
    print("---------------------- STEP 2: Filter Raw BED ----------------------")
    # state_key = 'tm' if state == 'Tumor' else 'bl'
    # bedfile = f"/data/maiziezhou_lab/CanLuo/FocalSV/GR_Revision/ComplexSV/{patient}_hg38/out_{dtype}_{state_key}/DUPs.bed"
    # cov_map = get_cov_map()
    # cov = cov_map[(patient, dtype, state.lower())]
    with open(out_dir+"/estimate_bam_cov/mean_cov",'r') as f:
        cov = eval(f.read())

    if dtype == 'HIFI':
        min_sup, min_mapq, min_size = cov * 0.1, 50, 200
    elif dtype == 'CLR':
        min_sup, min_mapq, min_size = cov * 0.14, 50, 500
    elif dtype == 'ONT':
        min_sup, min_mapq, min_size = cov * 0.1, 40, 500

    df = pd.read_csv(bedfile, sep='\t', header=None)
    df.insert(3,'size',df[2] - df[1])
    df = df[(df[3] >= min_sup) & (df[4] >= min_mapq) & (df['size'] >= min_size)].copy()
    # out_bed = os.path.join(out_dir, f"{patient}_{state}_{dtype}_filtered_dup.bed")
    df.to_csv(out_bed, sep='\t', header=False, index=False)
    return out_bed

# ---------------------- STEP 3: Add Coverage ----------------------

def add_coverage(bed_file, out_file,bam_file, threads=50):
    print("---------------------- STEP 3: Add Coverage ----------------------")
    # out_file = bed_file.replace('.bed', '_cov.bed')
    script = os.path.join(code_dir, "bed_avg_depth.py")
    cmd = f"python3 {script} -bed {bed_file} -bam {bam_file} -o {out_file} -t {threads}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return out_file




# ---------------------- STEP 4: Filter Based on Coverage Features ----------------------

def get_cov_map():
    df = pd.read_csv(os.path.join(code_dir, "data_cov.csv"))
    return {(row.Patient, row.Dtype, row.Sample): row.Cov for _, row in df.iterrows()}

def second_round_filter(df_path, out_dir, out_file,patient, state, dtype):
    print("---------------------- STEP 4: Filter Based on Coverage Features ----------------------")
    # cov_map = get_cov_map()
    # cov = cov_map[(patient, dtype, state.lower())]
    with open(out_dir+"/estimate_bam_cov/mean_cov",'r') as f:
        cov = eval(f.read())
    print(f"Using estimated mean cov: {cov}")
    df = pd.read_csv(df_path, sep='\t', header=None)
    df.columns = ['chrom','start','end','size','n_sup','mapq','std_left','std_right','cov_sv','cov_left','cov_right','cov_diff']
    df['rel_n_sup'] = df['n_sup'] * 2 / (df['cov_left'] + df['cov_right'])
    df['rel_cov_diff'] = df['cov_sv'] * 2 / (df['cov_left'] + df['cov_right'])
    df['rel_cov_sv'] = df['cov_sv'] / cov
    df['rel_cov_left'] = df['cov_left'] / cov
    df['rel_cov_right'] = df['cov_right'] / cov
    df['rel_std'] = df[['std_left','std_right']].min(axis=1) / df['n_sup']
    
    # --- Apply tuned filters
    if (dtype == 'HIFI') and (state == 'Tumor'):
        f = (df['rel_cov_diff'].between(1.1, 3)) & (df['rel_n_sup'].between(0.25, 1.5)) & \
            (df['mapq'] > 59.8) & (df['rel_cov_sv'].between(0.6, 5)) & \
            ((df['std_left'] < 1.4) | (df['std_right'] < 1.4))
    elif (dtype == 'HIFI') and (state == 'Normal'):
        f = (df['rel_cov_diff'].between(1.3, 4)) & (df['rel_n_sup'].between(0.25, 1.2)) & \
            (df['mapq'] > 59.5) & (df['rel_cov_sv'].between(1, 4)) & \
            ((df['std_left'] < 1.4) | (df['std_right'] < 1.4))
    elif (dtype == 'CLR') and (state == 'Tumor'):
        f = (df['rel_cov_diff'].between(1.15, 8)) & (df['rel_n_sup'].between(0.22, 4.6)) & \
            (df['mapq'] > 50) & (df['rel_std'] < 2) & \
            (df['size'].between(3000, 35e6)) & (df['rel_cov_sv'].between(0.7, 9)) & \
            ((df['std_left'] < 25) | (df['std_right'] < 25))
    elif (dtype == 'CLR') and (state == 'Normal'):
        f = (df['rel_cov_diff'].between(1.15, 8)) & (df['rel_n_sup'].between(0.22, 4.6)) & \
            (df['mapq'] > 50) & (df['rel_std'] < 0.5) & \
            (df['size'].between(3000, 35e6)) & (df['rel_cov_sv'].between(0.8, 4)) & \
            ((df['std_left'] < 15) | (df['std_right'] < 15))
    elif (dtype == 'ONT') and (state == 'Tumor'):
        f = (df['rel_cov_diff'].between(1.15, 8)) & (df['rel_n_sup'].between(0.22, 4.6)) & \
            (df['mapq'] > 50) & (df['rel_std'] < 2) & \
            (df['size'].between(3000, 35e6)) & (df['rel_cov_sv'].between(0.7, 9)) & \
            ((df['std_left'] < 25) | (df['std_right'] < 25))
    elif (dtype == 'ONT') and (state == 'Normal'):
        f = (df['rel_cov_diff'].between(1.15, 8)) & (df['rel_n_sup'].between(0.22, 4.6)) & \
            (df['mapq'] > 50) & (df['rel_std'] < 0.5) & \
            (df['size'].between(3000, 35e6)) & (df['rel_cov_sv'].between(0.8, 4)) & \
            ((df['std_left'] < 15) | (df['std_right'] < 15))

    df_filtered = df[f].copy()
    
    df_filtered.to_csv(out_file, sep='\t', index=False)
    print(f"[Done] Final SVs saved to: {out_file}")


# ---------------------- Run Pipeline ----------------------
bed_filtered =  os.path.join(out_dir, f"{patient}_{state}_{dtype}_filtered_DUP.bed")
bed_cov = bed_filtered.replace('.bed','_cov.bed')
bed_final = os.path.join(out_dir, f"{patient}_{state}_{dtype}_final_DUP.tsv")

estimate_avg_cov(bam_file, out_dir+"/estimate_bam_cov",t)
first_round_filter(bed_file, out_dir,bed_filtered, patient, dtype, state)
add_coverage(bed_filtered,bed_cov, bam_file,t)
second_round_filter(bed_cov, out_dir, bed_final, patient, state, dtype)


