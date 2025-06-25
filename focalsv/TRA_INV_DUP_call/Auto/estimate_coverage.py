import pysam
import numpy as np
import subprocess
from tqdm import tqdm
from joblib import Parallel, delayed
import os
def generate_genome_bins_from_bam(bam_path, bin_size):
    """
    Cut the genome (as defined in the BAM header) into fixed-size bins.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file.
    bin_size : int
        Desired size (bp) of each bin.

    Returns
    -------
    List[Tuple[str, int, int]]
        A list of (chrom, start, end) tuples covering the entire genome.
    """
    bins = []
    cnt = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, length in zip(bam.references, bam.lengths):
            cnt +=1
            if cnt>22:
                break
            for start in range(0, length, bin_size):
                end = min(start + bin_size, length)
                bins.append((chrom, start, end))
    return bins



def sample_unique_indices(n, k, seed=None):
    """
    Sample `n` unique integers from 1 to `k` inclusive using NumPy.
    
    Parameters
    ----------
    n : int
        Number of unique indices to sample.
    k : int
        Maximum index value (inclusive).
    seed : int or None
        Seed for reproducibility.
    
    Returns
    -------
    numpy.ndarray
        An array of shape (n,) containing unique values from 1 to k.
    """
    rng = np.random.default_rng(seed)
    # Generate n unique samples without replacement
    return rng.choice(np.arange(0, k ), size=n, replace=False)

# Example usage:
# indices = sample_unique_indices(n=5, k=100, seed=42)
# print(indices)



def write_bins_to_bed(bins, n,bed_path):
    """
    Write a list of bins to a BED‐format file.

    Parameters
    ----------
    bins : List[Tuple[str, int, int]]
        Output from generate_genome_bins_from_bam().
    bed_path : str
        Path to the BED file to create.
    """
    indices = sample_unique_indices(n, len(bins), seed=0)
    sampled_bins = []
    with open(bed_path, "w") as out:
        for i in indices:
            chrom, start, end = bins[i]
            sampled_bins.append(bins[i])
            out.write(f"{chrom}\t{start}\t{end}\n")
    return sampled_bins 


def compute_avg_depth(bam_path, chrom, start, end):
    """Run `samtools depth` on one interval and return (chrom, start, end, avg_depth)."""
    region = f"{chrom}:{start}-{end}"

    cmd = "samtools depth '{0}' -r {1} | awk '{{sum+=$3}} END {{ print sum/NR}}'".format(bam_path, region)
    try:
        output = subprocess.check_output(cmd, shell=True, text=True)
        avg_depth = float(output.strip())
        # print("Average depth:", avg_depth)
        return  avg_depth
    except subprocess.CalledProcessError as e:
        # print("Command failed:", e)
        return  -1
    

def estimate_bam_cov(bamfile, out_dir,t, ):

    os.system("mkdir -p "+out_dir )
    # 2) compute in parallel
    bins = generate_genome_bins_from_bam(bamfile, bin_size= 500000)
    chosen_bins = write_bins_to_bed(bins, 100, bed_path= out_dir+"/sampled_bin.bed" )
    results = Parallel(n_jobs=t)(
        delayed(compute_avg_depth)(bamfile, chrom, start, end)
        for chrom, start, end in tqdm(chosen_bins)
    )
    bed_file = out_dir+"/sampled_bin_cov.bed"
    covs = []
    with open(bed_file,'w') as f:
        for i in range(len(chosen_bins)):
            chrom, start, end = chosen_bins[i]
            cov = results[i]
            if cov!= -1:
                covs.append(cov)
            f.write(f"{chrom}\t{start}\t{end}\t{cov}\n")

    mean_cov = np.mean(covs)
    with open(out_dir+"/mean_cov",'w') as f:
        f.write(str(mean_cov))
    print(f"mean cov: {mean_cov} ({bamfile})")
    return mean_cov

if __name__ == '__main__':
    bam_list = [

        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395BL_Pacbio/minimap2/HCC1395BL_Pacbio.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_Pacbio/minimap2/HCC1395_Pacbio.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395BL_ONT/minimap2/HCC1395BL_ONT.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395_ONT/minimap2/HCC1395_ONT.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395T-N_HiFi/HCC1395-BL/minimap2_HCC1395-BL_Pacbio_CCS_hg38.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/HCC1395/HCC1395T-N_HiFi/HCC1395/minimap2_HCC1395_Pacbio_CCS_hg38.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/COLO829/HiFi/COLO829-BL/minimap2_COLO829-BL_Pacbio_CCS_hg38.bam",
        # "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/COLO829/HiFi/COLO829/minimap2_COLO829_Pacbio_CCS_hg38.bam",
        "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/COLO829/ONT/Normal/minimap2_COLO829_Normal_ONT_hg38.bam",
        "/lio/lfs/maiziezhou_lab/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/Datasets/Cancer_Data/COLO829/ONT/Tumor/minimap2_COLO829_Tumor_ONT_hg38.bam"
    ]
    for bamfile in bam_list:
        out_dir = os.path.dirname(bamfile)
        print(out_dir)
        # continue
        t = 50
        mean_cov = estimate_bam_cov(bamfile, out_dir,t)
        print(bamfile)
        print("mean cov:", mean_cov)
