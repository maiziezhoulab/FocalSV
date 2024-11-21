# FocalSV

## Overview:

FocalSV is a tool for region-based structural variant (SV) assembly and refinement using long-read sequencing data. It provides a targeted approach to SV detection, focusing on specific genomic regions using a phased diploid assembly method. This software is optimized for multiple data types, including HiFi, CLR, and ONT sequencing.

## Dependencies:

The following tools and libraries are required to run FocalSV:

- Python 3.x
  - numpy
  - pysam
- SAMtools
- Longshot
- Flye
- Hifiasm
- Shasta
- Minimap2

To install the necessary dependencies, you can either ensure the executables are in your PATH environment variable or use the provided installation script.

## Installation:

You can install FocalSV by cloning the repository from GitHub and running the installation script to check for dependencies.

```
git clone https://github.com/maiziezhoulab/FocalSV.git
cd FocalSV
chmod +x install.sh
./install.sh
```

## Running the Code:

To execute the code, either add `FocalSV/bin` to your `.bashrc` file or use the full path to `main.py`.

## Step 1: Generating Candidate SVs

### Parameters

#### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region or whole chromosome analysis.

#### Options for Region Selection:

##### Single Region (Start/End Positions)

- **--region_start/-S**: Start index of the target region (integer) (Required for single region mode).
- **--region_end/-E**: End index of the target region (integer) (Required for single region mode).

##### Mutli Region (BED file)

- **--target_bed/-target_bed**: BED file with multiple target regions (Optional for multiple region mode).

#### Optional Parameters:

- **--out_dir/-o**: Output directory to store results (default: `./RegionBased_results`).
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT) (default: `HIFI`).
- **--num_cpus/-cpu**: Number of CPUs to use (default: 10).
- **--num_threads/-thread**: Number of threads (default: 8).

### Examples

#### 1. Running for One Region (Start-End)

```
python3 FocalSV/focalsv/main.py \
--bam_file ./test/test_hifi_chr21.bam \
--ref_file ./test/test_chr21.fa \
--chr_num 21 \
--region_start 100000 \
--region_end 200000 \
--out_dir ./FocalSV_results \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

#### 2. Running with a BED File (Multiple Regions)

```
python3 FocalSV/focalsv/main.py \
--bam_file ./test/test_hifi_chr21.bam \
--ref_file ./test/test_chr21.fa \
--chr_num 21 \
--target_bed ./test/test_chr21_multiRegion.bed \
--out_dir ./FocalSV_results \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

### Output:

```
FocalSV_results/
  ├── results/
  │   ├── FocalSV_Candidate_SV.vcf
  │   ├── FocalSV_Candidate_SV_redundancy.vcf
  │   └── variants.vcf
  ├── regions/
  │   ├── Region_chr21_S100000_E200000/
  │   │   ├── results/
  │   │   ├── HP1.fa
  │   │   ├── HP2.fa
  │   │   ├── PSxxx_hp1.fa
  │   │   ├── PSxxx_hp2.fa
  │   │   ├── region.bam
  │   │   ├── region_phased.bam
  │   │   └── ...
  │   ├── Region_chr21_Sxxx_Exxx/
  │   └── ...
  └── logs/
```

#### `results/`

- **`results/FocalSV_Candidate_SV.vcf`**  
  Candidate structural variant (SV) results without redundancy (key output).

#### `regions/`

- **`Region_chr21_S100000_E200000/`**  
  Results for chromosome 21, positions 100000-200000:
  - **`HP1_xxx.fa`**: Assembled contigs for haplotype 1.
  - **`PSxxx_hp1.fa`**: Contigs for phase block `xxx` and haplotype 1.
  - **`region.bam`**: Cropped BAM file for the region.
  - **`region_phased.bam`**: Phased BAM file for haplotype-specific alignments.

**Note**:

- For **one region analysis** (given a start and end region), there will be **one region** under the `regions/` folder. For example, if specified with `--region_start 0 --region_end 200000 --chr_num 21`, the region will be named `Region_chr21_S0_E200000`.
- For **analysis using a BED file**, the number of regions equals the number of lines in the BED file, with each region named based on the start (S) and end (E) positions specified in each line.

#### `logs/`

- Contains log files for debugging and tracking pipeline steps.

## Step 2: Filtering and genotype correction


FocalSV incorporates a post-processing module to filter false positives and correct genotypes further. This step involves collecting reads-based signatures from the read-to-reference BAM file. You can either run it by chromosome or on a whole genome scale. Note that the minimum scale is per chromosome, not per region, because read depth is not so accurate on the edge of each region and we try to minimize the effect of read depth fluctuation.

### Parameters

#### Required Parameters:
- **--bam_file/-bam**: The input BAM file.
- **--vcf_file/-vcf**: Path to the input FocalSV Candidate SV VCF file.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region or whole chromosome analysis.

#### Optional Parameters:
- **--out_dir/-o**: Output directory to store results (default: `./FocalSV_Final_VCF`).
- **--num_threads/-thread**: Number of threads (default: 8).


### Examples

#### 1. Running for One Chromosome

To achieve the most computation efficiency, if you have multiple target regions in one chromosome, you should merge the VCF file in those regions and only run the post-processing once. Below is an example command.
\*Note that the minimum scale is per chromosome, not per region, because read depth is not so accurate on the edge of each region and we try to minimize the effect of read depth fluctuation. If you only run one region, make sure to put the chromosome number for this target region for post-processing.

```
python3 ./FocalSV/focalsv/post_processing/FocalSV_Filter_GT_Correct.py \
--bam_file ./test/test_hifi_chr21.bam \
--ref_file ./test/test_chr21.fa \
--vcf_file FocalSV_results/results/FocalSV_Candidate_SV.vcf  \
--chr_num 21 \
--out_dir ./FocalSV_Final_VCF \
--data_type HIFI \
--num_threads 8
```

This command runs for chromosome 21 using 8 threads.

#### 2. Running for Whole Genome

Below is an example command for running post-processing one whole genome scale.

```
python3 ./FocalSV/focalsv/post_processing/FocalSV_Filter_GT_Correct.py \
--bam_file ${wgs_bam} \
--ref_file ${wgs_reference} \
--vcf_file FocalSV_results/results/FocalSV_Candidate_SV.vcf  \
--chr_num wgs \
--out_dir ./FocalSV_Final_VCF \
--data_type HIFI \
--num_threads 8
```

This command runs for whole genome using 10 threads.

### Output:

```
FocalSV_Final_SV/FocalSV_Final_SV.vcf
```

#### result

- **`FocalSV_Final_SV/FocalSV_Final_SV.vcf`**  
  Final structural variant (SV) results.
  
## Step Three (Optional): Evaluation

The `evaluation.py` script runs Truvari to assess the accuracy of detected structural variants (SVs) by comparing them to a benchmark vcf file. It calculates precision and recall metrics and performs data cleanup to retain only high-confidence results, streamlining the output for further analysis.

```
python3 FocalSV/focalsv/evaluation.py \
--chr_num 21 \
--bench_vcf ./test/test_targetSV_chr21.vcf \
--reference ./test/test_chr21.fa \
--out_dir ./FocalSV_results \
--num_cpus 10 \
--num_threads 8 \
--data_type HIFI
```

### Output

```
FocalSV_results/
  ├── eval/
  │   ├── Truvari_xxx/
  │   │   ├── DEL_50_/
  │   │   ├── INS_50_/
  │   │   ├── Truvari_results.xls
  │   ├── FocalSV_Candidate_SV_DEL_xxx.vcf
  │   ├── FocalSV_Candidate_SV_INS_xxx.vcf
  │   ├── FocalSV_Candidate_SV_DEL_INS_xxx.vcf
  │   └── ...
  ├── results/
  ├── regions/
  └── logs/
```

#### `eval/`

- **`eval/Truvari_xxx/Truvari_results.xls`**  
  Final Truvari evaluation metrics for insertions (INS) and deletions (DEL), including precision, recall, and F1 score.

## Troubleshooting:

If you encounter issues, please submit them on the [FocalSV GitHub Issues](https://github.com/maiziezhoulab/FocalSV/issues) page.
