# FocalSV


FocalSV is a tool for region-based structural variant (SV) assembly and refinement using long-read sequencing data (PacBio HiFi, CLR, and ONT). It offers two modes - auto and target - to address diverse analytical goals. In target mode, users can specify regions of interest for focused SV detection. In auto mode, FocalSV autonomously detects and refines SV-rich regions by integrating population-level SV patterns with read-level signals from individual long-read data. 

## Table of Content 
- [Installation](#install-through-github)
- [Large INDEL detection](#Large-INDEL-detection)
    - [Step 0: Automatically detect target regions (for auto mode)](#Step-0-Automatically-detect-target-regions-for-auto-mode)
    - [Step 1: Generating Candidate SVs (for both modes)](#Step-1-Generating-Candidate-SVs-for-both-modes)
    - [Step 2: filtering and genotype correction (for both modes)](#Step-2-filtering-and-genotype-correction-for-both-modes)
- [TRA INV DUP detection](#TRA-INV-DUP-detection)
  -  [FocalSV-target mode](#FocalSVtarget-mode)
  -  [FocalSV-auto mode](#FocalSVauto-mode)
- [Troubleshooting](#Troubleshooting)


# Install through GitHub:

```
git clone  https://github.com/maiziezhoulab/FocalSV.git
```

## Dependencies for Github installation:
FocalSV utilizes **Python3.8.3**. To set up the environment, you need to have conda installed. Then, simply run
```
conda env create -f FocalSV/requirement.yaml
```

Then you will have a virtual environment called **`FocalSV`** created. **Before running any FocalSV commands, please activate this environment first**.
```
conda activate FocalSV
```


## Running the Code:

To execute the code, either add `FocalSV/focalsv` to your `.bashrc` file or use the full path.


## Before you start

- FocalSV is intended for use on autosomal chromosomes only, as phasing on sex chromosomes involves additional complexities.
- The input for each step can be a whole-genome sequencing (WGS) BAM file. In our test example, we used a chr21 BAM file solely for convenience due to its smaller size. However, manually splitting the WGS BAM file by chromosome is not required.
# Large INDEL detection
## Step 0: Automatically detect target regions (for auto mode)

If you already have target regions of interest, you can skip this step and proceed directly to Step 1. Otherwise, run this step to automatically detect potential SV regions in the auto mode.
By default, this script takes a whole-genome BAM file and scans for potential SV regions across the entire genome.

### Parameters

#### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--ref_file/-r**: Reference FASTA file.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).
- **--prior_file/-p**: population SV VCF file. either .vcf or .vcf.gz is accepted.
- **--out_dir/-o**:  Output directory to store results.
- **--lib/-l**: library name.


#### Optional Parameters:

- **--num_threads/-thread**: Number of threads (default: 8).

### Examples

By default, FocalSV(auto) uses a prior VCF generated from haplotype-resolved assemblies of nine individuals (HG01109, HG01243, HG02055, HG02080, HG02109, HG02145, HG02723, HG03098, and HG03492). SVs were identified from these assemblies and merged using the Pangenie vcf-merging pipeline. This multi-sample VCF serves as the default population reference for FocalSV(auto), though users may substitute any population-scale, pangenome graph-based SV catalog suitable for their specific studies. To use the default prior VCF, specify either `FocalSV/Population_SV/pangenie_hg19_SV_gt50.vcf.gz` or `FocalSV/Population_SV/pangenie_hg38_SV_gt50.vcf.gz`. Alternatively, you may provide your own VCF file (both .vcf and .vcf.gz formats are supported).

```
python3 FocalSV/focalsv/0_define_region.py \
--bam_file <wgs_bam> \
--ref_file <reference> \
--prior_file <popuplation SV file> \
--out_dir ./FocalSV_results/Define_Region \
--data_type HIFI \
--num_threads 8
```
The output file is `SV_Regions_<data_type>_<lib>.bed`.

## Step 1: Generating Candidate SVs (for both modes)

### Parameters

#### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).

#### Options for Region Selection:

##### For Single region analysis (Start/End Positions)

- **--region_start/-S**: Start index of the target region (integer) (Required for single region mode).
- **--region_end/-E**: End index of the target region (integer) (Required for single region mode).

##### For Multi-region analysis (BED file)

- **--target_bed/-target_bed**: BED file containing multiple target regions (optional for multi-region analysis in either target or auto mode. For auto mode, the BED file can be generated in Step 0)

#### Optional Parameters:

- **--out_dir/-o**: Output directory to store results (default: `./RegionBased_results`).
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

#### 2. Running for Multiple Regions with a BED File

Here is an example of running all target region on chromosome 21.
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

If you want to run FocalSV on a whole-genome scale (except sex chromosomes), here is an example

```
for i in {1..22}
do
python3 FocalSV/focalsv/main.py \
--bam_file <wgs_bam> \
--ref_file <reference> \
--chr_num $i \
--target_bed target_region_chr${i}.bed \
--out_dir ./FocalSV_results/chr${i} \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
done
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
  Candidate structural variant (SV) results without redundancy (key output). This VCF will be the input of step2 for filtering and genotype correction.

#### `regions/`

- **`Region_chr21_S100000_E200000/`**  
  Results for chromosome 21, positions 100000-200000:
  - **`HP1_xxx.fa`**: Assembled contigs for haplotype 1.
  - **`PSxxx_hp1.fa`**: Contigs for phase block `xxx` and haplotype 1.
  - **`region.bam`**: Cropped BAM file for the region.
  - **`region_phased.bam`**: Phased BAM file for haplotype-specific alignments.

**Note**:

- For **single region analysis** (specified by a start and end coordinate), a corresponding region folder will be created under the `regions/` directory. For example, if you specify `--region_start 0 --region_end 200000 --chr_num 21`, the region will be named `Region_chr21_S0_E200000`.
- For **multi-region analysis using a BED file**, the number of regions corresponds to the number of lines in the BED file. Each region is named according to the start (S) and end (E) positions specified in its respective line.

#### `logs/`

- Contains log files for debugging and tracking pipeline steps.

## Step 2: filtering and genotype correction (for both modes)


FocalSV incorporates a post-processing module to filter false positives and correct genotypes further. This step involves collecting read-based signatures from the read-to-reference BAM file. You can either run it by chromosome or on a whole-genome scale. 

### Parameters

#### Required Parameters:
- **--bam_file/-bam**: The input BAM file.
- **--vcf_file/-vcf**: Path to the input FocalSV Candidate SV VCF file generated in Step 1.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region or whole chromosome analysis.

#### Optional Parameters:
- **--out_dir/-o**: Output directory to store results (default: `./FocalSV_Final_VCF`).
- **--num_threads/-thread**: Number of threads (default: 8).


### Examples

#### 1. Running for One Chromosome

To achieve maximum computational efficiency, if you have multiple target regions on a single chromosome, you should merge the VCF files in those regions and run the post-processing only once. Below is an example command.
\*Note that the minimum scale is defined per chromosome, not per region, to reduce the impact of read depth fluctuations, which are less reliable at region boundaries. If you are analyzing a single region, be sure to specify the corresponding chromosome number for post-processing.

```
python3 ./FocalSV/focalsv/5_post_processing/FocalSV_Filter_GT_Correct.py \
--bam_file ./test/test_hifi_chr21.bam \
--ref_file ./test/test_chr21.fa \
--vcf_file FocalSV_results/results/FocalSV_Candidate_SV.vcf  \
--chr_num 21 \
--out_dir ./FocalSV_results/Final_VCF \
--data_type HIFI \
--num_threads 8
```

This command runs for chromosome 21 using 8 threads.

#### 2. Running for Whole Genome

Below is an example command for running post-processing on a whole-genome scale.

```
python3 ./FocalSV/focalsv/post_processing/FocalSV_Filter_GT_Correct.py \
--bam_file ${wgs_bam} \
--ref_file ${wgs_reference} \
--vcf_file FocalSV_results/results/FocalSV_Candidate_SV.vcf  \
--chr_num wgs \
--out_dir ./FocalSV_results/Final_VCF \
--data_type HIFI \
--num_threads 8
```

This command runs for the whole genome using 8 threads.

### Output:

```
FocalSV_results/
  |—— Final_SV/FocalSV_Final_SV.vcf
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

#### Result

- **`FocalSV_results/Final_SV/FocalSV_Final_SV.vcf`**  
  Final structural variant (SV) results.
  

# TRA INV DUP detection

## FocalSV(target) mode



### Parameters

#### Required Parameters:
- **--input_dir/-i**: FocalSV-target large indel call output folder
- **--bam_file/-bam**: The input BAM file.
- **--excel_file/-excel**: an Excel file that includes the ground truth SV (only for HCC1395)
- **--bed_file/-bed**: a bed file with interested target region. After each region, there should be SV type. For example, a row in the bed file can be "chr1 1000 5000 TRA", "chr2 1000 5000 INV", or "chr3 1000 5000 DUP"
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).
- **--ref_file/-r**: Reference FASTA file.
- **--out_dir/-o**: Output directory to store results.

#### Optional Parameters:
- **--num_threads/-thread**: Number of threads (default: 8).


### Examples

Here is an example of how to run FocalSV-target to get TRA INV and DUP on HCC1395.
```
python3 focalsv/TRA_INV_DUP_call/Target/FocalSV-target_TRA_INV_DUP_call.py \
--input_dir HCC1395_FocalSV-target_largeindel_output \
--bam_file HCC1395_Pacbio_hg38.bam \
--excel_file focalsv/TRA_INV_DUP_call/Target/High_confidence_callset.xlsx \
--data_type CLR \
--ref_file <hg38_reference> \
--out_dir HCC1395_FocalSV-target_tra_inv_dup_output
```
The output is `HCC1395_FocalSV-target_tra_inv_dup_output/FocalSV_TRA_INV_DUP.vcf`.

If you want to run it on another sample, you should provide a BED file with each region annotated with the SV type. For example, a row in the bed file can be "chr1 1000 5000 TRA", "chr2 1000 5000 INV", or "chr3 1000 5000 DUP"
```
python3 focalsv/TRA_INV_DUP_call/Target/FocalSV-target_TRA_INV_DUP_call.py \
--input_dir <sample>_FocalSV-target_largeindel_output \
--bam_file <sample>_<datatype>_hg38.bam \
--bed_file <target_bed> \
--data_type <datatype> \
--ref_file <hg38_reference> \
--out_dir <sample>_FocalSV-target_tra_inv_dup_output
```

## FocalSV(auto) mode

### Parameters

#### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).
- **--out_dir/-o**: Output directory to store results.
- **--patient/-p**: patient name.
- **--state/-s**: Tumor / Normal

#### Optional Parameters:
- **--num_threads/-thread**: Number of threads (default: 8).


### Examples

Here is an example of how to run FocalSV-auto to get TRA INV and DUP on HCC1395.
```
python3 focalsv/TRA_INV_DUP_call/Auto/FocalSV-auto_TRA_INV_DUP_call.py \
--bam_file HCC1395_Pacbio_hg38.bam \
--data_type CLR \
--out_dir HCC1395_FocalSV-auto_tra_inv_dup_output 
--patient HCC1395 \
--state Tumor 
```
The output is `HCC1395_FocalSV-auto_tra_inv_dup_output/FocalSV_TRA_INV_DUP.vcf`.


## Troubleshooting:

If you encounter issues, please submit them on the [FocalSV GitHub Issues](https://github.com/maiziezhoulab/FocalSV/issues) page.
