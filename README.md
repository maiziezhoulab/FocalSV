# FocalSV


FocalSV is a tool for region-based structural variant (SV) assembly and refinement using long-read sequencing data (PacBio HiFi, CLR, and ONT). It offers two modes - auto and target - to address diverse analytical goals. In target mode, users can specify regions of interest for SV detection, including deletions (DEL), insertions (INS), translocations (TRA), duplications (DUP), and inversions (INV). In auto mode, FocalSV autonomously detects and refines SV-rich regions by integrating population-level SV patterns with read-level signals from individual long-read data. 

## Table of Content 
- [Installation](#install-through-github)
- [Large INDEL detection](#Large-INDEL-detection)
    - [Step 0: Automatically detect target regions (for auto mode)](#Step-0-Automatically-detect-target-regions-for-auto-mode)
    - [Step 1: Generating Candidate SVs (for both modes)](#Step-1-Generating-Candidate-SVs-for-both-modes)
    - [Step 2: filtering and genotype correction (for both modes)](#Step-2-filtering-and-genotype-correction-for-both-modes)
- [TRA INV DUP detection](#TRA-INV-DUP-detection)
  -  [FocalSV-auto mode](#FocalSV-auto-mode)
  -  [FocalSV-target mode](#FocalSV-target-mode)
- [Troubleshooting](#Troubleshooting)


# Install through GitHub:

```
git clone  https://github.com/maiziezhoulab/FocalSV.git
```

## Dependencies for Github installation:
FocalSV integrated multiple state-of-the-art tools into the pipeline, including [longshot](https://github.com/pjedge/longshot), [hifiasm](https://github.com/chhylp123/hifiasm),[Flye](https://github.com/fenderglass/Flye),[Shasta](https://github.com/paoloshasta/shasta). To set up the environment, you need to have conda installed. Then, simply run
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
- You maybe download our test data from zenodo https://zenodo.org/records/15750913
  
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

By default, FocalSV(auto) uses a prior VCF generated from haplotype-resolved assemblies of nine individuals (HG01109, HG01243, HG02055, HG02080, HG02109, HG02145, HG02723, HG03098, and HG03492). SVs were identified from these assemblies and merged using the [Pangenie](https://github.com/eblerjana/pangenie) vcf-merging pipeline. This multi-sample VCF serves as the default population reference for FocalSV(auto), though users may substitute any population-scale, pangenome graph-based SV catalog suitable for their specific studies. To use the default prior VCF, specify either `FocalSV/Population_SV/pangenie_hg19_SV_gt50.vcf.gz` or `FocalSV/Population_SV/pangenie_hg38_SV_gt50.vcf.gz`. Alternatively, you may provide your own VCF file (both .vcf and .vcf.gz formats are supported).

The hg19 and hg38 references can be found on our zenodo link https://zenodo.org/records/15750913. You can also use your own matching BAM file and reference.

```
python3 FocalSV/focalsv/0_define_region.py \
--bam_file <wgs_bam> \
--ref_file <reference> \
--prior_file <popuplation SV file> \
--out_dir ./FocalSV_results/Define_Region \
--data_type HIFI \
--num_threads 8
--lib <lib_name>
```
The output file is `SV_Regions_<data_type>_<lib>.bed`.

## Step 1: Generating Candidate SVs (for both modes)

### Parameters

#### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region. Use 0 to select any chromosomes in the bed file.
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT).

#### Options for Region Selection:

##### For Single region analysis (Start/End Positions)

- **--region_start/-S**: Start coordinate of the target region (integer) (Required for single region mode).
- **--region_end/-E**: End coordinate of the target region (integer) (Required for single region mode).

##### For Multi-region analysis (BED file)

- **--target_bed/-target_bed**: BED file containing multiple target regions (optional for multi-region analysis in either target or auto mode. For auto mode, the BED file can be generated in Step 0)

#### Optional Parameters:

- **--out_dir/-o**: Output directory to store results (default: `./RegionBased_results`).
- **--num_cpus/-cpu**: Number of CPUs to use (default: 10).
- **--num_threads/-thread**: Number of threads (default: 8).

### Examples

#### 1. Running for One Region (Start-End)

```
python3 FocalSV/focalsv/focalsv.py \
--bam_file zenodo/HG002_HIFI_L1_chr21_hg19.bam  \
--ref_file zenodo/hg19_ref.fa \
--chr_num 21 \
--region_start 100000 \
--region_end 200000 \
--out_dir ./FocalSV_results/chr21_100000_200000 \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

#### 2. Running for Multiple Regions with a BED File
\*Note that you can provide a custom BED file based on your regions of interest, or directly use the whole-genome BED file generated in Step 0.

Here is an example of running multiple target regions on chromosome 21.
```
python3 FocalSV/focalsv/focalsv.py \
--bam_file zenodo/HG002_HIFI_L1_chr21_hg19.bam  \
--ref_file zenodo/hg19_ref.fa \
--chr_num 21 \
--target_bed zenodo/HG002_SV_rich_regions_chr21.bed \
--out_dir ./FocalSV_results/ \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

### Output:

```
FocalSV_results/chr21
  ├── results/
  │   └── FocalSV_Final_SV.vcf
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

- **`FocalSV_results/chr21/results/FocalSV_Final_SV.vcf`**  
  Final structural variant (SV) results for chr21.
  

#### `regions/`

- **`Region_chr21_S100000_E200000/`**  
  Results for chromosome 21, positions 100000-200000:
  - **`HP1.fa`**: Assembled contigs for haplotype 1.
  - **`PSxxx_hp1.fa`**: Contigs for phase block `xxx` and haplotype 1.
  - **`region.bam`**: Cropped BAM file for the region.
  - **`region_phased.bam`**: Phased BAM file for haplotype-specific alignments.

**Note**:

- For **single region analysis** (specified by a start and end coordinate), a corresponding region folder will be created under the `regions/` directory. For example, if you specify `--region_start 0 --region_end 200000 --chr_num 21`, the region will be named `Region_chr21_S0_E200000`.
- For **multi-region analysis using a BED file**, the number of regions corresponds to the number of lines in the BED file. Each region is named according to the start (S) and end (E) positions specified in its respective line.

#### `logs/`

- Contains log files for debugging and tracking pipeline steps.



#### 3. Running FocalSV on a whole-genome scale (except sex chromosomes)
\*Note that you can provide a custom BED file based on your regions of interest, or directly use the whole-genome BED file generated in Step 0.

```
python3 FocalSV/focalsv/focalsv.py \
--bam_file <wgs_bam> \
--ref_file <reference> \
--chr_num 0 \
--target_bed target_region_wgs.bed \
--out_dir ./FocalSV_results/ \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

#### Result

- **`FocalSV_results/FocalSV_Final_SV.vcf`**  
  Final structural variant (SV) results for the whole chromosome.

# TRA INV DUP detection

## FocalSV - auto mode

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

Here is an example of how to run FocalSV-auto to get TRA, INV, and DUP, on HCC1395.
```
python3 focalsv/TRA_INV_DUP_call/Auto/FocalSV-auto_TRA_INV_DUP_call.py \
--bam_file HCC1395_Pacbio_hg38.bam \
--data_type CLR \
--out_dir HCC1395_FocalSV-auto_tra_inv_dup_output 
--patient HCC1395 \
--state Tumor 
```
The output is `HCC1395_Tumor_CLR_final_TRA.tsv`, `HCC1395_Tumor_CLR_final_INV.tsv` and `HCC1395_Tumor_CLR_final_DUP.tsv`.


## FocalSV - target mode

For target mode,  you need to provide a BED file with special format. 
For DUP, each row is `chrom    start    end    DUP`.
For DUP, each row is `chrom    start    end    INV`.
For TRA, each row is `chrom1    start1    end1    chrom2    start2    end2    TRA`, where chrom1, start1, and end1 define the rough location for one breakend, chrom2, start2, and end2 define the rough location for the other breakend.

The target BED file for HCC1395 is provided in the zenode repo.


For target mode, you need to first perform FocalSV-target large indel call in the target regions of duplications, as we have a module of duplication recovery from insertions. You need to follow the step1 FocalSV-target mode to generate the large indel call result.

Here is an example of running the large indel call on HCC1395.
```
python3 FocalSV/focalsv/focalsv.py \
--bam_file HCC1395_Pacbio_hg38.bam \
--ref_file zenodo/hg38_ref.fa \
--chr_num 0 \
--bed_file zenodo/HCC1395_SV_rich_regions_DUP.bed \
--out_dir ./FocalSV_results_DUP \
--data_type CLR \
--num_cpus 10 \
--num_threads 8
```

If you want to run it on another sample, you need to provide the interested DUP target region BED file.
```
python3 FocalSV/focalsv/focalsv.py \
--bam_file<sample>_<dtype>_hg38.bam \
--ref_file zenodo/hg38_ref.fa \
--chr_num 0 \
--bed_file DUP_regions.bed \
--out_dir ./FocalSV_results_DUP \
--data_type CLR \
--num_cpus 10 \
--num_threads 8
```
`FocalSV_results_DUP/` will be the input for the next step.

Next, you can run the SV calling.
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
--input_dir FocalSV_results_DUP \
--bam_file HCC1395_Pacbio_hg38.bam \
--bed_file zenodo/HCC1395_SV_rich_regions.bed \
--data_type CLR \
--ref_file zenodo/hg38_ref.fa \
--out_dir HCC1395_FocalSV-target_tra_inv_dup_output
```
The output is `HCC1395_FocalSV-target_tra_inv_dup_output/FocalSV_TRA_INV_DUP.vcf`.

If you want to run it on another sample,
```
python3 focalsv/TRA_INV_DUP_call/Target/FocalSV-target_TRA_INV_DUP_call.py \
--input_dir FocalSV_results_DUP \
--bam_file <sample>_<datatype>_hg38.bam \
--bed_file <target_bed> \
--data_type <datatype> \
--ref_file zenodo/hg38_ref.fa \
--out_dir <sample>_FocalSV-target_tra_inv_dup_output
```



## Troubleshooting:

If you encounter issues, please submit them on the [FocalSV GitHub Issues](https://github.com/maiziezhoulab/FocalSV/issues) page.
