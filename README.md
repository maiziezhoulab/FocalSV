# FocalSV Assembly

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

### Example: Running for One Region (Start-End)

```
python3 FocalSV/main.py \
--bam_file sample.bam \
--ref_file genome.fa \
--chr_num 21 \
--region_start 100000 \
--region_end 200000 \
--out_dir ./FocalSV_results \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

### Example: Running with a BED File (Multiple Regions)

```
python3 FocalSV/main.py \
--bam_file sample.bam \
--ref_file genome.fa \
--chr_num 21 \
--target_bed regions.bed \
--out_dir ./FocalSV_results \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

### Example: Running Whole Chromosome Evaluation

```
python3 FocalSV/main.py \
--bam_file sample.bam \
--ref_file genome.fa \
--chr_num 21 \
--eval \
--bed_file regions.bed \
--vcf_file variants.vcf \
--svlen_threshold 50 \
--out_dir ./FocalSV_results \
--data_type HIFI \
--num_cpus 10 \
--num_threads 8
```

### Required Parameters:

- **--bam_file/-bam**: The input BAM file.
- **--ref_file/-r**: Reference FASTA file.
- **--chr_num/-chr**: Chromosome number for the target region or whole chromosome analysis.

### Options for Region Selection:

- **--region_start/-S**: Start index of the target region (integer) (Required for single region mode).
- **--region_end/-E**: End index of the target region (integer) (Required for single region mode).
- **--target_bed/-target_bed**: BED file with multiple target regions (Optional for multiple region mode).

### Optional Parameters:

- **--out_dir/-o**: Output directory to store results (default: `./RegionBased_results`).
- **--data_type/-d**: Type of sequencing data (HIFI, CLR, ONT) (default: `HIFI`).
- **--num_cpus/-cpu**: Number of CPUs to use (default: 10).
- **--num_threads/-thread**: Number of threads (default: 8).
- **--eval/-e**: Flag to perform whole chromosome evaluation (optional).
- **--bed_file/-bed**: Bed file for whole chromosome evaluation.
- **--svlen_threshold/-sv_thresh**: SV length threshold for filtering variants in whole chromosome mode (default: 50).
- **--vcf_file/-v**: VCF file called by FreeBayes (required for whole chromosome mode).
- **--flanking/-flank**: Length of the flanking region around the target region (default: 50000).

## Output

The output of the pipeline is stored in the specified `out_dir` and includes intermediate files such as cropped BAM files, phased results, and potential final variant calls.

### Output Directory:

The default directory for output is `./RegionBased_results`, but this can be modified using the `--out_dir` option.

### Example Output Directory Structure:

```
FocalSV_results/
  ├── regions/
  │   ├── Region_chr21_S100000_E200000/
  │   │   ├── results/
  │   │   │   ├── final_vcf/
  │   │   │   │   ├── eval/
  │   │   │   │   ├── dippav_variant_no_redundancy.vcf
  │   │   │   │   └── dippav_variant_redundancy.vcf
  │   │   │   ├── DipPAV_FULL.sorted.bam
  │   │   │   └── DipPAV_FULL.sorted.bam.bai
  │   │   ├── HP1_xxx.fa
  │   │   ├── HP2_xxx.fa
  │   │   ├── region.bam
  │   │   └── region.bam.bai
  ├── target_sv.vcf
  └── logs/
```

## Troubleshooting:

If you encounter issues, please submit them on the [FocalSV GitHub Issues](https://github.com/maiziezhoulab/FocalSV/issues) page.

This README file has been updated to reflect the three different modes of running FocalSV: processing a single region, processing multiple regions using a BED file, or performing a whole chromosome evaluation.
