# Data preparation
We provide the hg19 and hg38 reference on our [Zenodo repo](https://zenodo.org/records/15750913)

## HG002 Hifi_L1
```
# download fastq
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64011_190830_220126.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64011_190901_095311.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64012_190920_173625.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64012_190921_234837.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64015_190920_185703.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64015_190922_010918.fastq.gz

# alignment
ref=zenodo/hg19_ref.fa
minimap2 -t 30 --MD -Y -L -a -x map-hifi ${ref} ${reads} | samtools sort -o ${prefix}.bam
samtools index ${prefix}.bam
```
For HG002, we used hg19_ref.fa on our Zenodo link as the reference.

## HCC1395

```
# download fastq
prefetch $srid -O ./ --max-size 500G
fasterq-dump --split-3 -e 20 ./$srid
```
prefetch and fasterq-dump are included in the sra-tools package. You may use `conda install -c bioconda sra-tools` to install them.

Here are the srids for different libraries:
```
SRR8955953 pacbio tumor
SRR8955954 pacbio normal
SRR16005301 ONT tumor
SRR17096031 ONT normal
```
```
ref=zenodo/hg38_ref.fa
# pacbio alignment
minimap2 -t 30 --MD -Y -L -a -x map-pb ${ref} ${reads} | samtools sort -o ${prefix}.bam
samtools index ${prefix}.bam
# ONT alignment
minimap2 -t 30 --MD -Y -L -a -x map-ont ${ref} ${reads} | samtools sort -o ${prefix}.bam
samtools index ${prefix}.bam
```
For HCC1395, we used hg38_ref.fa on our zenodo link as the reference.
