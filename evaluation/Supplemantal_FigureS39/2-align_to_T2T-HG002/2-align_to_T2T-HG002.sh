# indexing with bwa, only need to run once
HG002-T2T=$1
bwa index ${HG002-T2T}
# /data/maiziezhou_lab/Datasets/Assemblies/HG002-T2T.v1.1.fasta

# Align strand seq to HG002 T2T
for cell in `cat ../0-strand-seq_reads/cells.list`; 
do 
bwa mem -t 20 ${HG002-T2T} ../0-strand-seq_reads/${cell}_1_sequence.fastq.gz ../0-strand-seq_reads/${cell}_2_sequence.fastq.gz\
| samtools sort -@4 -o ${cell}.bam; samtools index ${cell}.bam; 
done

# build cell bam list
for i in `cat ../0-strand-seq_reads/cells.list` ; do echo -e "${i}\t$(pwd)/${i}.bam" >> cell_bam.list; done

# collect informative cell:
python ../bin/collect_informative_cell_t2t.py -i cell_wc_ratio.pkl
