
ref=$1
# /data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa
seq=$2
# /data/maiziezhou_lab/Datasets/Assemblies/HG002-T2T.v1.1.fasta
outbam=HG002-T2T_hg19.bam

#NOTE: the minimap2 command was from HiCanu paper assemblies to ref alignment, it was also used in PAV pipeline
# Paper: https://genome.cshlp.org/content/30/9/1291.full.pdf
# Section: Commands for identifying contig ends

# eqx may make the downstream region parsing slow as it gives much more accurate information (the "X" CIGAR)
#minimap2 --secondary=no -a --eqx -Y -x asm20 -s 200000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 ${ref} ${seq} | samtools sort -o ${outbam}

minimap2 --secondary=no -a -Y -x asm20 -s 200000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 ${ref} ${seq} | samtools sort -o ${outbam}

samtools index ${outbam}
