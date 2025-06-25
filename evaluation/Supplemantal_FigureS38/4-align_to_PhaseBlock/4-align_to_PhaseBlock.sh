# cat HP1.fa | grep "^>" | awk -F "_" '{print $2 "\t" $3 "\t" $4}' | sort -u > SV_regions.bed
# why I got an additional space between fields with awk -F "_" '{print $2, "\t", $3, "\t", $4}' ?
# 因为在 awk 里，用逗号分隔多个表达式时，输出时会在它们之间自动插入“输出字段分隔符”（OFS），默认就是一个空格
ASM=$1
ref_hg19=$2 #hg19
ref_t2t=$3
for i in {1..22}
do
	python ../bin/fix_folded_fasta.py -i ${ASM}/chr${i}_HP1.fa -o chr${i}_HP1.fa --hap_tag HP1 --fai ${ref_hg19}.fai
	python ../bin/fix_folded_fasta.py -i ${ASM}/chr${i}_HP2.fa -o chr${i}_HP2.fa --hap_tag HP2 --fai ${ref_hg19}.fai
	cat chr${i}_HP*.fa > chr${i}_PhaseBlocks.fa; cat chr${i}_HP1.fa | grep "^>" | awk -F "_" '{print $2 "\t" $3 "\t" $4}' | sort -u > chr${i}_SV_regions.bed
	rm chr${i}_HP*.fa
done

cat chr*_SV_regions.bed > SV_regions.bed


# sbatch submit_read_extraction.slurm
python ../bin/extract_local_strand_seq.py \
--t2t_fai ${ref_t2t}.fai \
--sv_region_bed SV_regions.bed \
--informative_cell ../2-align_to_T2T-HG002/informative_cell.pkl \
--t2t_bam ../3-align_T2T-HG002_to_ref/HG002-T2T_hg19.bam \
--strand_to_t2t_dir ../2-align_to_T2T-HG002/ \
--strand_seq_dir ../1-strand-seq_reads \
--out_dir ./ \
--rewrite_t2t_region


# sbatch submit_bwa_index.slurm

for i in {1..22}
do
	bwa index chr${i}_PhaseBlocks.fa
done


for i in {1..22}
do
	cd chr${i}
	ls *.fastq | awk -F "_" '{print $1}' | sort -u > cell.list
	for cell in `cat cell.list`
	do
		bwa mem -t 20 ../chr${i}_PhaseBlocks.fa ${cell}_1.fastq ${cell}_2.fastq | samtools sort -@4 -o ${cell}.bam; samtools index ${cell}.bam
	done
	cd -
done

python ../bin/calculate_phasing_quality.py\
--informative_cell ../2-align_to_T2T-HG002/informative_cell.pkl \
--input_dir ./ \
> contig_phase_qual.log

