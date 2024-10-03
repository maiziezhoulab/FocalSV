#####################################################################################################################
#
#                       truvari evaluation
#
#####################################################################################################################


chr_num=$1
input_dir=$2
work_dir=$3
prefix=$4
dist=$5
p=$6
P=$7
S=$8
O=$9
source activate /home/liuy120/.conda/envs/bioinfo_temp
aligner=NGMLR
ref=/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa
bed=/data/maiziezhou_lab/CanLuo/long_reads_project/Benchmarks/HG002_SVs_Tier1_v0.6_chr_noXY.bed


mkdir $work_dir
output_dir=$work_dir/Truvari_${aligner}_aligned_results_p${p}_P${P}_dist${dist}_S${S}_O${O}_chr${chr_num}
rm -r $output_dir

mkdir $output_dir

python /data/maiziezhou_lab/CanLuo/long_reads_project/bin/vcf_filter.py -v ${input_dir}/${prefix}.vcf -o_dir $work_dir/  --chrs chr${chr_num}

vcf-sort $work_dir/${prefix}_DEL_INS_noXY.vcf > $work_dir/${prefix}_DEL_INS_noXY_sorted.vcf
bgzip -c $work_dir/${prefix}_DEL_INS_noXY_sorted.vcf > $work_dir/${prefix}_DEL_INS_noXY_sorted.vcf.gz
tabix -p vcf $work_dir/${prefix}_DEL_INS_noXY_sorted.vcf.gz

vcf-sort $work_dir/${prefix}_DEL_noXY.vcf > $work_dir/${prefix}_DEL_noXY_sorted.vcf
bgzip -c $work_dir/${prefix}_DEL_noXY_sorted.vcf > $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz
tabix -p vcf $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz

vcf-sort $work_dir/${prefix}_INS_noXY.vcf > $work_dir/${prefix}_INS_noXY_sorted.vcf
bgzip -c $work_dir/${prefix}_INS_noXY_sorted.vcf > $work_dir/${prefix}_INS_noXY_sorted.vcf.gz
tabix -p vcf $work_dir/${prefix}_INS_noXY_sorted.vcf.gz


#####################################################################################################################
#####################################################################################################################
#sizemin 50


truvari bench -b /data/maiziezhou_lab/CanLuo/long_reads_project/Benchmarks/HG002_SVs_Tier1_v0.6_chr${chr_num}_del.vcf.gz -c $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o $output_dir/DEL_50_ --includebed ${bed} -p $p -P $P -r $dist --passonly --sizemin 50 --prog -S ${S} -O ${O}

truvari bench -b /data/maiziezhou_lab/CanLuo/long_reads_project/Benchmarks/HG002_SVs_Tier1_v0.6_chr${chr_num}_ins.vcf.gz  -c $work_dir/${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o $output_dir/INS_50_  --includebed ${bed} -p $p -P $P -r $dist --passonly --sizemin 50 --prog -S ${S} -O ${O}

#####################################################################################################################
#       Parse result
#####################################################################################################################

#### parse excel
python3 /data/maiziezhou_lab/CanLuo/long_reads_project/bin/truvari_result_parser_indel.py -i $output_dir/
