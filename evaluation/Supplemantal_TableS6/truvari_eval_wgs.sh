#####################################################################################################################
#
#                       truvari evaluation
#
#####################################################################################################################
input_dir=$1
work_dir=$2
prefix=$3
dist=500
p=0.5
P=0.5
S=30
O=0.01
bed=$4
ref=$5
bench_dir=$6
dipcall=${7}



mkdir -p $work_dir
output_dir=$work_dir/Truvari_results_p${p}_P${P}_dist${dist}_S${S}_O${O}
rm -r $output_dir

mkdir -p $output_dir

python ./vcf_filter.py -v ${input_dir}/${prefix}.vcf -o_dir $work_dir/  $dipcall

vcf-sort $work_dir/${prefix}_DEL_noXY.vcf > $work_dir/${prefix}_DEL_noXY_sorted.vcf
bgzip -c $work_dir/${prefix}_DEL_noXY_sorted.vcf > $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz
tabix -p vcf $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz

vcf-sort $work_dir/${prefix}_INS_noXY.vcf > $work_dir/${prefix}_INS_noXY_sorted.vcf
bgzip -c $work_dir/${prefix}_INS_noXY_sorted.vcf > $work_dir/${prefix}_INS_noXY_sorted.vcf.gz
tabix -p vcf $work_dir/${prefix}_INS_noXY_sorted.vcf.gz


#####################################################################################################################
#     Evaluation
#####################################################################################################################


truvari bench -b ${bench_dir}/HG002_SVs_Tier1_v0.6_chr_ins_noXY.vcf.gz  -c $work_dir/${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o $output_dir/INS_50_  --includebed ${bed} -p $p -P $P -r $dist --passonly --sizemin 50  -S ${S}  -O ${O}

### DEL
truvari bench -b ${bench_dir}/HG002_SVs_Tier1_v0.6_chr_del_noXY.vcf.gz -c $work_dir/${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o $output_dir/DEL_50_ --includebed ${bed} -p $p -P $P -r $dist --passonly --sizemin 50    -S ${S}  -O ${O}

