#=========== need edit below before use
ref=<reference> 
vcf_dir=vcf_dir>
bench_dir=<bench_dir>
#=========== need edit above before use

prefix=TOOL
vcf=${vcf_dir}/TOOL.vcf
bench=${bench_dir}/HG002_SVs_Tier1_v0.6_chr_noXY.vcf.gz
bench_del=${bench_dir}/HG002_SVs_Tier1_v0.6_chr_del_noXY.vcf.gz
bench_ins=${bench_dir}/HG002_SVs_Tier1_v0.6_chr_ins_noXY.vcf.gz
out_dir=${prefix}_Truvari

p_default=0.5
P_default=0.5
r_default=500
O_default=0.01

rm -r ${out_dir}

mkdir ${out_dir}


python ./vcf_filter.py -v ${vcf} -o_dir . -p ${prefix}



vcf-sort ${prefix}_DEL_noXY.vcf > ${prefix}_DEL_noXY_sorted.vcf
bgzip -c ${prefix}_DEL_noXY_sorted.vcf > ${prefix}_DEL_noXY_sorted.vcf.gz
tabix -p vcf ${prefix}_DEL_noXY_sorted.vcf.gz

vcf-sort ${prefix}_INS_noXY.vcf > ${prefix}_INS_noXY_sorted.vcf
bgzip -c ${prefix}_INS_noXY_sorted.vcf > ${prefix}_INS_noXY_sorted.vcf.gz
tabix -p vcf ${prefix}_INS_noXY_sorted.vcf.gz


#####################################################################################################################
#####################################################################################################################

for O in 0.01 `seq 0.1 0.1 1.0`
do
	truvari bench -b ${bench_del} -c ${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/DEL_O${O}_50_ --includebed ${bed} -O ${O} -p ${p_default} -P ${P_default} -r ${r_default} --passonly --sizemin 50
	truvari bench -b ${bench_ins} -c ${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/INS_O${O}_50_  --includebed ${bed} -O ${O} -p ${p_default} -P ${P_default} -r ${r_default} --passonly --sizemin 50
done

for p in `seq 0 0.1 1.0`
do
	truvari bench -b ${bench_del} -c ${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/DEL_p${p}_50_ --includebed ${bed} -O ${O_default} -p ${p} -P ${P_default} -r ${r_default} --passonly --sizemin 50
	truvari bench -b ${bench_ins} -c ${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/INS_p${p}_50_  --includebed ${bed} -O ${O_default} -p ${p} -P ${P_default} -r ${r_default} --passonly --sizemin 50
done

for P in `seq 0 0.1 1.0`
do
	truvari bench -b ${bench_del} -c ${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/DEL_P${P}_50_ --includebed ${bed} -O ${O_default} -p ${p_default} -P ${P} -r ${r_default} --passonly --sizemin 50
	truvari bench -b ${bench_ins} -c ${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/INS_P${P}_50_  --includebed ${bed} -O ${O_default} -p ${p_default} -P ${P} -r ${r_default} --passonly --sizemin 50
done

#NOTE: for r loop, O is fixed to 0
for r in `seq 0 100 1000`
do
	truvari bench -b ${bench_del} -c ${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/DEL_r${r}_50_ --includebed ${bed} -O 0 -p ${p_default} -P ${P_default} -r ${r} --passonly --sizemin 50
	truvari bench -b ${bench_ins} -c ${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/INS_r${r}_50_  --includebed ${bed} -O 0 -p ${p_default} -P ${P_default} -r ${r} --passonly --sizemin 50
done

#Keep results for O=0.01
for r in `seq 0 100 1000`
do
	truvari bench -b ${bench_del} -c ${prefix}_DEL_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/DEL.O0.01_r${r}_50_ --includebed ${bed} -O ${O_default} -p ${p_default} -P ${P_default} -r ${r} --passonly --sizemin 50
	truvari bench -b ${bench_ins} -c ${prefix}_INS_noXY_sorted.vcf.gz -f ${ref} -o ${out_dir}/INS.O0.01_r${r}_50_  --includebed ${bed} -O ${O_default} -p ${p_default} -P ${P_default} -r ${r} --passonly --sizemin 50
done


#clean up
#rm ${prefix}*noXY*
