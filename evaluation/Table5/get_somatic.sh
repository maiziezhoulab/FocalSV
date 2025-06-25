#### input #####

vcffile_normal=${1}
vcffile_cancer=${2}
outdir=${3}/
pass_only=${4}
SURVIVOR=${5}

mkdir -p ${outdir}
prefix=variants
prefix_normal=${prefix}_normal
prefix_cancer=${prefix}_cancer


#### pass only?

if [ "$pass_only" -eq 1 ]; then
    awk '$0 ~ /^#/ || $7 == "PASS"' ${vcffile_normal} >  ${outdir}${prefix_normal}_passonly.vcf
    awk '$0 ~ /^#/ || $7 == "PASS"' ${vcffile_cancer} >  ${outdir}${prefix_cancer}_passonly.vcf
    vcffile_normal=${outdir}${prefix_normal}_passonly.vcf
    vcffile_cancer=${outdir}${prefix_cancer}_passonly.vcf
fi


######### split normal vcf #########

grep "#" ${vcffile_normal}>${outdir}header
cat ${vcffile_normal} | grep -w "SVTYPE=BND\|SVTYPE=TRA" > ${outdir}${prefix_normal}_BND.body
cat ${outdir}header ${outdir}${prefix_normal}_BND.body > ${outdir}${prefix_normal}_BND.vcf
cat ${vcffile_normal} | grep -w -v "SVTYPE=BND\|SVTYPE=TRA\|SVTYPE=INS\|SVTYPE=DEL" > ${outdir}${prefix_normal}_noBND.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_normal}_noBND.vcf NA 50 100 0 -1  ${outdir}${prefix_normal}_noBND_50_100.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_normal}_noBND.vcf NA 101 500 0 -1  ${outdir}${prefix_normal}_noBND_101_500.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_normal}_noBND.vcf NA 501 1000 0 -1  ${outdir}${prefix_normal}_noBND_501_1000.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_normal}_noBND.vcf NA 1001 30000 0 -1  ${outdir}${prefix_normal}_noBND_1001_30000.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_normal}_noBND.vcf NA   30000 -1 0 -1  ${outdir}${prefix_normal}_noBND_30000_.vcf


######### split cancer vcf #########

grep "#" ${vcffile_cancer}>${outdir}header
cat ${vcffile_cancer} | grep -w "SVTYPE=BND\|SVTYPE=TRA" > ${outdir}${prefix_cancer}_BND.body
cat ${outdir}header ${outdir}${prefix_cancer}_BND.body > ${outdir}${prefix_cancer}_BND.vcf
cat ${vcffile_cancer} | grep -w -v "SVTYPE=BND\|SVTYPE=TRA\|SVTYPE=INS\|SVTYPE=DEL" > ${outdir}${prefix_cancer}_noBND.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_cancer}_noBND.vcf NA 50 100 0 -1  ${outdir}${prefix_cancer}_noBND_50_100.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_cancer}_noBND.vcf NA 101 500 0 -1  ${outdir}${prefix_cancer}_noBND_101_500.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_cancer}_noBND.vcf NA 501 1000 0 -1  ${outdir}${prefix_cancer}_noBND_501_1000.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_cancer}_noBND.vcf NA 1001 30000 0 -1  ${outdir}${prefix_cancer}_noBND_1001_30000.vcf
${SURVIVOR}/Debug/SURVIVOR filter ${outdir}${prefix_cancer}_noBND.vcf NA   30000 -1 0 -1  ${outdir}${prefix_cancer}_noBND_30000_.vcf


######### merge normal and cancer #########

rm ${outdir}vcflist_50_100.txt
echo ${outdir}${prefix_normal}_noBND_50_100.vcf >> ${outdir}vcflist_50_100.txt
echo ${outdir}${prefix_cancer}_noBND_50_100.vcf >> ${outdir}vcflist_50_100.txt
${SURVIVOR}/Debug/SURVIVOR merge ${outdir}vcflist_50_100.txt 50 1 1 0 0 50 ${outdir}${prefix}_noBND_50_100_merged.vcf


rm ${outdir}vcflist_101_500.txt
echo ${outdir}${prefix_normal}_noBND_101_500.vcf >> ${outdir}vcflist_101_500.txt
echo ${outdir}${prefix_cancer}_noBND_101_500.vcf >> ${outdir}vcflist_101_500.txt
${SURVIVOR}/Debug/SURVIVOR merge ${outdir}vcflist_101_500.txt 101 1 1 0 0 101 ${outdir}${prefix}_noBND_101_500_merged.vcf


rm ${outdir}vcflist_501_1000.txt
echo ${outdir}${prefix_normal}_noBND_501_1000.vcf >> ${outdir}vcflist_501_1000.txt
echo ${outdir}${prefix_cancer}_noBND_501_1000.vcf >> ${outdir}vcflist_501_1000.txt
${SURVIVOR}/Debug/SURVIVOR merge ${outdir}vcflist_501_1000.txt 501 1 1 0 0 501 ${outdir}${prefix}_noBND_501_1000_merged.vcf


rm ${outdir}vcflist_1001_30000.txt
echo ${outdir}${prefix_normal}_noBND_1001_30000.vcf >> ${outdir}vcflist_1001_30000.txt
echo ${outdir}${prefix_cancer}_noBND_1001_30000.vcf >> ${outdir}vcflist_1001_30000.txt
${SURVIVOR}/Debug/SURVIVOR merge ${outdir}vcflist_1001_30000.txt 1001 1 1 0 0 1001 ${outdir}${prefix}_noBND_1001_30000_merged.vcf

rm ${outdir}vcflist_30000_.txt
echo ${outdir}${prefix_normal}_noBND_30000_.vcf >> ${outdir}vcflist_30000_.txt
echo ${outdir}${prefix_cancer}_noBND_30000_.vcf >> ${outdir}vcflist_30000_.txt
${SURVIVOR}/Debug/SURVIVOR merge ${outdir}vcflist_30000_.txt 10000 1 1 0 0 10000 ${outdir}${prefix}_noBND_30000_merged.vcf


rm ${outdir}vcflist_BND.txt
echo ${outdir}${prefix_normal}_BND.vcf >> ${outdir}vcflist_BND.txt
echo ${outdir}${prefix_cancer}_BND.vcf >> ${outdir}vcflist_BND.txt
${SURVIVOR}/Debug/SURVIVOR merge  ${outdir}vcflist_BND.txt 1000 1 1 0 0 0 ${outdir}${prefix}_BND_merged.vcf


######### merge all size #########


cat ${outdir}${prefix}_noBND_50_100_merged.vcf \
    ${outdir}${prefix}_noBND_101_500_merged.vcf \
    ${outdir}${prefix}_noBND_501_1000_merged.vcf \
    ${outdir}${prefix}_noBND_1001_30000_merged.vcf \
    ${outdir}${prefix}_noBND_30000_merged.vcf \
     |grep -v "#"> ${outdir}${prefix}_noBND_merged.body
cat ${outdir}header ${outdir}${prefix}_noBND_merged.body > ${outdir}${prefix}_noBND_merged.vcf

grep "SUPP_VEC=01" ${outdir}${prefix}_noBND_merged.vcf > ${outdir}${prefix}_noBND_somatic.body
cat ${outdir}header ${outdir}${prefix}_noBND_somatic.body > ${outdir}${prefix}_noBND_somatic.vcf

grep "SUPP_VEC=01"  ${outdir}${prefix}_BND_merged.vcf > ${outdir}${prefix}_BND_somatic.body
cat ${outdir}header ${outdir}${prefix}_BND_somatic.body > ${outdir}${prefix}_BND_somatic.vcf 

cat ${outdir}header ${outdir}${prefix}_noBND_somatic.body  ${outdir}${prefix}_BND_somatic.body > ${outdir}${prefix}_all_somatic.vcf



### clean file
# find ${outdir} -type f ! -name '*somatic*.vcf' -delete
find "${outdir}" -type f ! \( -name '*somatic*.vcf' -o -name '*passonly.vcf' \) -delete
