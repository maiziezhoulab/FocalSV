

import pandas as pd
from subprocess import Popen
from match_sv import match_union_ins
import os
def vcf_to_bed(input_vcf, output_bed, flank=100, svlen_threshold=30):
	with open(input_vcf, 'r') as vcf, open(output_bed, 'w') as bed:
		for line in vcf:
			if line.startswith('#'):
				continue
			cols = line.strip().split('\t')
			chrom, pos, info = cols[0], int(cols[1]), cols[7]
			if chrom.startswith('chr'):
				chrnum = chrom[3:]
				if chrnum.isdigit() and 1 <= int(chrnum) <= 22:
					info_fields = {key: val for key, val in (field.split('=') for field in info.split(';') if '=' in field)}
					svlen = int(info_fields.get('SVLEN', 0))
					if abs(svlen) >= svlen_threshold:
						start = max(pos - flank, 0)
						end = pos + flank
						bed.write(f"{chrom}\t{start}\t{end}\n")

def filter_vcf_by_bed_del(invcf, bedfile, ):
	invcf = invcf.replace(".vcf",'')
	cmd = f'''
	grep '#\|DEL'  {invcf}.vcf |bgzip > {invcf}_del.vcf.gz
	tabix -p vcf {invcf}_del.vcf.gz 
	bcftools view -R {bedfile} -o {invcf}_del_filter.vcf  {invcf}_del.vcf.gz 
	'''
	print(cmd)
	Popen(cmd, shell= True).wait()

def final_process_ont(infile,reads_draft_vcf, outfile):
	ins_vcf = infile.replace(".vcf",'_ins_union.vcf')
	del_vcf = infile.replace(".vcf",'_del_filter.vcf')
	bed_file = reads_draft_vcf.replace(".vcf","_chr1_22_gte30.bed")
	match_union_ins(infile,reads_draft_vcf,ins_vcf)

	vcf_to_bed(reads_draft_vcf, bed_file )
	filter_vcf_by_bed_del(infile, bed_file)

	cmd = f"(cat {ins_vcf}; grep -v '^#' {del_vcf}) | vcf-sort > {outfile}"
	print(cmd)
	Popen(cmd, shell= True).wait()

		
		
		
