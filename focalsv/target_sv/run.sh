awk 'BEGIN{OFS="\t"} !/^#/ {
len_ref = length($4);
len_alt = length($5);
len_diff = (len_ref > len_alt) ? (len_ref - len_alt) : 0;
start = $2 - 50000;
end = $2 + len_diff + 50000;
if (start < 0) start = 0;
	print $1, start, end
}' bench |vcf-sort> HG002_SV_rich_regions.bed
grep chr21 HG002_SV_rich_regions.bed > HG002_SV_rich_regions_chr21.bed

