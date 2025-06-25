
libs=("hifi_l1" "clr_l1" "ont_l1")
ranks=("top10" "bottom10" "top50" "bottom50")

highc_bed="$1/HG002_SVs_Tier1_v0.6_chr_noXY.bed"

outdir="./intersect_results"
rank_bed_dir=./Phasing_Rank_BED
mkdir -p $outdir


export highc_bed outdir rank_bed_dir
parallel --will-cite -j 6 '
    lib={1}
    rank={2}
    phase_bed=${rank_bed_dir}/${lib}_${rank}_regions.bed
    out_bed=${outdir}/${lib}_${rank}_highconf_overlap.bed
    bedtools intersect -a $highc_bed -b $phase_bed > $out_bed
' ::: "${libs[@]}" ::: "${ranks[@]}"
