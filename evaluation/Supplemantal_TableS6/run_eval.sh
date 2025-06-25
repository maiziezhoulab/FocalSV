libs=("Hifi_L1" "CLR_L1" "ONT_L1")
ranks=("top10" "bottom10" "top50" "bottom50")


script=./truvari_eval_wgs.sh
indir=$1


for lib in libs
do
for rank in ranks
do
prefix=${lib}_FocalSV-auto
lib_lower=$(echo $lib | tr "[:upper:]" "[:lower:]")
outdir=./eval_${lib_lower}_${rank}
bed=./intersect_results/${lib_lower}_${rank}_highconf_overlap.bed
$script $indir $outdir $prefix $bed $2 $3
done
done
