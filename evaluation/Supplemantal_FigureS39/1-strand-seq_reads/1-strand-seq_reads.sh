# NOTE: This is just for changing the suffix to .fastq.gz and shorten then file name for downstream processing. File content not changed
#/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/Strand-seq/HG002_UCSC-hpp-v1_12162019/2019-12-16-HWVTJAFXY/
strandseq_download_dir=$1


for i in `ls ${strandseq_download_dir}`; do j=${i##HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1}; j=${j%%.txt.gz}.fastq.gz; cp ${strandseq_download_dir}/${i} $j; done

ls *.fastq.gz | cut -d'_' -f1 | sort -u > cells.list

# convert to BGZIP for biopyhon indexing 
for i in *.gz; do zcat ${i} | bgzip -c > ${i%%.gz}.bgz; done

# overwrite the copied gz with bgz to save place
for i in *.bgz; do mv ${i} ${i%%.bgz}.gz; done

