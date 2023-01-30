#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=6


bam_dir=/data/nym/20210339_LGG/hisat2_res

cd /data/nym/20210339_LGG/all_call_peak/mergeBed_res
bam_ls=`ls ${bam_dir}/*.bam`
merged_bed=/data/nym/20210339_LGG/all_call_peak/mergeBed_res/mergeBed_res.txt
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,"",$4}' $merged_bed >data_for_counts.bed
multiBamSummary BED-file --BED data_for_counts.bed --bamfiles ${bam_ls} -o merged_peak_counts.npz --outRawCounts merged_peak_counts.tab





