#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

cd /data/nym/20210339_LGG/all_call_peak/mergeBed_res
peak_dir=/data/nym/20210339_LGG/all_call_peak
#for i in `ls $peak_dir/*.xls`
#do
#	base_name=`basename ${i} .xls`
#	sample_name=`echo ${base_name%%_*}`
#	cut -f 1-4,8,9 ${i}|awk -v name=$sample_name 'BEGIN{OFS="\t"}{print $0,name}' > ${sample_name}.bed
#done

rm peak_gtf_intersect.bed
rm mergeBed_res.txt
touch mergeBed_res.txt
touch peak_gtf_intersect.bed
peak_ls=`ls *[NT].bed`
cat $peak_ls|sort -k1,1 -k2,2n >peak_cated.bed
cut -f 4 peak_cated.bed |sort|uniq >gene.txt 
num=`less -S gene.txt | wc -l`
for i in `seq 1 $num`
do
	gene=`head -$i gene.txt | tail -1`;grep -r "$gene" peak_cated.bed >tmp.bed 
	~/tools/bedtools2-2.26.0/bin/mergeBed -i tmp.bed -s -c 4,6,7 -o distinct,distinct,distinct >>mergeBed_res.txt
done


num=`less -S gene.txt | wc -l`
ref_file=/data/nym/190715-AD-human/diff_3sample/add_exon/merge_GRCh38.85_2.bed
for i in `seq 1 $num`
do
	gene=`head -$i gene.txt | tail -1`
	if [ `grep -c "$gene" $ref_file` -ne '0' ];then
		grep -r "$gene" /data/nym/190715-AD-human/diff_3sample/add_exon/merge_GRCh38.85_2.bed > tmp_ref_gtf.bed
		grep -r "$gene" mergeBed_res.txt > tmp_peak.bed
		#intersect
		~/tools/bedtools2-2.26.0/bin/intersectBed -s -wo -a tmp_peak.bed -b tmp_ref_gtf.bed >> peak_gtf_intersect.bed
	fi
done
perl peak_block.pl peak_gtf_intersect.bed peak_out.bed
~/tools/bedtools2-2.26.0/bin/sortBed -i peak_out.bed > peak_out_sort.bed




