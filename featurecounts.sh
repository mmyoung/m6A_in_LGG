#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

cd /data/nym/20210339_LGG/all_call_peak/src
R CMD BATCH featurecounts.R

