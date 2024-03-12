#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=60G
#$ -l h_rt=20:00:00 
#$ -o logs/BGEN2BED_info_chr6.out
#$ -e logs/BGEN2BED_info_chr6.err
#$ -N bgen2bed_info_chr6

. /etc/profile.d/modules.sh

module load igmm/apps/qctool/2.0.8

qctool -g /exports/igmm/eddie/ponting-lab/breeshey/data/genomicc/imputed/genomic-isaric/GenOMICC_ISARIC_chr6.bgen -s /exports/igmm/eddie/ponting-lab/breeshey/data/genomicc/imputed/genomic-isaric/GenOMICC_ISARIC_chr6.sample -excl-rsids output/chr6_rsids2exclude_info_score0.6.txt -og "output/bed_files/GenOMICC_ISARIC_info_score_0.6_chr6.bed"

