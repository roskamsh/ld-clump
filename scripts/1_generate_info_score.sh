#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=40G
#$ -l h_rt=10:00:00 
#$ -o logs/SNPSTATS_chr6.out
#$ -e logs/SNPSTATS_chr6.err
#$ -N snpstats_chr6

. /etc/profile.d/modules.sh

module load igmm/apps/qctool/2.0.8

qctool -g /exports/igmm/eddie/ponting-lab/breeshey/data/genomicc/imputed/genomic-isaric/GenOMICC_ISARIC_chr6.bgen -s /exports/igmm/eddie/ponting-lab/breeshey/data/genomicc/imputed/genomic-isaric/GenOMICC_ISARIC_chr6.sample -snp-stats -osnp output/chr6_imputation_stats.snpstats
