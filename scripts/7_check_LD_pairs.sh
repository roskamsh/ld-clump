#!/usr/bin/env bash
conda activate plink2

# Test code for checking LD between two SNPs (bQTL first, eQTL second)
var1="chr6:89985527:A:G"
var2="chr6:151689119:A:G"
chr1="chr6"
chr2="chr6"

if [ $chr1 == $chr2 ]; then
    plink2 --bfile work/cc/d7ef07e9bb8102c2373da7bab8c473/GenOMICC_ISARIC_chr_info_score_0.6_chr6 --ld chr6:89985527:A:G chr6:151689119:A:G > plink2_out.log
    r2=$( grep "r^2" plink2_out.log | cut -f 5 -d " " )

    if [ $r2 >= 0.8 ]; then
        echo "$var1 and $var2 are in LD with one another, exclude this interaction"
    else
        echo "$var1 and $var2 are not in LD with one another, keep this interaction"
else
    echo "$var1 and $var2 are on different chromosomes, keep this interaction"
fi

