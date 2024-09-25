conda activate plink

bqtl_list="bqtl_list_ESR1.txt"
eqtl_list="eqtl_list_ESR1.txt"

merged_bed="work/50/ce2002899b431aa351b547ea6335f4/merged_dataset.bed"

# Unfortunately, going to have to restrict the BED file to only the SNPs we are interested in
cat $bqtl_list <(echo) $eqtl_list > snp_list.txt 
bed_prefix=$( echo $merged_bed | sed 's/.bed//g' )

plink --bfile $bed_prefix --extract snp_list.txt --make-bed --out extracted_data

plink --bfile extracted_data --r2 inter-chr --out ld_results
