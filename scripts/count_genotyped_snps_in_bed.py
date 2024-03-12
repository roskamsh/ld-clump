import numpy as np
import pandas as pd

if __name__ == '__main__':
    snp_list = pd.read_csv("input/genotyped_snps_converted_hg19_hg38.csv", delimiter = "\t")
    snp_list = snp_list.convert_dtypes()
    snp_list['POS_grch38'] = snp_list['POS_grch38'].astype(str)
    snp_list["BGEN_ID_grch38"] = snp_list[["CHR","POS_grch38","REF","ALT"]].apply(lambda row: ':'.join(row), axis=1) 

    # Remove NA values
    snp_list_clean = snp_list.dropna(subset=['RSID'])

    ## Now read in the .bim files into pandas and count the number of genotyped SNPs in each .bim
    imp0_8 = pd.read_csv("output/bed_files/GenOMICC_ISARIC_imp0.8_chr6.bim", delimiter = " ")
    imp0_9 = pd.read_csv("output/bed_files/GenOMICC_ISARIC_imp0.9_chr6.bim", delimiter = " ")

    imp0_8.columns = ["CHR","BGEN_ID","GENETIC_DISTANCE","POS","REF","ALT"]
    imp0_9.columns = ["CHR","BGEN_ID","GENETIC_DISTANCE","POS","REF","ALT"]

    # loop through and find matches
    matches_imp0_8 = []
    matches_imp0_9 = []
    for index, row in snp_list_clean.iterrows():
        print(index, "/", snp_list_clean.shape[0])
        if row['BGEN_ID_grch38'] in imp0_8['BGEN_ID'].values:
            matches_imp0_8.append(row)
            print(row["BGEN_ID_grch38"])
        if row['BGEN_ID_grch38'] in imp0_9['BGEN_ID'].values:
            matches_imp0_9.append(row) 
            print(row["BGEN_ID_grch38"])

    matches_imp0_8_df = pd.DataFrame(matches_imp0_8)
    matches_imp0_9_df = pd.DataFrame(matches_imp0_9)

    matches_imp0_8_df.to_csv("output/genotyped_snps_that_match_imp0.8_bed.csv", index = False)
    matches_imp0_9_df.to_csv("output/genotyped_snps_that_match_imp0.9_bed.csv", index = False)
