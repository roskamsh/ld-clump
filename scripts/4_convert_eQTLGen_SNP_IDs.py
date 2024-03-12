import pandas as pd
import numpy as np
from collections import set

eqtls = pd.read_csv("input/ESR1_cis_eQTLs_hg38.assoc", delimiter="\t")
bed_snps = pd.read_csv("output/bed_files/GenOMICC_ISARIC_info_score_0.6_chr6.bim", delimiter = " ", header=None)
bed_snps.columns = ["CHR","SNP_ID","GENETIC_DISTANCE","POS","REF","ALT"]

# Loop through eqtls and identify where in the bed_snps there contains a match
# Then match this SNP_ID to the row in eqtls
new_ids = [0]*eqtls.shape[0]
for i, snp in enumerate(eqtls.SNP.values):
    pos = eqtls.iloc[i]['BP']
    allele1 = eqtls.iloc[i]['A1']
    allele2 = eqtls.iloc[i]['A2']

    row = bed_snps.index[
                (bed_snps['POS']==pos)
                ]
    sub = bed_snps.loc[row].copy()

    # Loop through in case there are multiple SNPs at that position
    for j, matched_snp in enumerate(sub.SNP_ID.values):
        ref_allele = sub.iloc[j]['REF']
        alt_allele = sub.iloc[j]['ALT']
        # Check for a matching of alleles in either order
        bed_snp_alleles = [ref_allele, alt_allele]
        eqtlgen_alleles = [allele1, allele2]

        # Checks if elements match, regardless of ordering
        if set(bed_snp_alleles) == set(eqtlgen_alleles):
            new_ids[i] = matched_snp

# Add to eqtls
eqtls['SNP'] = new_ids

eqtls.to_csv("input/ESR1_cis_eQTLs_hg38_formatted_IDs.assoc", index = False, sep = "\t")