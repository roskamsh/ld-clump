import pandas as pd
import numpy as np
import collections

all_eqtls = pd.read_csv("input/all_genomicc_associated_eQTLs_hg38.csv")
tf = "ESR1"
bed_snps = pd.read_csv("output/bed_files/GenOMICC_ISARIC_info_score_0.6_chr6.bim", delimiter = " ", header=None)
bed_snps.columns = ["CHR","SNP_ID","GENETIC_DISTANCE","POS","REF","ALT"]

# Filter all_eqtls for only TF-of-interest
eqtls = all_eqtls[all_eqtls.GeneSymbol==tf].copy()

# Loop through eqtls and identify where in the bed_snps there contains a match
# Then match this SNP_ID to the row in eqtls
nrow = eqtls.shape[0]
new_ids = [0]*nrow
ref_alleles = [0]*nrow
alt_alleles = [0]*nrow
for i, snp in enumerate(eqtls.SNP.values):
    pos = eqtls.iloc[i]['SNPPos']
    allele1 = eqtls.iloc[i]['AssessedAllele']
    allele2 = eqtls.iloc[i]['OtherAllele']

    row = bed_snps.index[
                (bed_snps['POS']==pos)
                ]
    sub = bed_snps.loc[row].copy()

    if sub.shape[0]==0:
        new_ids[i] = pd.NA
        ref_alleles[i] = pd.NA
        alt_alleles[i] = pd.NA

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
            ref_alleles[i] = ref_allele
            alt_alleles[i] = alt_allele

# Add to eqtls
eqtls['BGEN_ID'] = new_ids
eqtls['A1'] = ref_allele
eqtls['A2'] = alt_allele

# Remove missing values
eqtls = eqtls.dropna(subset = ["BGEN_ID"])

# Reformat to match what is expected in PLINK for LD-clumping
eqtls = eqtls[["BGEN_ID","SNP","SNPChr","SNPPos","A1","A2","Pvalue"]]
eqtls.columns = ["SNP","RSID","CHR","BP","A1","A2","P"]

eqtls.to_csv(f"input/{tf}_ciseQTLs_hg38.assoc", index = False, sep = "\t")