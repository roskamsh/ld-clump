import pandas as pd
import numpy as np
from argparse import ArgumentParser

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--snps", type = str, help = "SNPs from eQTLGen for which you want to generate an .assoc file for LD clumping"
    )
    parser.add_argument(
        "--tf", type = str, help = "Gene Symbol for the TF you wish to generate this file for"
    )
    parser.add_argument(
        "--prefix", type = str, help = "bed file prefix for chromosome which contains SNPs you will run LD clumping on"
    )
    return parser.parse_args()

if __name__ == '__main__':
    snps = get_input().snps
    tf = get_input().tf
    prefix = get_input().prefix
    
    bim_file = prefix + ".bim"
    all_eqtls = pd.read_csv(snps, delimiter = "\t")
    bed_snps = pd.read_csv(bim_file, delimiter = " ", header=None)
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

        # I think a bunch of the eQTLs are not in this bed file, so we need to catch these cases
        # This is occuring when the eQTLs fall in the "exclusion regions"
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

    eqtls['BGEN_ID'] = new_ids
    eqtls['A1'] = ref_alleles
    eqtls['A2'] = alt_alleles

    # Remove missing values
    eqtls = eqtls.dropna(subset = ["BGEN_ID"])

    # Reformat to match what is expected in PLINK for LD-clumping
    eqtls = eqtls[["BGEN_ID","SNP","SNPChr","SNPPos","A1","A2","Pvalue"]]
    eqtls.columns = ["SNP","RSID","CHR","BP","A1","A2","P"]
    
    if eqtls.empty:
        print("cis-eQTLs lie within exclusion regions, no SNPs found in processed .bim file. Exiting...")
        exit(0)
    else:
        eqtls.to_csv(f"{tf}_ciseQTLs_hg38.assoc", index = False, sep = "\t")