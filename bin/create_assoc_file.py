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
    transactors = pd.read_csv(snps)
    bed_snps = pd.read_csv(bim_file, delim_whitespace=True, header=None)
    bed_snps.columns = ["CHR","SNP_ID","GENETIC_DISTANCE","POS","A1","A2"]

    # Ensure only TF-of-interest
    transactors = transactors[transactors.GeneSymbol==tf].copy()

    # Loop through eqtls and identify where in the bed_snps there contains a match
    # Then match this SNP_ID to the row in eqtls
    nrow = transactors.shape[0]
    new_ids = [0]*nrow
    first_alleles = [0]*nrow
    second_alleles = [0]*nrow

    for i, snp in enumerate(transactors.SNP.values):
        pos = transactors.iloc[i]['SNPPos']
        allele1 = transactors.iloc[i]['AssessedAllele']
        allele2 = transactors.iloc[i]['OtherAllele']

        # Some transactors are not in this bed file, so we need to catch these cases
        # This is occuring when the eQTLs fall in the "exclusion regions"
        row = bed_snps.index[
                    (bed_snps['POS']==pos)
                    ]
        sub = bed_snps.loc[row].copy()

        if sub.shape[0]==0:
            new_ids[i] = pd.NA
            first_alleles[i] = pd.NA
            second_alleles[i] = pd.NA

        # Loop through in case there are multiple SNPs at that position
        for j, matched_snp in enumerate(sub.SNP_ID.values):
            first_allele = sub.iloc[j]['A1']
            second_allele = sub.iloc[j]['A2']
            # Check for a matching of alleles in either order
            bed_snp_alleles = [first_allele, second_allele]
            eqtlgen_alleles = [allele1, allele2]

            # Checks if elements match, regardless of ordering
            if set(bed_snp_alleles) == set(eqtlgen_alleles):
                new_ids[i] = matched_snp
                first_alleles[i] = first_allele
                second_alleles[i] = second_allele

    transactors['BGEN_ID'] = new_ids
    transactors['A1'] = first_alleles
    transactors['A2'] = second_alleles

    # Remove missing values
    transactors = transactors.dropna(subset = ["BGEN_ID"])

    # Reformat to match what is expected in PLINK for LD-clumping
    transactors = transactors[["BGEN_ID","SNP","SNPChr","SNPPos","A1","A2","Pvalue"]]
    transactors.columns = ["SNP","RSID","CHR","BP","A1","A2","P"]
    
    if transactors.empty:
        print("cis-eQTLs lie within exclusion regions, no SNPs found in processed .bim file. Exiting...")
        exit(0)
    else:
        transactors.to_csv(f"{tf}_ciseQTLs_hg38.assoc", index = False, sep = "\t")