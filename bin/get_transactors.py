import pandas as pd
from argparse import ArgumentParser

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--eqtls", type = str, help = "Full input matrix, usually eQTLGen, and must match format expected."
    )
    parser.add_argument(
        "--transactors", type = str, help = """Full input matrix, for additional transactor QTLs, and must match format expected. 
        This means it must contain the following columns: []"SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","GeneSymbol","Pvalue"]"""
    )
    parser.add_argument(
        "--tf", type = str, help = "Gene Symbol for the TF you want to find the CHR for."
    )
    return parser.parse_args()

if __name__ == '__main__':
    eqtls = get_input().eqtls
    transactors = get_input().transactors
    tf = get_input().tf

    eqtls = pd.read_csv(eqtls, delimiter="\t")
    eqtls = eqtls[eqtls.GeneSymbol == tf]
    eqtls = eqtls[["SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","GeneSymbol","Pvalue"]]

    # If duplicate SNPs (GTEx), keep only 1, prioritizing lower p-value
    eqtls = eqtls.sort_values(by="Pvalue", ascending=True).drop_duplicates(subset="SNP", keep="first")
    eqtls = eqtls.reset_index(drop=True)

    if transactors == "NO_ADDITIONAL_TRANSACTORS":
        data = eqtls
    else:
        transactors = pd.read_csv(transactors, delimiter="\t")
        transactors = transactors[transactors.GeneSymbol == tf]
        transactors = transactors[["SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","GeneSymbol","Pvalue"]]
        data = pd.concat([eqtls,transactors], ignore_index=True)

    # Check if TF is present in eQTLGen database
    if data.shape[0] == 0:
        print(f"No transactor QTLs provided or available in eQTLGen for {tf}")
    else:
        data.to_csv(f"{tf}_transactor_QTLs.csv", index = False)
