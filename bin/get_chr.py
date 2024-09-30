import pandas as pd
from argparse import ArgumentParser

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--data", type = str, help = "eQTLGen data, full matrix"
    )
    parser.add_argument(
        "--tf", type = str, help = "Gene Symbol for the TF you want to find the CHR for"
    )
    return parser.parse_args()

if __name__ == '__main__':
    data = get_input().data
    tf = get_input().tf

    eqtl_data = pd.read_csv(data, delimiter="\t")
    eqtl_data = eqtl_data[eqtl_data.GeneSymbol == tf]

    # Check if TF is present in eQTLGen database
    if eqtl_data.shape[0] == 0:
        print("Not found in eQTLGen")
    else:
        chr = list(eqtl_data.SNPChr.unique())[0]
        print(chr)
