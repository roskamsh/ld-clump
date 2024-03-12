import numpy as np
import pandas as pd

df = pd.read_csv("output/chr6_imputation_stats.snpstats", delimiter="\t", skiprows=9)
df = df[df["info"] < 0.6]

snps = df[["alternate_ids","rsid","chromosome","position","alleleA","alleleB"]]
snps.columns = ["SNPID","rsid","chromosome","position","alleleA","alleleB"] 

rsids = " ".join(snps.rsid.values)
snps.to_csv("output/chr6_variants2exclude_info_score0.6.txt", index = False, sep="\t")

# Write rsids
with open("output/chr6_rsids2exclude_info_score0.6.txt", "w") as text_file:
    text_file.write(rsids)

