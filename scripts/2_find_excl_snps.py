import pandas as pd

df = pd.read_csv("chr16.snpstats", delimiter="\t", skiprows=10)

# Remove any multiallelic SNPs
biallelic_df = df.drop_duplicates(subset = "rsid", keep = False)

# Only keep SNPs passing QC thresholds
rsidsPassQC = list(biallelic_df[(biallelic_df["info"] >= 0.9) | (biallelic_df["minor_allele_frequency"] >= 0.01)].rsid.values) 

rsids = " ".join(rsidsPassQC)

# Write rsids
with open("output/chr16_rsids2include_info_score0.9.txt", "w") as text_file:
    text_file.write(rsids)
