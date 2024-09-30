import pandas as pd

df = pd.read_csv("chr16.snpstats", delimiter="\t", skiprows=9)
rsidsfailQC = list(df[(df["info"] < 0.9) | (df["minor_allele_frequency"] < 0.01)].rsid.values)

## Identify SNPs to include after imposing these filter
df2include = df[~df['rsid'].isin(rsidsfailQC)]

# Identify any remaining multiallelc SNPs
duplicated_mask = df2include["rsid"].duplicated(keep=False)  # keep=False marks all duplicates as True
# Filter the DataFrame to show only the duplicated rows
duplicated_rows = df2include[duplicated_mask]

# Add multiallelic SNPs remaining to rsids2exclude
multiallelic_snps = list(duplicated_rows.rsid.unique())
rsids2exclude = rsidsfailQC + multiallelic_snps

rsids = " ".join(rsids2exclude)

# Write rsids
with open("output/chr16_rsids2exclude_info_score0.9.txt", "w") as text_file:
    text_file.write(rsids)
