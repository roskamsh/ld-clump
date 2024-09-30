import os
import pandas as pd

dir = "test/baal-nf-out"
concordant_value = bool("True")

# Read in all CSV files, and filter by Concordant value
files = os.listdir(dir)

# Initialize an empty list to store DataFrames
dataframes = []

# Loop through the CSV files and read them into DataFrames
for file in files:
    file_path = os.path.join(dir, file)
    df = pd.read_csv(file_path)
    if concordant_value != "":
        df = df[df.Concordant == concordant_value] 
    dataframes.append(df)

results = pd.concat(dataframes, ignore_index=True)

# Define TF list from results
tfs = pd.DataFrame(list(results.tf.unique()))
filt = results[["CHROM","tf"]].copy()
filt = filt.drop_duplicates()

results.to_csv("all_bQTLs.csv", index = False)
tfs.to_csv("TF_list.csv", index = False, header=False)
filt.to_csv("bQTL_CHR_list.csv", index = False)