import pandas as pd
import os

results_dir = "output-ukbb/bqtl_ld_results"
files = [f for f in os.listdir(results_dir) if f.endswith('_ld_results.ld')]

compiled_dataframes = []

# Iterate through the list of files
for file in files:
    # Read the file into a DataFrame, skip if empty
    try:
        df = pd.read_csv(os.path.join(results_dir,file), delim_whitespace=True, header=0)
        if df.empty:
            print(f"Skipping empty file: {file}. No variants in LD with one another.")
            continue
        # Add a new column with the file name (without '_ld_results.ld')
        df['TF'] = file.replace('_ld_results.ld', '')
        df = df[df['R2'] >= 0.8]
        compiled_dataframes.append(df)
    except pd.errors.EmptyDataError:
        print(f"File is empty or unreadable: {file}")
        continue

if compiled_dataframes:
    compiled_df = pd.concat(compiled_dataframes, ignore_index=True)
    print("Number of bQTL pairs in LD across all TFs inspected:", compiled_df.shape[0])
    print("Saving bQTLs in LD to file: bQTLs_in_LD.csv")
    compiled_df.to_csv("bQTLs_in_LD.csv", index = False)
else:
    compiled_df = pd.DataFrame()
    print("No bQTLs found to be in LD with one another. Exiting...")
