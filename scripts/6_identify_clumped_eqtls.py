import pandas as pd

clumps = pd.read_csv("output/clumps/RAD21_clumped_r0.8.clumped", delim_whitespace=True)

final_eqtls = clumps[["SNP"]]
tf = "RAD21"

final_eqtls.to_csv(f"{tf}_independent_eQTLs.csv", index = False)
