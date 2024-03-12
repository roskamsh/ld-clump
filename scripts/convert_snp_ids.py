import numpy as np
import pandas as pd

genotyped_snps = pd.read_csv("input/genotyped_snps.csv")
GSA_conversion = pd.read_csv("input/GSA-24v3-0_A1_b151_rsids.txt",  delimiter="\t")
GSA_conversion.columns = ['Illumina_ID','RSID']
genotyped_snps.columns = ['Illumina_ID',"BGEN_ID"]
# downloaded from https://emea.support.illumina.com/downloads/infinium-global-screening-array-v3-0-support-files.html

merged = pd.merge(genotyped_snps, GSA_conversion, on = 'Illumina_ID', how = 'left')

merged.to_csv("input/genotyped_snps_converted.csv", index = False)