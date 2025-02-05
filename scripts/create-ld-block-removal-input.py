import pandas as pd
import os

final_bqtls = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_bQTLs.csv")
final_transactors = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_transactors.csv")
input_bqtls = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/f4/f4d35fd5d5b5556307ce65e87b52ce/all_bQTLs.csv")

transactor_list = []

for filename in os.listdir("."):
    if "transactor_QTLs.csv" in filename: 
                df = pd.read_csv(filename)  
                transactor_list.append(df) 

input_transactors = pd.concat(transactor_list, ignore_index=True)

annot_bqtls = final_bqtls.merge(input_bqtls, 
                                left_on = ['SNP','TF'], 
                                right_on = ['ID','tf'], 
                                how = 'inner')

annot_bqtls["CHR_num"] = annot_bqtls['CHROM'].str.replace(r'^CHR', '', case=False, regex=True)
annot_bqtls = annot_bqtls[["SNP","CHR_num","POS"]]
annot_bqtls.columns = ["RSID","CHR","POS"]

annot_transactors = final_transactors.merge(input_transactors,
                                            left_on = ['SNP','TF'],
                                            right_on = ['SNP','GeneSymbol'],
                                            how = 'inner')
annot_transactors = annot_transactors[["SNP","SNPChr","SNPPos"]]
annot_transactors.columns = ["RSID","CHR","POS"] 

final = pd.concat([annot_bqtls, annot_transactors], ignore_index=True)
final.to_csv("ld-block-removal-input.csv", index = False)