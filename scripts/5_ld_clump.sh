# Load: Modules and Environments
module load igmm/apps/plink/1.90b4

# Command
plink --bfile output/bed_files/GenOMICC_ISARIC_info_score_0.6_chr6 --clump input/ESR1_ciseQTLs_hg38.assoc --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 0.2 --out output/clumps/ESR1_clumped_r0.2
