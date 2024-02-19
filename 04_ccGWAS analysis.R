
library(data.table)
library(R.utils)#unzip file
library(devtools)
library(CCGWAS)

# ccGWAS analyze 

setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas")
# luminal-B vs luminal-A 
CCGWAS( outcome_file = "lumB_lumA.out" , A_name = "LMB" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.014 , K_A1A0_high = 0.021 , K_A1A0_low = 0.007 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094 , K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.531 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.7765 , intercept_A1A0_B1B0 = 0.1148 , m = 5000 ,  
        N_A1 = 12414 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )
# Definition: Luminal B-like, ER+ and/or PR+, HER2+
# Definition: Luminal A-like, ER+ and/or PR+, HER2-, grade 1 and 2
# K_A1A0/B1B0: the most likely lifetime disorder prevalences of disorder A and disorder B in the population, following the disorder definition used in the respective case-control GWAS
# K_A1A0_low/K_B1B0_high/low: higher/lower bound of disorder prevalences.
# h2l_A1A0/h2l_B1B0: SNP-based heritability on the liability scale, as estimated with e.g. stratified LD Score regression
# intercept_A1A0_B1B0: intercept from cross-trait LD score regression
# m: When estimates of genome-wide polygenicity are not available, our recommendation is to specify m=1,000 for traits that are expected to have relatively sparse architectures (e.g. autoimmune diseases), m=10,000 for traits that are expected to have highly polygenic architectures (e.g. psychiatric disorders), and m=5,000 for traits with no clear expectation.

# The lifetime prevalence of breast cancer in women in the US is about 1/8 (12.83% from SEER 2014-2016)
# For Luminal B-like BC, ER+ and/or PR+, HER2+
# ER+ and/or PR+, HER2+ BC is about 11% (no range is available, just assume +/- 50%)
# K_A1A0 = 0.1283*0.11 = 0.014
# K_A1A0_high = 0.1283*0.11*1.5 = 0.021
# K_A1A0_low = 0.1283*0.11*0.5 = 0.007

# For luminal-A BC, ER+ and/or PR+, HER2- BC is about 73%
# Assume grade 1 and 2 accounts for 70% of BC (max 100%, min 50%)
# K_B1B0 = 0.1283*0.73*0.7 = 0.066
# K_B1B0_high = 0.1283*0.73*1 = 0.094
# K_B1B0_low = 0.1283*0.73*0.5 = 0.047

# h2l_A1A0 = 0.0844; h2l_B1B0 = 0.0859; rg_A1A0_B1B0 = 0.7845; intercept_A1A0_B1B0 = 0.1245; m= 5000
# h2l_A1A0/B1B0, rg_A1A0_B1B0 and intercepts are estimated by cross-trait ldsc (e.g. check LumA_LumB_liability.log) 


# her2 enrich vs lumA 
CCGWAS( outcome_file = "her2_enrich_lumA.out" , A_name = "HER" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.005 , K_A1A0_high = 0.008 , K_A1A0_low = 0.003 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.9839 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.6752 , intercept_A1A0_B1B0 = 0.1213 , m = 5000 ,  
        N_A1 = 5859 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )

# For Her2 enriched BC, ER- and PR-, HER2+
# ER- and PR-, HER2+ BC is about 4% of total BC (no range is available, just assume +/- 50%)
# K_B1B0 = 0.1283*0.04=0.005132
# K_B1B0_high = 0.1283*0.04*1.5 = 0.008
# K_B1B0_low = 0.1283*0.04*0.5 = 0.003


# triple negative vs lumA 
CCGWAS( outcome_file = "trip_lumA.out" , A_name = "TRP" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.5305 , intercept_A1A0_B1B0 = 0.2123 , m = 5000 ,  
        N_A1 = 16499 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data = TRUE )

# For triple negative BC, ER- and PR-, HER2-
# ER- and PR-, HER2- BC is about 12% of total BC (no range is available, just assume +/- 50%)
# K_B1B0 = 0.1283*0.12=0.015
# K_B1B0_high = 0.1283*0.12*1.5 = 0.023
# K_B1B0_low = 0.1283*0.12*0.5 = 0.008


# her2 enriched vs lumB
CCGWAS( outcome_file = "her2_enrich_lumB.out" , A_name = "HER" , B_name = "LMB" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.005 , K_A1A0_high = 0.008 , K_A1A0_low = 0.003 , 
        K_B1B0 = 0.014 , K_B1B0_high = 0.021 , K_B1B0_low = 0.007 ,  
        h2l_A1A0 = 0.9839 , h2l_B1B0 = 0.531 , rg_A1A0_B1B0 = 0.6006 , intercept_A1A0_B1B0 =0.0559 , m = 5000 ,  
        N_A1 = 5859 , N_B1 = 12414 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data = TRUE )

#  triple neg vs luminal-B 
CCGWAS( outcome_file = "trip_lumB.out" , A_name = "TRP" , B_name = "LMB" ,
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
        K_B1B0 = 0.014 , K_B1B0_high = 0.021 , K_B1B0_low = 0.007 ,
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.531, rg_A1A0_B1B0 = 0.7101 , intercept_A1A0_B1B0 = 0.0734, m = 5000 ,
        N_A1 = 16499 , N_B1 = 12414 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )

# triple negative vs her2 enriched 
CCGWAS( outcome_file = "trip_her2.out" , A_name = "TRP" , B_name = "HER" ,
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
        K_B1B0 = 0.005 , K_B1B0_high = 0.008 , K_B1B0_low = 0.003 ,
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.9839 , rg_A1A0_B1B0 = 0.7994 , intercept_A1A0_B1B0 = -0.0871, m = 5000 ,
        N_A1 = 16499, N_B1 = 5859 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )


# triple negative vs luminalB-Her2-neg
CCGWAS( outcome_file = "trip_lumB_her2_neg.out" , A_name = "TRP" , B_name = "LBH" ,
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_her2_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
        K_B1B0 = 0.028 , K_B1B0_high = 0.047 , K_B1B0_low = 0.009 ,
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.3466 , rg_A1A0_B1B0 = 0.5171 , intercept_A1A0_B1B0 = 0.0944, m = 5000 ,
        N_A1 = 16499 , N_B1 = 16943 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )

# For luminalB Her2 negative BC, ER+ and/or PR+, HER2-, grade 3
# ER+ and/or PR+, HER2- BC is about 73%
# Assume grade 1 and 2 accounts for 30% of BC (max 50%, min 10%)
# K_B1B0 = 0.1283*0.73*0.3 = 0.028
# K_B1B0_high = 0.1283*0.73*0.5 = 0.047
# K_B1B0_low = 0.1283*0.73*0.1 = 0.009

