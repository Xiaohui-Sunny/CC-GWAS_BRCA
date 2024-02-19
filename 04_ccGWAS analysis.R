# introduction
#https://github.com/wouterpeyrot/CCGWAS

# bug fixed
#https://github.com/wouterpeyrot/CCGWAS/issues/2


#preparation of packages
library(data.table)
library(R.utils)#unzip file
library(devtools)
#install_github("wouterpeyrot/CCGWAS")
library(CCGWAS)

# ccGWAS analyze 
#Definition of the variables in the code
#outcome_file : name of the file where the outcome should be saved 
#A_name/B_name: disorder names
#sumstas_fileA1A0/sumstas_fileB1B0: summary data(!!zip-file) for two traits. 
#column names are SNP, CHR, BP, EA(effect allele),NEA, FRQ(frequency of the EA in controls), OR SE,P, Neff(effective sample size=4/{(1/N_case)+(1/N_control)}),
#K_A1A0/K_B1B0: the most likely lifetime disorder prevalences of disorder A and disorder B in the population
#K_A1A0_high/K_B1B0_high: upper bound of disorder prevalences. 
#K_A1A0_low/K_B1B0_low: lower bound of disorder prevalences.
#h2l_A1A0/h2l_B1B0: SNP-based heritability on the liability scale, as estimated with e.g. stratified LD Score regression
#rg_A1A0_B1B0: genetic correlation between disorders A and B, as estimated with e.g. cross-trait LD score regression
#intercept_A1A0_B1B0: intercept from cross-trait LD score regression
#m: approximation of number of independent effective loci
#N_A1/N_B1: total number of cases in the respective input case-control GWAS
#N_A0/N_B0: total number of controls in the respective input case-control GWAS
#N_overlap_A0B0: confirmed number of overlapping controls between the A1A0 and B1B0 GWAS samples.
#subtype_data: set to FALSE when comparing two different disorders, set TRUE when comparing subtypes of a disorder
#sumstats_fileA1B1: set to NA when applying CC-GWAS based on case-control GWAS results only.
#N_A1_inA1B1/N_B1_inA1B1: set to NA when applying CC-GWAS based on case-control GWAS results only
#intercept_A1A0_A1B1/intercept_B1B0_A1B1: set to NA when applying CC-GWAS based on case-control GWAS results only.

setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas")
# luminal-B vs luminal-A 
#CCGWAS( outcome_file = "lumB_lumA.out" , A_name = "LMB" , B_name = "LMA" , 
 #       sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_ccgwas_input_changeneff.txt.gz" ,
 #       sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input.txt.gz" ,
 #       K_A1A0 = 0.014 , K_A1A0_high = 0.021 , K_A1A0_low = 0.007 , 
 #       K_B1B0 = 0.066 , K_B1B0_high = 0.094 , K_B1B0_low = 0.047 ,  
 #       h2l_A1A0 = 1.4388 , h2l_B1B0 = 0.2019 , rg_A1A0_B1B0 = 0.7815 , intercept_A1A0_B1B0 = 0.1141 , m = 5000 ,
 #       N_A1 = 6659 , N_B1 = 45436 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )
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

# N_A1 = 6659 , N_B1 = 45436 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477
# N_A1/B1 are Neff (e.g. check data_extr.R, median of Neff for SNPs)
# N_A0/B0 is number of controls, which is fixed at 91477 (check the original BC subtype GWAS paper)
# N_overlap_A0B0 is the number of controls since we are doing case-case comparisons for subtypes of cancer

# Error messages: LumA
# ...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.6304 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Error in CCGWAS(outcome_file = "lumA_lumB.out", A_name = "LMA", B_name = "LMB",  :
#  ...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.6304 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Solution: check https://github.com/wouterpeyrot/CCGWAS/issues/2
# Change Neff to 1.6304^2*Neff_old
# BCAC_LuminaA_meta_results_050321.txt_v2.gz

# Error messages: LumB
#...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.9317 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
#Error in CCGWAS(outcome_file = "lumA_lumB.out", A_name = "LMA", B_name = "LMB",  :
#  ...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.9317 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Change Neff to 1.9317^2*Neff_old
# BCAC_LuminaB_meta_results_050321.txt_v2.gz


# luminal-B vs luminal-A_(changed Neff)
CCGWAS( outcome_file = "lumB_lumA.out" , A_name = "LMB" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.014 , K_A1A0_high = 0.021 , K_A1A0_low = 0.007 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094 , K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.531 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.7765 , intercept_A1A0_B1B0 = 0.1148 , m = 5000 ,  
        N_A1 = 12414 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )

# lumB-her2-neg vs lumA , rg=0.8709 > 0.8
# skipped

# her2 enrich vs lumA 
#CCGWAS( outcome_file = "her2_enrich_lumA.out" , A_name = "HER" , B_name = "LMA" , 
#       sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input.txt.gz" ,
#      sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input_changeneff.txt.gz" ,
#        K_A1A0 = 0.005 , K_A1A0_high = 0.008 , K_A1A0_low = 0.003 , 
#       K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
#      h2l_A1A0 = 0.0316 , h2l_B1B0 = 0.0405 , rg_A1A0_B1B0 = 0.775 , intercept_A1A0_B1B0 = 0.13 , m = 5000 ,  
#      N_A1 = 3026 , N_B1 = 45436 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )
# For Her2 enriched BC, ER- and PR-, HER2+
# ER- and PR-, HER2+ BC is about 4% of total BC (no range is available, just assume +/- 50%)
# K_B1B0 = 0.1283*0.04=0.005132
# K_B1B0_high = 0.1283*0.04*1.5 = 0.008
# K_B1B0_low = 0.1283*0.04*0.5 = 0.003

# Error in CCGWAS(outcome_file = "lumA_her2_enrich.out", A_name = "LMA",  :
#   ...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.9683 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Change Neff to 1.9683^2*Neff_old

# Her2-enriched vs luminal-A_(changed Neff)
CCGWAS( outcome_file = "her2_enrich_lumA.out" , A_name = "HER" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.005 , K_A1A0_high = 0.008 , K_A1A0_low = 0.003 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.9839 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.6752 , intercept_A1A0_B1B0 = 0.1213 , m = 5000 ,  
        N_A1 = 5859 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data=TRUE )

# triple negative vs lumA 
#CCGWAS( outcome_file = "trip_lumA.out" , A_name = "TRP" , B_name = "LMA" , 
#        sumstats_fileA1A0 = "/home/nfs/sunx3/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input.txt.gz" ,
#        sumstats_fileB1B0 = "/home/nfs/sunx3/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input_change neff.txt.gz" ,
#        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 , 
#        K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
#        h2l_A1A0 = 0.0322 , h2l_B1B0 = 0.0404 , rg_A1A0_B1B0 = 0.5483 , intercept_A1A0_B1B0 = 0.2198 , m = 5000 ,  
#        N_A1 = 9067 , N_B1 = 45436 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data = TRUE )
# For triple negative BC, ER- and PR-, HER2-
# ER- and PR-, HER2- BC is about 12% of total BC (no range is available, just assume +/- 50%)
# K_B1B0 = 0.1283*0.12=0.015
# K_B1B0_high = 0.1283*0.12*1.5 = 0.023
# K_B1B0_low = 0.1283*0.12*0.5 = 0.008

#Error in CCGWAS(outcome_file = "lumA_triple.out", A_name = "LMA", B_name = "TRP",  :
#  ...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.9088 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Change Neff to 1.9088^2*Neff_old

CCGWAS( outcome_file = "trip_lumA.out" , A_name = "TRP" , B_name = "LMA" , 
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaA_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 , 
        K_B1B0 = 0.066 , K_B1B0_high = 0.094, K_B1B0_low = 0.047 ,  
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.1422 , rg_A1A0_B1B0 = 0.5305 , intercept_A1A0_B1B0 = 0.2123 , m = 5000 ,  
        N_A1 = 16499 , N_B1 = 60715 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477 , subtype_data = TRUE )




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

# luminalB-Her2-neg vs Her2 enriched 0.8185>0.8


# triple negative vs her2 enriched 
CCGWAS( outcome_file = "trip_her2.out" , A_name = "TRP" , B_name = "HER" ,
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Her2_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
        K_B1B0 = 0.005 , K_B1B0_high = 0.008 , K_B1B0_low = 0.003 ,
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.9839 , rg_A1A0_B1B0 = 0.7994 , intercept_A1A0_B1B0 = -0.0871, m = 5000 ,
        N_A1 = 16499, N_B1 = 5859 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )


# triple negative vs luminalB-Her2-neg
#CCGWAS( outcome_file = "trip_lumB_her2_neg.out" , A_name = "TRP" , B_name = "LBH" ,
#        sumstats_fileA1A0 = "/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_Triple_ccgwas_input_change neff.txt.gz" ,
#        sumstats_fileB1B0 = "/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_her2_ccgwas_input.txt.gz" ,
#        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
#        K_B1B0 = 0.028 , K_B1B0_high = 0.047 , K_B1B0_low = 0.009 ,
#        h2l_A1A0 = 0.0338 , h2l_B1B0 = 0.0424 , rg_A1A0_B1B0 = 0.529 , intercept_A1A0_B1B0 = 0.0983, m = 5000 ,
#        N_A1 = 9067 , N_B1 = 9336 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )

# For luminalB Her2 negative BC, ER+ and/or PR+, HER2-, grade 3
# ER+ and/or PR+, HER2- BC is about 73%
# Assume grade 1 and 2 accounts for 30% of BC (max 50%, min 10%)
# K_B1B0 = 0.1283*0.73*0.3 = 0.028
# K_B1B0_high = 0.1283*0.73*0.5 = 0.047
# K_B1B0_low = 0.1283*0.73*0.1 = 0.009

#Error in CCGWAS(outcome_file = "lumB_her2_neg_trip.out", A_name = "LBH",  :
#...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = 1.9062 > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.
# Change Neff to 1.9062^2*Neff_old

CCGWAS( outcome_file = "trip_lumB_her2_neg.out" , A_name = "TRP" , B_name = "LBH" ,
        sumstats_fileA1A0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_Triple_ccgwas_input0428_change neff.txt.gz" ,
        sumstats_fileB1B0 = "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ccgwas/BCAC_LuminaB_her2_ccgwas_input0428_change neff.txt.gz" ,
        K_A1A0 = 0.015 , K_A1A0_high = 0.023 , K_A1A0_low = 0.008 ,
        K_B1B0 = 0.028 , K_B1B0_high = 0.047 , K_B1B0_low = 0.009 ,
        h2l_A1A0 = 0.2984 , h2l_B1B0 = 0.3466 , rg_A1A0_B1B0 = 0.5171 , intercept_A1A0_B1B0 = 0.0944, m = 5000 ,
        N_A1 = 16499 , N_B1 = 16943 , N_A0 = 91477 , N_B0 = 91477 , N_overlap_A0B0 = 91477, subtype_data = TRUE )
