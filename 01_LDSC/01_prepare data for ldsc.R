####preparation of data
# heretability, genetic correlation and intercept were generated from LDSC (python)

#prepare data for ldsc
library(data.table)
library(R.utils)
setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data")
bc <- read.table("icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
dim(bc) #9965310 53

# extract SNPs included in HampMap3 with maf >= 0.01 & info > 0.3 in the OncoArray data( following the NG)
hapmap3 <- fread("/home/nfs/sunx3/data/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map")  #1490422 The version is not hg19
hapmap3 <- fread("/home/nfs/sunx3/data/hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.map") #1437409 This one is transformed by using liftover
dim(hapmap3)
#hapmap3$chr <- substr(hapmap3$V1,start=4,stop=nchar(hapmap3$V1))
hapmap3$chr_bp <- paste0(hapmap3$V1,"_", hapmap3$V4)

bc.sub.pass <- subset(bc, r2.Onco >= 0.8 & EAFcontrols.Onco >= 0.01 & EAFcontrols.Onco <=0.99)
dim(bc.sub.pass) # 8709696      53    

# according to LDSC, remove INDELs
bc.sub.pass.2 <- subset(bc.sub.pass, Effect.Meta %in% c('A', 'C', 'G', 'T') & Baseline.Meta %in% c('A', 'C', 'G', 'T'))
dim(bc.sub.pass.2) # 7745309      53
table(bc.sub.pass.2$Effect.Meta, bc.sub.pass.2$Baseline.Meta)

# remove A-T/T-A, C-G/G-C
bc.sub.pass.3 <- subset(bc.sub.pass.2, !(Effect.Meta =='A' & Baseline.Meta =='T' | Effect.Meta =='T' & Baseline.Meta =='A' |
                                           Effect.Meta =='C' & Baseline.Meta =='G' | Effect.Meta =='G' & Baseline.Meta =='C'))
dim(bc.sub.pass.3)
# [1] 6532284      53
table(bc.sub.pass.3$Effect.Meta, bc.sub.pass.3$Baseline.Meta)

rm(bc.sub.pass.2)
rm(bc.sub.pass)
# output results by subtype
# Luminal-A
bc.sub.pass.3$ngt <- 0

bc.sub.pass.3$LumA_or <- exp(bc.sub.pass.3$Luminal_A_log_or_meta)
bc.sub.pass.3$LumA_z <- bc.sub.pass.3$Luminal_A_log_or_meta/bc.sub.pass.3$Luminal_A_se_meta
bc.sub.pass.3$LumA_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumA_z)))
bc.sub.pass.3$LumA_neff <- 1/(bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*bc.sub.pass.3$Luminal_A_se_meta^2)
summary(bc.sub.pass.3$LumA_neff) #median neff = 60715


# Luminal-B
bc.sub.pass.3$LumB_or <- exp(bc.sub.pass.3$Luminal_B_log_or_meta)
bc.sub.pass.3$LumB_z <- bc.sub.pass.3$Luminal_B_log_or_meta/bc.sub.pass.3$Luminal_B_se_meta
bc.sub.pass.3$LumB_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_z)))
bc.sub.pass.3$LumB_neff <- 1/(bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*bc.sub.pass.3$Luminal_B_se_meta^2)
summary(bc.sub.pass.3$LumB_neff) #median neff = 12414.5


# Luminal-B-Her2-neg
bc.sub.pass.3$LumB_her2_neg_or <- exp(bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta)
bc.sub.pass.3$LumB_her2_neg_z <- bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta/bc.sub.pass.3$Luminal_B_HER2Neg_se_meta
bc.sub.pass.3$LumB_her2_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_her2_neg_z)))
bc.sub.pass.3$LumB_her2_neff <- 1/(bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*bc.sub.pass.3$Luminal_B_HER2Neg_se_meta^2)
summary(bc.sub.pass.3$LumB_her2_neff) #median neff = 16943

# Her2-enrich
bc.sub.pass.3$Her2_enrich_or <- exp(bc.sub.pass.3$HER2_Enriched_log_or_meta)
bc.sub.pass.3$Her2_enrich_z <- bc.sub.pass.3$HER2_Enriched_log_or_meta/bc.sub.pass.3$HER2_Enriched_se_meta
bc.sub.pass.3$Her2_enrich_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Her2_enrich_z)))
bc.sub.pass.3$Her2_neff <- 1/(bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*bc.sub.pass.3$HER2_Enriched_se_meta^2)
summary(bc.sub.pass.3$Her2_neff) #median neff = 5859

# Triple-neg
bc.sub.pass.3$Trip_neg_or <- exp(bc.sub.pass.3$Triple_Neg_log_or_meta)
bc.sub.pass.3$Trip_neg_z <- bc.sub.pass.3$Triple_Neg_log_or_meta/bc.sub.pass.3$Triple_Neg_se_meta
bc.sub.pass.3$Trip_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Trip_neg_z)))
bc.sub.pass.3$Trip_neff <- 1/(bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*bc.sub.pass.3$Triple_Neg_se_meta^2)
summary(bc.sub.pass.3$Trip_neff) #median neff = 16499

# match snp names with the hapmap3, since the ldsc use hapmap3 to generate the ldsc input data
bc.sub.pass.3$chr_bp <- paste0(bc.sub.pass.3$chr.Onco, "_", bc.sub.pass.3$Position.Onco)


bc.sub.pass.m <- merge(bc.sub.pass.3, hapmap3, by.x="chr_bp", by.y="chr_bp")

# remove duplicates
var <- c('chr.Onco', 'Position.Onco')
sum(duplicated(bc.sub.pass.m[var])) # 0

#write down the input file for ldsc
bc.sub.pass.m <- as.data.frame(bc.sub.pass.m)

bc_LumA <- bc.sub.pass.m[,c("V2","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumA_or","Luminal_A_se_meta","LumA_p","r2.Onco","ngt","EAFcontrols.Onco","LumA_neff")]
colnames(bc_LumA) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf','N')
write.table(bc_LumA, '/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaA_ldsc_input_neff.txt', row.names=F, quote=F)

bc_LumB <- bc.sub.pass.m[,c("V2","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumB_or","Luminal_B_se_meta","LumB_p","r2.Onco","ngt","EAFcontrols.Onco","LumB_neff")]
colnames(bc_LumB) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf','N')
write.table(bc_LumB, '/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_ldsc_input_neff.txt', row.names=F, quote=F)

bc_LumB_her2 <- bc.sub.pass.m[,c("V2","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumB_her2_neg_or","Luminal_B_HER2Neg_se_meta","LumB_her2_neg_p","r2.Onco","ngt","EAFcontrols.Onco","LumB_her2_neff")]
colnames(bc_LumB_her2) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf','N')
write.table(bc_LumB_her2, '/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_Her2_ldsc_input_neff.txt', row.names=F, quote=F)

bc_Her2_enrich <- bc.sub.pass.m[,c("V2","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","Her2_enrich_or","HER2_Enriched_se_meta","Her2_enrich_p","r2.Onco","ngt","EAFcontrols.Onco","Her2_neff")]
colnames(bc_Her2_enrich) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf','N')
write.table(bc_Her2_enrich, '/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Her2_enrich_ldsc_input_neff.txt', row.names=F, quote=F)

bc_Trip <- bc.sub.pass.m[,c("V2","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","Trip_neg_or","Triple_Neg_se_meta","Trip_neg_p","r2.Onco","ngt","EAFcontrols.Onco","Trip_neff")]
colnames(bc_Trip) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf','N')
write.table(bc_Trip, '/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Trip_ldsc_input_neff.txt', row.names=F, quote=F)
