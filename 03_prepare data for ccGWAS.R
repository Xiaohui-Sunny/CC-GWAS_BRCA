# Prepare data for ccGWAS analyses
library(data.table)
library(R.utils)

setwd("/home/nfs/sunx3/bra_subtypes_ccgwas/data")
bc <- fread("icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
bc <- fread("G:/My work/BRA subtypes_ccGWAS/results/try.txt",head=T)
dim(bc) #9965310 53

# extract SNPs with maf >= 0.01 & info >= 0.8
bc.sub.pass <- subset(bc, r2.Onco >= 0.8 & EAFcontrols.Onco >= 0.01 & EAFcontrols.Onco <=0.99)
dim(bc.sub.pass) # 8709696      53

# remove INDELs
bc.sub.pass.2 <- subset(bc.sub.pass, Effect.Meta %in% c('A', 'C', 'G', 'T') & Baseline.Meta %in% c('A', 'C', 'G', 'T'))
dim(bc.sub.pass.2) # 7745309      53
table(bc.sub.pass.2$Effect.Meta, bc.sub.pass.2$Baseline.Meta)

# remove A-T/T-A, C-G/G-C
bc.sub.pass.3 <- subset(bc.sub.pass.2, !(Effect.Meta =='A' & Baseline.Meta =='T' | Effect.Meta =='T' & Baseline.Meta =='A' |
                                           Effect.Meta =='C' & Baseline.Meta =='G' | Effect.Meta =='G' & Baseline.Meta =='C'))
rm(bc)
rm(bc.sub.pass)
rm(bc.sub.pass.2)

# output the summary data for breast subtypes
## Luminal-A
bc.sub.pass.3$LumA_or <- exp(bc.sub.pass.3$Luminal_A_log_or_meta)
bc.sub.pass.3$LumA_z <- bc.sub.pass.3$Luminal_A_log_or_meta/bc.sub.pass.3$Luminal_A_se_meta
bc.sub.pass.3$LumA_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumA_z)))
LumA_Neff <- 91477/(2*bc.sub.pass.3$Luminal_A_se_meta^2*bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*91477-1)
summary(LumA_Neff) # Median: N_eff= 45436
bc.sub.pass.3$LumA_Neff <- LumA_Neff
# remove those with Neff out of +/- 2sd of Neff
bc.LumA <- bc.sub.pass.3[bc.sub.pass.3$LumA_Neff>=median(bc.sub.pass.3$LumA_Neff)-2*sd(bc.sub.pass.3$LumA_Neff) & bc.sub.pass.3$LumA_Neff<=median(bc.sub.pass.3$LumA_Neff)+2*sd(bc.sub.pass.3$LumA_Neff), c('SNP.Onco', 'chr.Onco', 'Position.Onco', 'Effect.Onco', 'Baseline.Onco', 'EAFcontrols.Onco', 'LumA_or_meta', 'Luminal_A_se_meta', 'LumA_p', 'LumA_Neff')]
colnames(bc.LumA) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.LumA, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaA_ccgwas_input.txt', row.names=F, quote=F)
rm(bc.LumA)

## Luminal-B
bc.sub.pass.3$LumB_or <- exp(bc.sub.pass.3$Luminal_B_log_or_meta)
bc.sub.pass.3$LumB_z <- bc.sub.pass.3$Luminal_B_log_or_meta/bc.sub.pass.3$Luminal_B_se_meta
bc.sub.pass.3$LumB_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_z)))
LumB_Neff <- 91477/(2*bc.sub.pass.3$Luminal_B_se_meta^2*bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*91477-1)
summary(LumB_Neff) # Median: N_eff= 6659 
bc.sub.pass.3$LumB_Neff <- LumB_Neff
bc.LumB <- bc.sub.pass.3[bc.sub.pass.3$LumB_Neff>=median(bc.sub.pass.3$LumB_Neff)-2*sd(bc.sub.pass.3$LumB_Neff) & bc.sub.pass.3$LumB_Neff<=median(bc.sub.pass.3$LumB_Neff)+2*sd(bc.sub.pass.3$LumB_Neff), c('SNP.Onco', 'chr.Onco', 'Position.Onco', 'Effect.Onco', 'Baseline.Onco', 'EAFcontrols.Onco', 'LumB_or_meta', 'Luminal_B_se_meta', 'LumB_p', 'LumB_Neff')]
colnames(bc.LumB) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.LumB, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_ccgwas_input.txt', row.names=F, quote=F)
rm(bc.LumB)

## Luminal-B-Her2-neg
bc.sub.pass.3$LumB_her2_neg_or <- exp(bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta)
bc.sub.pass.3$LumB_her2_neg_z <- bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta/bc.sub.pass.3$Luminal_B_HER2Neg_se_meta
bc.sub.pass.3$LumB_her2_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_her2_neg_z)))
LumB_her2_neg_Neff <- 91477/(2*bc.sub.pass.3$Luminal_B_HER2Neg_se_meta^2*bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*91477-1)
summary(LumB_her2_neg_Neff) # Median: N_eff= 9336.2
bc.sub.pass.3$LumB_her2_neg_Neff <- LumB_her2_neg_Neff
bc.LumB.her2.neg <- bc.sub.pass.3[bc.sub.pass.3$LumB_her2_neg_Neff>=median(bc.sub.pass.3$LumB_her2_neg_Neff)-2*sd(bc.sub.pass.3$LumB_her2_neg_Neff) & bc.sub.pass.3$LumB_her2_neg_Neff<=median(bc.sub.pass.3$LumB_her2_neg_Neff)+2*sd(bc.sub.pass.3$LumB_her2_neg_Neff), c('SNP.Onco', 'chr.Onco', 'Position.Onco', 'Effect.Onco', 'Baseline.Onco', 'EAFcontrols.Onco', 'LumB_her2_neg_or_meta', 'Luminal_B_HER2Neg_se_meta', 'LumB_her2_neg_p', 'LumB_her2_neg_Neff')]
colnames(bc.LumB.her2.neg) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.LumB.her2.neg, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_her2_ccgwas_input.txt', row.names=F, quote=F)
rm(bc.LumB.her2.neg)

## Her2-enrich
bc.sub.pass.3$Her2_enrich_or <- exp(bc.sub.pass.3$HER2_Enriched_log_or_meta)
bc.sub.pass.3$Her2_enrich_z <- bc.sub.pass.3$HER2_Enriched_log_or_meta/bc.sub.pass.3$HER2_Enriched_se_meta
bc.sub.pass.3$Her2_enrich_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Her2_enrich_z)))
Her2_enrich_Neff <- 91477/(2*bc.sub.pass.3$HER2_Enriched_se_meta^2*bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*91477-1)
summary(Her2_enrich_Neff) # Median: N_eff= 3026.5bc.sub.pass.3$Her2_enrich_Neff <- Her2_enrich_Neff
bc.her2 <- bc.sub.pass.3[bc.sub.pass.3$Her2_enrich_Neff>=median(bc.sub.pass.3$Her2_enrich_Neff)-2*sd(bc.sub.pass.3$Her2_enrich_Neff) & bc.sub.pass.3$Her2_enrich_Neff<=median(bc.sub.pass.3$Her2_enrich_Neff)+2*sd(bc.sub.pass.3$Her2_enrich_Neff), c('SNP.Onco', 'chr.Onco', 'Position.Onco', 'Effect.Onco', 'Baseline.Onco', 'EAFcontrols.Onco', 'Her2_enrich_or_meta', 'HER2_Enriched_se_meta', 'Her2_enrich_p', 'Her2_enrich_Neff')]
colnames(bc.her2) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.her2, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_Her2_ccgwas_input.txt', row.names=F, quote=F)
rm(bc.her2)

## Triple-neg
bc.sub.pass.3$Trip_neg_or <- exp(bc.sub.pass.3$Triple_Neg_log_or_meta)
bc.sub.pass.3$Trip_neg_z <- bc.sub.pass.3$Triple_Neg_log_or_meta/bc.sub.pass.3$Triple_Neg_se_meta
bc.sub.pass.3$Trip_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Trip_neg_z)))
Trip_neg_Neff <- 91477/(2*bc.sub.pass.3$Triple_Neg_se_meta^2*bc.sub.pass.3$EAFcontrols.Onco*(1-bc.sub.pass.3$EAFcontrols.Onco)*91477-1)
summary(Trip_neg_Neff) # Median: N_eff= 9067.7
bc.sub.pass.3$Trip_neg_Neff <- Trip_neg_Neff
bc.trip <- bc.sub.pass.3[bc.sub.pass.3$Trip_neg_Neff>=median(bc.sub.pass.3$Trip_neg_Neff)-2*sd(bc.sub.pass.3$Trip_neg_Neff) & bc.sub.pass.3$Trip_neg_Neff<=median(bc.sub.pass.3$Trip_neg_Neff)+2*sd(bc.sub.pass.3$Trip_neg_Neff), c('SNP.Onco', 'chr.Onco', 'Position.Onco', 'Effect.Onco', 'Baseline.Onco', 'EAFcontrols.Onco', 'Trip_neg_or', 'Triple_Neg_se_meta', 'Trip_neg_p', 'Trip_neg_Neff')]
colnames(bc.trip) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.trip, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_Triple_ccgwas_input.txt', row.names=F, quote=F)
rm(bc.trip)

bc.LumB.her2.neg <- fread(gunzip("/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_her2_ccgwas_input.txt.gz",remove=F))
bc.LumB.her2.neg$Neff2 <- 1.9062^2*bc.LumB.her2.neg$Neff
bc.LumB.her2.neg <- bc.LumB.her2.neg[,c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff2')]
colnames(bc.LumB.her2.neg) <- c('SNP','CHR','BP','EA','NEA','FRQ','OR','SE','P','Neff')
write.table(bc.LumB.her2.neg, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_her2_ccgwas_input_change neff.txt", row.names=F, quote=F)
gzip("/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_her2_ccgwas_input_change neff.txt")
rm(bc.LumB.her2.neg)
