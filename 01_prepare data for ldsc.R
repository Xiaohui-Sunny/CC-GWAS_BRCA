
####preparation of data
# heretability, genetic correlation and intercept were generated from LDSC (python)
#prepare data for ldsc
library(data.table)
library(R.utils)
setwd("/home/nfs/sunx3/bra_subtypes_ccgwas/data")
bc <- fread("icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
dim(bc) #9965310 53

# extract SNPs with maf >= 0.01 & info >= 0.8
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
bc.sub.pass.3$LumA_or <- exp(bc.sub.pass.3$Luminal_A_log_or_meta)
bc.sub.pass.3$LumA_z <- bc.sub.pass.3$Luminal_A_log_or_meta/bc.sub.pass.3$Luminal_A_se_meta
bc.sub.pass.3$LumA_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumA_z)))
bc.sub.pass.3$ngt <- 0

# Luminal-B
bc.sub.pass.3$LumB_or <- exp(bc.sub.pass.3$Luminal_B_log_or_meta)
bc.sub.pass.3$LumB_z <- bc.sub.pass.3$Luminal_B_log_or_meta/bc.sub.pass.3$Luminal_B_se_meta
bc.sub.pass.3$LumB_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_z)))

# Luminal-B-Her2-neg
bc.sub.pass.3$LumB_her2_neg_or <- exp(bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta)
bc.sub.pass.3$LumB_her2_neg_z <- bc.sub.pass.3$Luminal_B_HER2Neg_log_or_meta/bc.sub.pass.3$Luminal_B_HER2Neg_se_meta
bc.sub.pass.3$LumB_her2_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$LumB_her2_neg_z)))

# Her2-enrich
bc.sub.pass.3$Her2_enrich_or <- exp(bc.sub.pass.3$HER2_Enriched_log_or_meta)
bc.sub.pass.3$Her2_enrich_z <- bc.sub.pass.3$HER2_Enriched_log_or_meta/bc.sub.pass.3$HER2_Enriched_se_meta
bc.sub.pass.3$Her2_enrich_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Her2_enrich_z)))

# Triple-neg
bc.sub.pass.3$Trip_neg_or <- exp(bc.sub.pass.3$Triple_Neg_log_or_meta)
bc.sub.pass.3$Trip_neg_z <- bc.sub.pass.3$Triple_Neg_log_or_meta/bc.sub.pass.3$Triple_Neg_se_meta
bc.sub.pass.3$Trip_neg_p <- 2*(1-pnorm(abs(bc.sub.pass.3$Trip_neg_z)))

# match snp names with the ones in the ldscore files
# unzip the eru_w_ld data
for (i in 1:22){
  gunzip(paste0("/home/nfs/sunx3/software/ldsc/eur_w_ld_chr/",i,".l2.ldscore.gz"),remove=F)
}
#combine to a file
ldscore <- list()
for ( i in 1:22){
  ldscore [[i]] <- fread (paste0("/home/nfs/sunx3/software/ldsc/eur_w_ld_chr/",i,".l2.ldscore"))
}
ldscore <- rbindlist(ldscore, use.names=T)

bc.sub.mg <- merge(bc.sub.pass.3, ldscore, by.x=c('chr.Onco', 'Position.Onco'), by.y=c('CHR', 'BP'))
dim(bc.sub.mg) # 1129962    73
bc.sub.mg <- data.frame(bc.sub.mg)

# remove duplicates
var <- c('chr.Onco', 'Position.Onco')
sum(duplicated(bc.sub.mg[var])) # 252, remove them all
bc.sub.mg$chr_bp <- paste0(bc.sub.mg$chr.Onco, '_', bc.sub.mg$Position.Onco)
dup <- bc.sub.mg$chr_bp[duplicated(bc.sub.mg[var])]
bc.sub.f <- subset(bc.sub.mg, !chr_bp %in% dup)
dim(bc.sub.f) # 1129458      74

#write down the input file for ldsc
bc_LumA <- bc.sub.f[,c("SNP","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumA_or","Luminal_A_se_meta","LumA_p","r2.Onco","ngt","EAFcontrols.Onco")]
colnames(bc_LumA) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf')
write.table(bc_LumA, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaA_ldsc_input.txt', row.names=F, quote=F)

bc_LumB <- bc.sub.f[,c("SNP","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumB_or","Luminal_B_se_meta","LumB_p","r2.Onco","ngt","EAFcontrols.Onco")]
colnames(bc_LumB) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf')
write.table(bc_LumB, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_ldsc_input.txt', row.names=F, quote=F)

bc_LumB_her2 <- bc.sub.f[,c("SNP","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","LumB_her2_neg_or","Luminal_B_HER2Neg_se_meta","LumB_her2_neg_p","r2.Onco","ngt","EAFcontrols.Onco")]
colnames(bc_LumB_her2) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf')
write.table(bc_LumB_her2, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_LuminaB_Her2_ldsc_input.txt', row.names=F, quote=F)

bc_Her2_enrich <- bc.sub.f[,c("SNP","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","Her2_enrich_or","HER2_Enriched_se_meta","Her2_enrich_p","r2.Onco","ngt","EAFcontrols.Onco")]
colnames(bc_Her2_enrich) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf')
write.table(bc_Her2_enrich, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_Her2_enrich_ldsc_input.txt', row.names=F, quote=F)

bc_Trip <- bc.sub.f[,c("SNP","chr.iCOGs","Position.iCOGs","Effect.iCOGs","Baseline.iCOGs","Trip_neg_or","Triple_Neg_se_meta","Trip_neg_p","r2.Onco","ngt","EAFcontrols.Onco")]
colnames(bc_Trip) <- c('snpid','hg18chr','bp','a1','a2','or','se','pval','info','ngt','CEUaf')
write.table(bc_Trip, '/home/nfs/sunx3/bra_subtypes_ccgwas/data/BCAC_Trip_ldsc_input.txt', row.names=F, quote=F)


###############################################################################
##### correlation matrix plot
library(corrplot)
library(RColorBrewer)
res <- read.csv("G:/My work/Project/BRA subtypes_ccGWAS/results/ldsc/figure.csv",head=T)
rownames(res) <- res$ï..
res2 <- as.matrix(res[,-1])
corrplot(res2, method="circle",col.lim = c(0.5, 1), is.corr=FALSE, col = COL1('Blues',10),addCoef.col = 'grey55',addgrid.col = 'white',type = 'lower',
         tl.col = 'black',tl.srt = 25)
corrplot(res2, method="circle",col.lim = c(0.4, 1), is.corr=FALSE, col = COL2('RdYlBu'),addgrid.col = 'grey80',
         tl.col = 'black',tl.srt = 35,addCoef.col = 'black')
.
