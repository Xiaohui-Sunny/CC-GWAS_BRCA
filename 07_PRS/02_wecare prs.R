library(data.table)
library(haven)
library(ggpubr)
library(ggplot2)

setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/prs_weare")

# prepare wecare covariate data
dat.cov <- read_sas("/home/nfs/sunx3/project/cbc_rt_interaction/data/covdata/covar_xs3_gwas_tx.sas7bdat")
dat.cov2 <- read_sas("/home/nfs/sunx3/project/cbc_rt_interaction/data/covdata/covar_xs2_overall_gwas.sas7bdat") # add nondanish in wecare 1 and wecare 2, respectivey
dat.cov2 <- dat.cov2[,c(25,46:55)]

dat.cov.m <- merge(dat.cov, dat.cov2, by.x="patient_id", by.y="patient_id")
dat.cov.nonda <-subset(dat.cov.m, dat.cov.m$sub_registry!="DE" & dat.cov.m$sub_registry!="de") 

# prepare geno data
dat.geno <- read.csv("WECARE.dosage.raw.csv",head=T)
rownames(dat.geno) <- dat.geno$IID
dat.geno <- select(dat.geno, -c("IID","PAT","MAT","SEX","PHENOTYPE")) 

#### Outcome: ER status######
### Association between each SNPs identified from CC-GWAS and ER risk 
dat <- merge(dat.cov.nonda, dat.geno, by.x="patient_id", by.y="row.names")#2887
## only included 0=non-Hispanic white :2483
table(dat$race)
dat <- subset(dat,dat$race==0)

## only included er  != NA
table(dat$er1_cat)
#negative positive  unknown 
#712     1451      320 
dat <- subset(dat, dat$er1_cat != "unknown")
dat$ER <- ifelse(dat$er1_cat == 'negative', 1, 0)


SNP_all_effect <- matrix(nrow = ncol(dat.geno), ncol = 3)
rownames(SNP_all_effect) <- colnames(dat.geno)
colnames(SNP_all_effect) <- c("beta", "se", "P")

for ( i in 1:27){
  glmfunc <- as.formula(paste0("ER ~ ", colnames(dat.geno)[i], "+ PC1_nonde + PC2_nonde + PC3_nonde + PC4_nonde + PC5_nonde + sub_dx_age + Phase + stage_fd"))
  re <- glm(glmfunc, data=dat, family = binomial)
  SNP_all_effect[i,] <- summary(re)$coefficients[2,c(1,2,4)]
  
}

write.csv(SNP_all_effect, "wecare_snp.csv", row.names = F)


##################################################################################
# Association between PRS and ER risk
prs.weights <- read.csv("PRS.ccgwas with proxy.csv", head=T)

prsSNPInfo <- as.data.frame(t(t(colnames(dat.geno))))
prsSNPInfo$EA.wecare <- sapply(strsplit(prsSNPInfo$V1,"_"),'[',2)
prsSNPInfo$chr <- sapply(strsplit(prsSNPInfo$V1,"[.]"),'[',1)
prsSNPInfo$chr <- lapply(prsSNPInfo$chr, function(x) as.numeric(gsub("[X,]", "", x)))
prsSNPInfo$BP37 <- sapply(strsplit(prsSNPInfo$V1,"[.]"),'[',2)
prsSNPInfo$chr_bp <- paste0(prsSNPInfo$chr, "_", prsSNPInfo$BP37)
colnames(prsSNPInfo) <- c("SNP","EA.wecare","chr","BP37","chr_bp")

prs.weights$chr_bp <- paste0(prs.weights$CHR,"_", prs.weights$BP..hg19.) #24


prs.weights <- merge(prs.weights, prsSNPInfo, by = "chr_bp")

# change beta based on the effect allele
for (i in 1:nrow(prs.weights)){
  if(prs.weights$EA.wecare[i] == prs.weights$EA[i]) {
    prs.weights$beta.change[i] <- prs.weights$CC.GWAS_beta[i]
  } else {
   prs.weights$beta.change[i] <- -prs.weights$CC.GWAS_beta[i]
}
}

# calculate weighted prs
dat.geno[] <- sapply(dat.geno, as.numeric)
dat.geno <- t(dat.geno)
prs.ccgwas <- t(dat.geno[prs.weights$SNP.y,]) %*% prs.weights$beta.change


dat <- merge(prs.ccgwas, dat.cov.nonda, by.x = "row.names", by.y = "patient_id")
colnames(dat)[2] <- c("prs.ccgwas")

dat$prs.ccgwas <- as.numeric(dat$prs.ccgwas)

dat <- subset(dat,dat$race==0) #2483
dat <- subset(dat, dat$er1_cat != "unknown")
dat$ER <- ifelse(dat$er1_cat == 'negative', 1, 0) #2163
dat$group <- ifelse(dat$er1_cat == "negative", "ER negative", "ER postive")

dat$prs.ccgwas.scale <- dat$prs.ccgwas/sd(dat$prs.ccgwas)


dat.plot <- dat[,c("Row.names","prs.ccgwas.scale","group")]

p <- ggplot(dat.plot, aes(x = prs.ccgwas.scale, fill=group, color = group)) +
  geom_density(alpha=0.4) +
  xlab("Standardized PRS") +
  ylab("Density") +
  theme_bw() +
  theme(legend.title=element_blank())

p

### Association between PRS and ER risk
dat$prs.ccgwas.tertile <- with(dat, cut(dat$prs.ccgwas, quantile(dat$prs.ccgwas[dat$ER == 0],c(0:3/3)),
                               include.lowest = T, labels=c("Tertile1","Tertile2","Tertile3")))

table(dat$prs.ccgwas.tertile, dat$ER)
re <- glm(ER ~ prs.ccgwas.tertile + PC1_nonde + PC2_nonde + PC3_nonde + PC4_nonde + PC5_nonde + sub_dx_age + as.factor(Phase) + as.factor(stage_fd) + sub_registry, data=dat, family = binomial)
summary(re)
