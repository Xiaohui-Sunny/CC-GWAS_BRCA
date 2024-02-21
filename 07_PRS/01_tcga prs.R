library(data.table)
library(survival)
library(ggplot2)
library(reshape2)

setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/prs_tcga")

dat.cov <- read.csv("TCGA_BRCA_Clinical_With_GenomicRace.csv", head=T)

dat.pca <- fread("TCGA_BRCA_EUR_Female.PC.eigenvec", head=F)
dat.pca <- dat.pca[,c(1,3:7)] #extract the top five CPAs
colnames(dat.pca) <- c("ID","PCA1","PCA2","PCA3","PCA4","PCA5")

### Association between each SNPs identified from CC-GWAS
dat.geno <- read.table("TCGA.PRS.dosage.raw.txt", head=T)
row.names(dat.geno) <- dat.geno$IID
dat.geno <- select(dat.geno,-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")) 
dat.geno <- as.data.frame(t(dat.geno))
dat.geno$SNP <- sapply(strsplit(row.names(dat.geno),"_"),'[',1) 

ccgwasSNPInfo <- read.csv("SNP.ccgwas.csv", head=T)
dat.geno <- dat.geno[dat.geno$SNP %in% ccgwasSNPInfo$Lead.SNP,]#25
dat.geno <- select(dat.geno,-c("SNP"))

dat.geno <- as.data.frame(t(dat.geno))
dat.geno[] <- sapply(dat.geno, as.numeric)

dat.cov$ID <- paste0(dat.cov$bcr_patient_barcode,"-10")
dat <- merge(dat.geno, dat.cov, by.x = "row.names", by.y = "ID")
dat <- merge(dat, dat.pca, by.x = "Row.names", by.y="ID")

table(dat$Genomic.Race)
dat <- subset(dat, dat$Genomic.Race == "WHITE") #743

dat$stage <- ifelse(dat$pathologic_stage=="Stage_I",1,
                    ifelse(dat$pathologic_stage=="Stage_II",2,
                           ifelse(dat$pathologic_stage=="Stage_III",3,
                                  ifelse(dat$pathologic_stage=="Stage_IV",4, 9))))

#### Association between each SNPs identified from CC-GWAS and basal-like risk
## basal vs. lumA
basal.lumA <- subset(dat ,dat$BRCA_Subtype_PAM50_from_mmc4=="Basal" | dat$BRCA_Subtype_PAM50_from_mmc4=="LumA") #527
basal.lumA$outcome <- ifelse(basal.lumA$BRCA_Subtype_PAM50_from_mmc4=="Basal",1,0)
table(basal.lumA$outcome)


SNP_basal.lumA_effect <- matrix(nrow = ncol(dat.geno), ncol = 3)
rownames(SNP_basal.lumA_effect) <- colnames(dat.geno)
colnames(SNP_basal.lumA_effect) <- c("beta", "se", "P")

for ( i in 1:25){
  glmfunc <- as.formula(paste0("outcome ~ ", colnames(dat.geno)[i], "+ PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + age_at_initial_pathologic_diagnosis + as.factor(stage)"))
  re <- glm(glmfunc, data = basal.lumA, family = binomial)
  SNP_basal.lumA_effect[i,] <- summary(re)$coefficients[2,c(1,2,4)]
  
}

write.csv(SNP_basal.lumA_effect,"basal&lumA_snps.csv", row.names = F)


## basal vs. lumB
basal.lumB <- subset(dat, dat$BRCA_Subtype_PAM50_from_mmc4=="Basal" | dat$BRCA_Subtype_PAM50_from_mmc4=="LumB")
basal.lumB$outcome <- ifelse(basal.lumB$BRCA_Subtype_PAM50_from_mmc4=="Basal",1,0)
table(basal.lumB$outcome)

SNP_basal.lumB_effect <- matrix(nrow = ncol(dat.geno), ncol = 3)
rownames(SNP_basal.lumB_effect) <- colnames(dat.geno)
colnames(SNP_basal.lumB_effect) <- c("beta", "se", "P")

for ( i in 1:25){
  glmfunc <- as.formula(paste0("outcome ~ ", colnames(dat.geno)[i], "+ PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + age_at_initial_pathologic_diagnosis + as.factor(stage)"))
  re <- glm(glmfunc, data = basal.lumB, family = binomial)
  SNP_basal.lumB_effect[i,] <- summary(re)$coefficients[2,c(1,2,4)]
  
}

write.csv(SNP_basal.lumB_effect,"basal&lumB_snps.csv", row.names = F)

##################################################################################
# Association between PRS and basal-like risk
dat.geno <- read.table("TCGA.PRS.dosage.raw.txt", head=T)
row.names(dat.geno) <- dat.geno$IID
dat.geno <- select(dat.geno,-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")) 
dat.geno <- as.data.frame(t(dat.geno))


prsSNPInfo <- as.data.frame(t(t(row.names(dat.geno))))
prsSNPInfo$EA.tcga <- sapply(strsplit(prsSNPInfo$V1,"_"),'[',2)
prsSNPInfo$SNP <- sapply(strsplit(prsSNPInfo$V1,"_"),'[',1)

prs.ccgwas.weights <- read.csv("PRS.ccgwas.csv", head=T)
prs.TNBC25.weights <- read.csv("PRS.TNBC25.csv", head=T)
prs.TNBC41.weights <- read.csv("PRS.TNBC41.csv", head=T)

## CCGWAS PRS
prs.ccgwas.weights <- merge(prs.ccgwas.weights, prsSNPInfo, by.x = "SNP", by.y = "SNP")
for (i in 1:nrow(prs.ccgwas.weights)){
  if(prs.ccgwas.weights$EA.tcga[i] == prs.ccgwas.weights$EA[i]) {
    prs.ccgwas.weights$beta.change[i] <- prs.ccgwas.weights$beta[i]
  } else {
    prs.ccgwas.weights$beta.change[i] <- -prs.ccgwas.weights$beta[i]
  }
}

## TNBC25 PRS
prs.TNBC25.weights <- merge(prs.TNBC25.weights, prsSNPInfo, by.x = "SNP", by.y = "SNP")
for (i in 1:nrow(prs.TNBC25.weights)){
  if(prs.TNBC25.weights$EA.tcga[i] == prs.TNBC25.weights$EA[i]) {
    prs.TNBC25.weights$beta.change[i] <- prs.TNBC25.weights$beta[i]
  } else {
    prs.TNBC25.weights$beta.change[i] <- -prs.TNBC25.weights$beta[i]
  }
}


## TNBC41 PRS
prs.TNBC41.weights <- merge(prs.TNBC41.weights, prsSNPInfo, by.x = "SNP", by.y = "SNP")
for (i in 1:nrow(prs.TNBC41.weights)){
  if(prs.TNBC41.weights$EA.tcga[i] == prs.TNBC41.weights$EA[i]) {
    prs.TNBC41.weights$beta.change[i] <- prs.TNBC41.weights$beta[i]
  } else {
    prs.TNBC41.weights$beta.change[i] <- -prs.TNBC41.weights$beta[i]
  }
}


# calculate weighted prs
dat.geno[] <- sapply(dat.geno, as.numeric)

prs.ccgwas <- t(dat.geno[prs.ccgwas.weights$V1,]) %*% prs.ccgwas.weights$beta.change
prs.TNBC25 <- t(dat.geno[prs.TNBC25.weights$V1,]) %*% prs.TNBC25.weights$beta.change
prs.TNBC41 <- t(dat.geno[prs.TNBC41.weights$V1,]) %*% prs.TNBC41.weights$beta.change

allprs <- merge(prs.ccgwas, prs.TNBC25, by = "row.names")
allprs <- merge(allprs, prs.TNBC41, by.x="Row.names", by.y = "row.names")
colnames(allprs) <- c("ID", "prs.ccgwas", "prs.TNBC25","prs.TNBC41")

dat.cov$ID <- paste0(dat.cov$bcr_patient_barcode,"-10")
dat <- merge(allprs, dat.cov, by.x = "ID", by.y = "ID")
dat <- merge(dat, dat.pca, by.x = 'ID', by.y = "ID")
dat <- subset(dat, dat$Genomic.Race == "WHITE") #743

dat$stage <- ifelse(dat$pathologic_stage=="Stage_I",1,
                            ifelse(dat$pathologic_stage=="Stage_II",2,
                                   ifelse(dat$pathologic_stage=="Stage_III",3,
                                          ifelse(dat$pathologic_stage=="Stage_IV",4, 9))))

#################
### Distribution of PRS
basal.lumCa <- subset(dat, dat$BRCA_Subtype_PAM50_from_mmc4=="Basal" | dat$BRCA_Subtype_PAM50_from_mmc4=="LumA" | dat$BRCA_Subtype_PAM50_from_mmc4=="LumB")
basal.lumCa$outcome <- ifelse(basal.lumCa$BRCA_Subtype_PAM50_from_mmc4=="Basal",1,0)

basal.lumCa$group <- ifelse(basal.lumCa$outcome == 0,"LumCa", "Basal")

aggregate(basal.lumCa$prs.ccgwas, list(basal.lumCa$group), FUN=mean) 
aggregate(basal.lumCa$prs.ccgwas, list(basal.lumCa$group), FUN=sd) 
t.test(basal.lumCa$prs.ccgwas ~ basal.lumCa$group)

t.test(basal.lumCa$prs.TNBC25 ~ basal.lumCa$group)
t.test(basal.lumCa$prs.TNBC41 ~ basal.lumCa$group)



basal.lumCa$prs.ccgwas.scale <- basal.lumCa$prs.ccgwas/sd(basal.lumCa$prs.ccgwas)
basal.lumCa$prs.TNBC25.scale <- basal.lumCa$prs.TNBC25/sd(basal.lumCa$prs.TNBC25)
basal.lumCa$prs.TNBC41.scale <- basal.lumCa$prs.TNBC41/sd(basal.lumCa$prs.TNBC41)

dat.plot <- basal.lumCa[,c("group","prs.ccgwas.scale","prs.TNBC25.scale","prs.TNBC41.scale")]
dat.plot <- melt(dat.plot, "group")


p <- ggplot(dat.plot, aes(x = value, fill=group, color = group)) +
  geom_density(alpha=0.4) +
  facet_wrap( ~ variable, labeller = as_labeller(c(prs.ccgwas.scale = "PRS of cc-GWAS",
                                                   prs.TNBC25.scale = "PRS25 from Zhang et al.",
                                                   prs.TNBC41.scale = "PRS41 from Zhang et al."))) +
  xlab("Standardized PRS") +
  ylab("Density") + 
  theme_bw() +
  theme(legend.title=element_blank()) +
  theme(strip.text = element_text(size = 12, face = "bold"))
p



basal.lumCa$prs.ccgwas.tertile <- with(basal.lumCa, cut(basal.lumCa$prs.ccgwas, quantile(basal.lumCa$prs.ccgwas[basal.lumCa$outcome == 0],c(0:3/3)),
                                  include.lowest =T, labels=c("Tertile1","Tertile2","Tertile3")))



re.dat <- basal.lumCa %>%
  tabyl(outcome, prs.ccgwas.tertile) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns(position = "front")
re.dat


basal.lumCa$prs.TNBC25.tertile <- with(basal.lumCa, cut(basal.lumCa$prs.TNBC25, quantile(basal.lumCa$prs.TNBC25[basal.lumCa$outcome == 0],c(0:3/3)),
                                       include.lowest = T, labels= c("Tertile1","Tertile2","Tertile3")))

re.dat <- basal.lumCa %>%
  tabyl(outcome, prs.TNBC25.tertile) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns(position = "front")
re.dat



basal.lumCa$prs.TNBC41.tertile <- with(basal.lumCa, cut(basal.lumCa$prs.TNBC41, quantile(basal.lumCa$prs.TNBC41[basal.lumCa$outcome == 0],c(0:3/3)),
                                       include.lowest = T, labels= c("Tertile1","Tertile2","Tertile3")))


re.dat <- basal.lumCa %>%
  tabyl(outcome, prs.TNBC41.tertile) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns(position = "front")
re.dat


#########################
### Association between PRS and TNBC risk
#### Model 1: adjusted age, stage and top five PCAs
PRS_model1_effect <- matrix(nrow = 6, ncol = 3)
rownames(PRS_model1_effect) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)",
                              "Model 1(prs.TNBC25.tertile2)","Model 1(prs.TNBC25.tertile3)",
                              "Model 1(prs.TNBC41.tertile2)","Model 1(prs.TNBC41.tertile3)")

colnames(PRS_model1_effect) <- c("beta", "se", "P")

re1 <- glm(outcome ~ prs.ccgwas.tertile + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + as.factor(stage), data = basal.lumCa, family = binomial)
PRS_model1_effect[1:2,] <- summary(re1)$coefficients[2:3,c(1,2,4)]

re2 <- glm(outcome ~ prs.TNBC25.tertile + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + as.factor(stage), data = basal.lumCa, family = binomial)
PRS_model1_effect[3:4,] <- summary(re2)$coefficients[2:3,c(1,2,4)]

re3 <- glm(outcome ~ prs.TNBC41.tertile + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + as.factor(stage), data = basal.lumCa, family = binomial)
PRS_model1_effect[5:6,] <- summary(re3)$coefficients[2:3,c(1,2,4)]

#### Model 2: adjusted age, stage,top five PCAs and PRS.TNBC25/PRS.TNBC41
PRS_model2_effect <- matrix(nrow = 8, ncol = 3)
rownames(PRS_model2_effect) <- c("Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)",
                              "Model 2(prs.TNBC25.tertile2)","Model 2(prs.TNBC25.tertile3)",
                              "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)",
                              "Model 2(prs.TNBC41.tertile2)","Model 2(prs.TNBC41.tertile3)")

colnames(PRS_model2_effect) <- c("beta", "se", "P")

re1 <- glm(outcome ~ prs.ccgwas.tertile + prs.TNBC25.tertile + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + as.factor(stage), data = basal.lumCa, family = binomial)
PRS_model2_effect[1:4,] <- summary(re1)$coefficients[2:5,c(1,2,4)]

re2 <- glm(outcome ~ prs.ccgwas.tertile + prs.TNBC41.tertile + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5 + as.factor(stage), data = basal.lumCa, family = binomial)
PRS_model2_effect[5:8,] <- summary(re2)$coefficients[2:5,c(1,2,4)]


###########################
### Association between PRS and survival(OS and DSS)
#### All subjects
dat$prs.ccgwas.tertile <- with(dat, cut(dat$prs.ccgwas, quantile(dat$prs.ccgwas, c(0:3/3)),
                                   include.lowest = T, labels=c("Tertile1","Tertile2","Tertile3")))

dat$prs.ccgwas.tertile.comb <- ifelse(dat$prs.ccgwas.tertile == "Tertile1", "Tertile1","Tertile2&3")

##### OS
dat$OS.time <- as.numeric(dat$OS.time)
dat$age_at_initial_pathologic_diagnosis <- as.numeric(dat$age_at_initial_pathologic_diagnosis)
Cox_OS_all <- matrix(nrow = 9, ncol = 4)
rownames(Cox_OS_all) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                          "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)",
                          "Model 3(prs.ccgwas.tertile2)","Model 3(prs.ccgwas.tertile3)","Model 3(prs.ccgwas.tertile2&3)")

colnames(Cox_OS_all) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
re2.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
Cox_OS_all[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_OS_all[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_OS_all[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_OS_all[3,4] <- summary(re2.model1)$coefficients[1,5] 



#model2: adjusted for age, stage, PCAs, PAM50, PRS.TNBC25
re1.model3 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + BRCA_Subtype_PAM50_from_mmc4 + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
re2.model3 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + BRCA_Subtype_PAM50_from_mmc4 + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
Cox_OS_all[7:8,1:3] <- summary(re1.model3)$conf.int[1:2,c(1,3,4)] 
Cox_OS_all[7:8,4] <- summary(re1.model3)$coefficients[1:2,5]
Cox_OS_all[9,1:3] <- summary(re2.model3)$conf.int[1,c(1,3,4)] 
Cox_OS_all[9,4] <- summary(re2.model3)$coefficients[1,5] 

##### DSS
dat$DSS.time <- as.numeric(dat$DSS.time)
dat$DSS <- as.numeric(dat$DSS)
Cox_DSS_all <- matrix(nrow = 9, ncol = 4)
rownames(Cox_DSS_all) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                          "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)",
                          "Model 3(prs.ccgwas.tertile2)","Model 3(prs.ccgwas.tertile3)","Model 3(prs.ccgwas.tertile2&3)")

colnames(Cox_DSS_all) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
re2.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
Cox_DSS_all[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_all[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_DSS_all[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_DSS_all[3,4] <- summary(re2.model1)$coefficients[1,5] 


#model2: adjusted for age, stage, PCAs, PAM50, PRS.TNBC25
re1.model3 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + BRCA_Subtype_PAM50_from_mmc4 + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
re2.model3 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + BRCA_Subtype_PAM50_from_mmc4 + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat)
Cox_DSS_all[7:8,1:3] <- summary(re1.model3)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_all[7:8,4] <- summary(re1.model3)$coefficients[1:2,5]
Cox_DSS_all[9,1:3] <- summary(re2.model3)$conf.int[1,c(1,3,4)] 
Cox_DSS_all[9,4] <- summary(re2.model3)$coefficients[1,5] 


#### Basal subjects
dat.basal <- subset(dat, dat$BRCA_Subtype_PAM50_from_mmc4=="Basal") #108
dat.basal$prs.ccgwas.tertile <- with(dat.basal, cut(dat.basal$prs.ccgwas, quantile(dat.basal$prs.ccgwas, c(0:3/3)),
                                        include.lowest = T, labels=c("Tertile1","Tertile2","Tertile3")))

dat.basal$prs.ccgwas.tertile.comb <- ifelse(dat.basal$prs.ccgwas.tertile == "Tertile1", "Tertile1","Tertile2&3")

##### OS
dat.basal$OS.time <- as.numeric(dat.basal$OS.time)
dat.basal$age_at_initial_pathologic_diagnosis <- as.numeric(dat.basal$age_at_initial_pathologic_diagnosis)
Cox_OS_basal <- matrix(nrow = 6, ncol = 4)
rownames(Cox_OS_basal) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                          "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)")

colnames(Cox_OS_basal) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
re2.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
Cox_OS_basal[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_OS_basal[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_OS_basal[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_OS_basal[3,4] <- summary(re2.model1)$coefficients[1,5] 


#model2: adjusted for age, stage, PCAs, PRS.TNBC25
re1.model2 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
re2.model2 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
Cox_OS_basal[4:5,1:3] <- summary(re1.model2)$conf.int[1:2,c(1,3,4)] 
Cox_OS_basal[4:5,4] <- summary(re1.model2)$coefficients[1:2,5]
Cox_OS_basal[6,1:3] <- summary(re2.model2)$conf.int[1,c(1,3,4)] 
Cox_OS_basal[6,4] <- summary(re2.model2)$coefficients[1,5] 

##### DSS
Cox_DSS_basal <- matrix(nrow = 6, ncol = 4)
rownames(Cox_DSS_basal) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                            "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)")

colnames(Cox_DSS_basal) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
re2.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
Cox_DSS_basal[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_basal[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_DSS_basal[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_DSS_basal[3,4] <- summary(re2.model1)$coefficients[1,5] 


#model2: adjusted for age, stage, PCAs, PRS.TNBC25
re1.model2 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
re2.model2 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.basal)
Cox_DSS_basal[4:5,1:3] <- summary(re1.model2)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_basal[4:5,4] <- summary(re1.model2)$coefficients[1:2,5]
Cox_DSS_basal[6,1:3] <- summary(re2.model2)$conf.int[1,c(1,3,4)] 
Cox_DSS_basal[6,4] <- summary(re2.model2)$coefficients[1,5] 


#### LumA+lumB subjects
dat.lum <- subset(dat, dat$BRCA_Subtype_PAM50_from_mmc4=="LumA" |dat$BRCA_Subtype_PAM50_from_mmc4=="LumB") #565
dat.lum$prs.ccgwas.tertile <- with(dat.lum, cut(dat.lum$prs.ccgwas, quantile(dat.lum$prs.ccgwas, c(0:3/3)),
                                                    include.lowest = T, labels=c("Tertile1","Tertile2","Tertile3")))

dat.lum$prs.ccgwas.tertile.comb <- ifelse(dat.lum$prs.ccgwas.tertile == "Tertile1", "Tertile1","Tertile2&3")
##### OS
Cox_OS_lum <- matrix(nrow = 6, ncol = 4)
rownames(Cox_OS_lum) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                            "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)")

colnames(Cox_OS_lum) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
re2.model1 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
Cox_OS_lum[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_OS_lum[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_OS_lum[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_OS_lum[3,4] <- summary(re2.model1)$coefficients[1,5] 


#model2: adjusted for age, stage, PCAs, PRS.TNBC25
re1.model2 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
re2.model2 <- coxph(Surv(OS.time, OS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
Cox_OS_lum[4:5,1:3] <- summary(re1.model2)$conf.int[1:2,c(1,3,4)] 
Cox_OS_lum[4:5,4] <- summary(re1.model2)$coefficients[1:2,5]
Cox_OS_lum[6,1:3] <- summary(re2.model2)$conf.int[1,c(1,3,4)] 
Cox_OS_lum[6,4] <- summary(re2.model2)$coefficients[1,5] 

##### DSS
Cox_DSS_lum <- matrix(nrow = 6, ncol = 4)
rownames(Cox_DSS_lum) <- c("Model 1(prs.ccgwas.tertile2)","Model 1(prs.ccgwas.tertile3)","Model 1(prs.ccgwas.tertile2&3)",
                          "Model 2(prs.ccgwas.tertile2)","Model 2(prs.ccgwas.tertile3)","Model 2(prs.ccgwas.tertile2&3)")

colnames(Cox_DSS_lum) <- c("HR", "95.l","95.u","P")

#model1: adjusted for age, stage ,PCAs
re1.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
re2.model1 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
Cox_DSS_lum[1:2,1:3] <- summary(re1.model1)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_lum[1:2,4] <- summary(re1.model1)$coefficients[1:2,5]
Cox_DSS_lum[3,1:3] <- summary(re2.model1)$conf.int[1,c(1,3,4)] 
Cox_DSS_lum[3,4] <- summary(re2.model1)$coefficients[1,5] 


#model2: adjusted for age, stage, PCAs, PRS.TNBC25
re1.model2 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
re2.model2 <- coxph(Surv(DSS.time, DSS) ~ prs.ccgwas.tertile.comb + as.factor(stage) + age_at_initial_pathologic_diagnosis + prs.TNBC25 + PCA1 + PCA2 + PCA3 +PCA4 + PCA5, data = dat.lum)
Cox_DSS_lum[4:5,1:3] <- summary(re1.model2)$conf.int[1:2,c(1,3,4)] 
Cox_DSS_lum[4:5,4] <- summary(re1.model2)$coefficients[1:2,5]
Cox_DSS_lum[6,1:3] <- summary(re2.model2)$conf.int[1,c(1,3,4)] 
Cox_DSS_lum[6,4] <- summary(re2.model2)$coefficients[1,5] 
