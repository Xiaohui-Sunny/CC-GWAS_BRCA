# prepare input data for GCTA-COJO
#Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, 
#p-value and sample size. The headers are not keywords and will be omitted by the program. 
#Important: "A1" needs to be the effect allele with "A2" being the other allele and "freq" should be the frequency of "A1".
#example:
#SNP A1 A2 freq b se p N 
#rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
#rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
#rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830

library(data.table)
library(R.utils)
# mernge with original gwas data for freq
bc <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/data/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
bc <- bc[,c("SNP.iCOGs","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs")]
bc$chr_bp <- paste0(bc$chr.iCOGs,"_", bc$Position.iCOGs)

#merge with 1000 genomes obtain the rs for SNP
dbsnp <- fread("/home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur.bim", head=F)
dim(dbsnp)#22665064
dbsnp$chr_bp <- paste0(dbsnp$V1,"_", dbsnp$V4)


#Her2 enriched and LumA
Her2.LumA <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas_result/her2_enrich_lumA.out.results")
Her2.LumA$n <- 
 #cases of her2 and cases of lumA
Her2.LumA$chr_bp <- paste0(Her2.LumA$CHR, '_',  Her2.LumA$BP)
#merge and extract freq
Her2.LumA.m <- merge(bc,Her2.LumA,by.x="chr_bp", by.y="chr_bp")
dim(Her2.LumA.m)#5169082
#merge and extract rs no.
Her2.LumA.m <- merge(Her2.LumA.m, dbsnp,by.x="chr_bp",by.y="chr_bp")
dim(Her2.LumA.m)#5112630
Her2.LumA.m <-Her2.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Her2.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Her2.LumA.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_her2_lumA_gcta_input.ma",row.names = F, quote = F)
rm(Her2.LumA.m)
rm(Her2.LumA) 

#LumB and LumA
LumB.LumA <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas_result/lumB_lumA.out.results")
LumB.LumA$n <- 6659+45436  
LumB.LumA$chr_bp <- paste0(LumB.LumA$CHR, '_',  LumB.LumA$BP)   
LumB.LumA.m <- merge(bc,LumB.LumA,by.x="chr_bp", by.y="chr_bp")
LumB.LumA.m <- merge(LumB.LumA.m, dbsnp,by.x="chr_bp",by.y="chr_bp")
LumB.LumA.m <-LumB.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(LumB.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')   
write.table(LumB.LumA.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_lumB_lumA_gcta_input.ma",row.names = F, quote = F)
rm(LumB.LumA.m)
rm(LumB.LumA) 

#Triple and LumA
Trip.LumA <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas/trip_lumA.out.results")
Trip.LumA$n <- 9067+45436   
Trip.LumA$chr_bp <- paste0(Trip.LumA$CHR, '_',  Trip.LumA$BP)   
Trip.LumA.m <- merge(bc,Trip.LumA,by.x="chr_bp", by.y="chr_bp")
Trip.LumA.m <- merge(Trip.LumA.m, dbsnp,by.x="chr_bp",by.y="chr_bp")
Trip.LumA.m <-Trip.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')   
write.table(Trip.LumA.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumA_gcta_input.ma",row.names = F, quote = F)
rm(Trip.LumA.m)
rm(Trip.LumA) 


#Her2 and LumB
Her2.LumB <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas_result/her2_enrich_lumB.out.results")
Her2.LumB$n <- 3026+6659
Her2.LumB$chr_bp <- paste0(Her2.LumB$CHR, '_',  Her2.LumB$BP)
#merge and extract freq
Her2.LumB.m <- merge(bc,Her2.LumB,by.x="chr_bp", by.y="chr_bp")
Her2.LumB.m <- merge(Her2.LumB.m, dbsnp,by.x="chr_bp",by.y="chr_bp")

Her2.LumB.m <-Her2.LumB.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Her2.LumB.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Her2.LumB.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_her2_lumB_gcta_input.ma",row.names = F, quote = F)
dim(Her2.LumB.m) #5890817
rm(Her2.LumB.m)
rm(Her2.LumB) 


#Triple and LumB
Trip.LumB <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas_result/trip_lumB.out.results")
Trip.LumB$n <- 9067+6659
Trip.LumB$chr_bp <- paste0(Trip.LumB$CHR, '_',  Trip.LumB$BP)
#merge and extract freq
Trip.LumB.m <- merge(bc,Trip.LumB,by.x="chr_bp", by.y="chr_bp")
Trip.LumB.m <- merge(Trip.LumB.m, dbsnp,by.x="chr_bp",by.y="chr_bp")

Trip.LumB.m <-Trip.LumB.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumB.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Trip.LumB.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_gcta_input.ma",row.names = F, quote = F)
dim(Trip.LumB.m) #5938455
rm(Trip.LumB.m)
rm(Trip.LumB) 

# Triple and LumB-Her2-negative
Trip.LumB_her2 <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/ccgwas_result/trip_lumB_her2_neg.out.results")
Trip.LumB_her2$n <- 9067+3026
Trip.LumB_her2$chr_bp <- paste0(Trip.LumB_her2$CHR, '_',  Trip.LumB_her2$BP)
#merge and extract freq
Trip.LumB_her2.m <- merge(bc,Trip.LumB_her2,by.x="chr_bp", by.y="chr_bp")
Trip.LumB_her2.m <- merge(Trip.LumB_her2.m, dbsnp,by.x="chr_bp",by.y="chr_bp")

Trip.LumB_her2.m <-Trip.LumB_her2.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumB_her2.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Trip.LumB_her2.m, "/home/nfs/sunx3/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_her2_neg_gcta_input.ma",row.names = F, quote = F)
dim(Trip.LumB_her2.m) #5974168
rm(Trip.LumB_her2.m)
rm(Trip.LumB_her2) 

#prepare input SNP test list
library(tidyr)
library(dplyr)
#
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_her2_lumA.csv",head=T)
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_her2_lumB.csv",head=T)
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_lumB_lumA.csv",head=T)
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_trip_lumA.csv",head=T)
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_trip_lumB.csv",head=T)
snp.list <-read.csv("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/index and nearby loci_trip_lumB_her2_neg.csv",head=T)
# exclude subjects without pair
#extract index SNP identified from ccgwas
snp.list.index <- filter(snp.list, !duplicated(snp.list$V2))
snp.list.index <- as.data.frame(snp.list.index$V2)

#extract nearby loci(+/- 500kb) associated with overall breast cancer
temp<-as.data.frame(NA)
for(i in 1:nrow(snp.list.index)){
  temp<-as.data.frame(snp.list[which(snp.list$V2==snp.list.index[i,]),c(2,6)])
  snp.list.loci <-c(temp[1,1],temp$SNP)
  write.table(snp.list.loci,paste("G:/My work/Project/BRA subtypes_ccGWAS/data/gcta_test_list2(subtypes)/trip_lumB_her2_neg/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
}
