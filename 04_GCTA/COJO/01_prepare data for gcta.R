# prepare input data for GCTA-COJO
#Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, 
#p-value and sample size. The headers are not keywords and will be omitted by the program. 

library(data.table)
library(R.utils)

# mernge with original gwas data for freq
bc <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
bc <- bc[,c("SNP.iCOGs","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs")]
bc$chr_bp <- paste0(bc$chr.iCOGs,"_", bc$Position.iCOGs)

# merge with 1000 genomes to obtain the rs for SNP
geno1000 <- fread("/home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur.bim", head=F)
dim(geno1000)#22665064
geno1000$chr_bp <- paste0(geno1000$V1,"_", geno1000$V4)


## Her2 enriched and LumA
Her2.LumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/her2_enrich_lumA.out.results.gz")
Her2.LumA$n <- 5859+60715 #cases of her2 and cases of lumA
Her2.LumA$chr_bp <- paste0(Her2.LumA$CHR, '_',  Her2.LumA$BP)
# merge and extract freq
Her2.LumA.m <- merge(bc,Her2.LumA,by.x="chr_bp", by.y="chr_bp")
dim(Her2.LumA.m)#5982435
# merge and extract rs no.
Her2.LumA.m <- merge(Her2.LumA.m, geno1000,by.x="chr_bp",by.y="chr_bp")
dim(Her2.LumA.m)#5915848
Her2.LumA.m <-Her2.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Her2.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Her2.LumA.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_her2_lumA_gcta_input0428.ma",row.names = F, quote = F)
rm(Her2.LumA.m)
rm(Her2.LumA) 

## LumB and LumA
LumB.LumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/lumB_lumA.out.results.gz")
LumB.LumA$n <- 12414+60715  
LumB.LumA$chr_bp <- paste0(LumB.LumA$CHR, '_',  LumB.LumA$BP)   
LumB.LumA.m <- merge(bc,LumB.LumA,by.x="chr_bp", by.y="chr_bp")
LumB.LumA.m <- merge(LumB.LumA.m, geno1000,by.x="chr_bp",by.y="chr_bp")
LumB.LumA.m <-LumB.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(LumB.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')   
write.table(LumB.LumA.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_lumB_lumA_gcta_input0428.ma",row.names = F, quote = F)
rm(LumB.LumA.m)
rm(LumB.LumA) 

## Triple and LumA
Trip.LumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/trip_lumA.out.results.gz")
Trip.LumA$n <- 16499+60715  
Trip.LumA$chr_bp <- paste0(Trip.LumA$CHR, '_',  Trip.LumA$BP)   
Trip.LumA.m <- merge(bc,Trip.LumA,by.x="chr_bp", by.y="chr_bp")
Trip.LumA.m <- merge(Trip.LumA.m, geno1000,by.x="chr_bp",by.y="chr_bp")
Trip.LumA.m <-Trip.LumA.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumA.m) <- c('SNP','A1','A2','freq','b','se','p','N')   
write.table(Trip.LumA.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumA_gcta_input0428.ma",row.names = F, quote = F)
rm(Trip.LumA.m)
rm(Trip.LumA) 

## Her2 and LumB
Her2.LumB <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/her2_enrich_lumB.out.results.gz")
Her2.LumB$n <- 5859+12414
Her2.LumB$chr_bp <- paste0(Her2.LumB$CHR, '_',  Her2.LumB$BP)
#merge and extract freq
Her2.LumB.m <- merge(bc,Her2.LumB,by.x="chr_bp", by.y="chr_bp")
Her2.LumB.m <- merge(Her2.LumB.m, geno1000,by.x="chr_bp",by.y="chr_bp")
Her2.LumB.m <- Her2.LumB.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Her2.LumB.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Her2.LumB.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_her2_lumB_gcta_input0428.ma",row.names = F, quote = F)
dim(Her2.LumB.m) #5892359
rm(Her2.LumB.m)
rm(Her2.LumB) 


## Triple and LumB
Trip.LumB <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/trip_lumB.out.results.gz")
Trip.LumB$n <- 16499+12414
Trip.LumB$chr_bp <- paste0(Trip.LumB$CHR, '_',  Trip.LumB$BP)
#merge and extract freq
Trip.LumB.m <- merge(bc,Trip.LumB,by.x="chr_bp", by.y="chr_bp")
Trip.LumB.m <- merge(Trip.LumB.m, geno1000,by.x="chr_bp",by.y="chr_bp")
Trip.LumB.m <-Trip.LumB.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumB.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Trip.LumB.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_gcta_input0428.ma",row.names = F, quote = F)
dim(Trip.LumB.m) #5938363
rm(Trip.LumB.m)
rm(Trip.LumB) 

## Triple and LumB-Her2-negative
Trip.LumB_her2 <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/trip_lumB_her2_neg.out.results.gz")
Trip.LumB_her2$n <- 16499+16943
Trip.LumB_her2$chr_bp <- paste0(Trip.LumB_her2$CHR, '_',  Trip.LumB_her2$BP)
#merge and extract freq
Trip.LumB_her2.m <- merge(bc,Trip.LumB_her2,by.x="chr_bp", by.y="chr_bp")
Trip.LumB_her2.m <- merge(Trip.LumB_her2.m, geno1000,by.x="chr_bp",by.y="chr_bp")
Trip.LumB_her2.m <-Trip.LumB_her2.m[,c("V2","EA","NEA","EAFcontrols.iCOGs","OLS_beta","OLS_se","OLS_pval","n")]
colnames(Trip.LumB_her2.m) <- c('SNP','A1','A2','freq','b','se','p','N')
write.table(Trip.LumB_her2.m, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_her2_neg_gcta_input0428.ma",row.names = F, quote = F)
dim(Trip.LumB_her2.m) #5972322
rm(Trip.LumB_her2.m)
rm(Trip.LumB_her2) 


#Triple and Her2: 0 snp has been identified 


#Her2 and LumB-Her2-neg


###################################################################################################################
##################################################################################################################
#prepare input SNP test list
library(tidyr)
library(dplyr)

## prepare the snps associated with overall snps from 
snp.list.overallBC <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/pre_loci.csv",head=T)
snp.list.overallBC$chr_bp <- paste0(snp.list.overallBC$Chr, '_',  snp.list.overallBC$Position)
snp.list.overallBC <- filter(snp.list.overallBC,!duplicated(snp.list.overallBC$chr_bp)) #remove the duplicated SNP
dim(snp.list.overallBC) #209
colnames(snp.list.overallBC) <- c('SNP', 'CHR', 'BP','chr_bp')

bc <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")
bc$chr_bp <- paste0(bc$chr.iCOGs,"_", bc$Position.iCOGs)
### merge with the overall BC summary statistics to extract the meta P
snp.list.overallBC.m <- merge(snp.list.overallBC, bc, by.x="chr_bp",by.y="chr_bp")
snp.list.overallBC.m <- filter(snp.list.overallBC.m,!duplicated(snp.list.overallBC.m$chr_bp)) #remove the duplicated SNPs
dim(snp.list.overallBC.m) #206


overallBC.notmerge <- filter(snp.list.overallBC, !(snp.list.overallBC$SNP %in% snp.list.overallBC.m$SNP))
### find proxy for these 3 SNPs
# rs554219 no proxy
## 	rs11571833: 13_32968550 proxy can not find in the summary data
##  rs17879961: 22_29008888, 22_29098376 proxy can not find in the summary data
a <- bc[which(bc$chr_bp=='22_29098376'),]


##################################################################################
##################################################################################
# extract nearby loci(+/- 500kb) associated with overall breast cancer
## Step 1 : extract nearby SNPs associated with overall Bca (378 SNPs) (+/- 500kb) for each index SNP 
## Step 2: examine if independent nearby snps can be merged with the summary statistics
## Step 3: find a proxy for the missing SNP
## Step 4: generate the .txt file for each index snp and its nearby(+/- 500 kb) SNPs associated with overall breast cancer
## Step 5: add the proxy to the snp list and generate the input data for gcta-cojo

snp.list.ccgwas <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/BCAC_subtype_ccgwas_sugg_LD_prune_0605.csv",head=T)
snp.list.ccgwas <- subset(snp.list.ccgwas,snp.list.ccgwas$index==1,) #29
dim(snp.list.ccgwas)
snp.list.ccgwas.m <- merge(snp.list.ccgwas,geno1000,by.x="chr_bp", by.y="chr_bp") #extract the snp id for each index snp

snp.list.ccgwas <- snp.list.ccgwas.m[,c('V2.x','CHR','BP','cc')]
colnames(snp.list.ccgwas)[1] <- c('SNP')

### triple and lumA
trip.lumA.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='trip_lumA.out.results.gz') #21index SNPs
trip.lumA.snp.BC.all <- c()
snp.list.overallBC.m <- as.data.frame(snp.list.overallBC.m)

for (i in 1:nrow(trip.lumA.snp.list.ccgwas)){
  snp <- subset(snp.list.overallBC.m, CHR==trip.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumA.snp.list.ccgwas$BP[i]-500000)
  snp$ccgwas.snp <- trip.lumA.snp.list.ccgwas$SNP[i]
  snp <- snp[,c('SNP','CHR','BP','ccgwas.snp')]
  trip.lumA.snp.BC.all <- rbind(trip.lumA.snp.BC.all, snp)
}


dim(trip.lumA.snp.BC.all) # 27 independent snps associated overall breast cancer have been found for 21 index SNPs 
trip.lumA.snp.BC.all <- as.data.frame(trip.lumA.snp.BC.all)
trip.lumA.snp.BC.undup <- filter(trip.lumA.snp.BC.all, !duplicated(trip.lumA.snp.BC.all$SNP)) #19

trip.lumA.sum <- fread ("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumA_gcta_input0428.ma")
trip.lumA.snp.BC.m <- merge(trip.lumA.snp.BC.undup, trip.lumA.sum, by.x="SNP", by.y="SNP") # examine if the snps can be found in the summary statistics
dim(trip.lumA.snp.BC.m) #11 index  in summary data

trip.lumA.nomerge <- filter(trip.lumA.snp.BC.undup, !(trip.lumA.snp.BC.undup$SNP %in% trip.lumA.snp.BC.m$SNP)) 
#13 can not be found in the ccgwas summary statistics
write.csv(trip.lumA.snp.BC.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumA.snp.overall.0605.csv",row.names = F)

write.csv(trip.lumA.nomerge,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumA.snp.overall.nofind.0605.csv",row.names = F)

trip.lumA.proxy <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/triple.lumA.proxy.csv")
trip.lumA.proxy.m <- merge(trip.lumA.proxy, trip.lumA.sum, by.x="Proxy", by.y="SNP")

# rs11814448 proxy can not be found
# rs11374964 no proxy
# rs74911261: 
proxy0 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs149934734"),]
# rs67397200
proxy1 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs61494113"),]
# rs6678914:
proxy2 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs58320062"),]
# rs10069690: 
proxy3 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs2242652"),]
# rs3215401
proxy5 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs3824987"),]
# rs514192: 
proxy4 <- trip.lumA.sum[which(trip.lumA.sum$SNP=="rs578355"),]




trip.lumA.snp.BC.all[nrow(trip.lumA.snp.BC.all)+1,] <- c('rs61494113',19,17401859,'rs67397200')
trip.lumA.snp.BC.all[nrow(trip.lumA.snp.BC.all)+1,] <- c('rs58320062',1,202181341,'rs6678914')
trip.lumA.snp.BC.all[nrow(trip.lumA.snp.BC.all)+1,] <- c('rs578355',8,102483098,'rs514192')
trip.lumA.snp.BC.all[nrow(trip.lumA.snp.BC.all)+1,] <- c('rs2242652',5,1280028,'rs10069690')
trip.lumA.snp.BC.all[nrow(trip.lumA.snp.BC.all)+1,] <- c('rs3824987',11,108346103,'rs10069690')



trip.lumA.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='trip_lumA.out.results.gz') #39 index SNPs

snp.list.all <- c()
for(i in 1:nrow(trip.lumA.snp.list.ccgwas)){
  snp.list <- c()
  snp <- c()
  snp <- subset(trip.lumA.snp.BC.all, CHR==trip.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumA.snp.list.ccgwas$BP[i]-500000)
  snp <- filter(snp, !duplicated(snp$SNP))
  colnames(snp)[1] <- c("SNP")
  snp.list <- rbind(trip.lumA.snp.list.ccgwas[i,c(1:3)],snp[,c(1:3)])
  snp.list$no.test <- paste0("test",i)
  snp.list$ccgwas.snp <- trip.lumA.snp.list.ccgwas[i,1]
  write.table(snp.list[,c(1:3)],paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumA/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
  snp.list.all <- rbind(snp.list.all,snp.list)
  
}
write.table(snp.list.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumA/test.all.txt")


#################################################################################
##### triple and lumB
trip.lumB.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='trip_lumB.out.results.gz')#4
trip.lumB.snp.BC.all <- c()
for (i in 1:nrow(trip.lumB.snp.list.ccgwas)){
  snp <- subset(snp.list.overallBC.m, CHR==trip.lumB.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumB.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumB.snp.list.ccgwas$BP[i]-500000)
  snp$ccgwas.snp <- trip.lumB.snp.list.ccgwas$SNP[i]
  snp <- snp[,c('SNP','CHR','BP','ccgwas.snp')]
  trip.lumB.snp.BC.all <- rbind(trip.lumB.snp.BC.all, snp)
}

dim(trip.lumB.snp.BC.all) #8 independent snps have been found for 4 index SNPs
trip.lumB.snp.BC.all <- as.data.frame(trip.lumB.snp.BC.all)
trip.lumB.snp.BC.undup <- filter(trip.lumB.snp.BC.all, !duplicated(trip.lumB.snp.BC.all$SNP)) #8
dim(trip.lumB.snp.BC.undup)
trip.lumB.sum <- fread ("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_gcta_input0428.ma")
trip.lumB.snp.BC.m <- merge(trip.lumB.snp.BC.undup, trip.lumB.sum, by.x="SNP", by.y="SNP") # examine if the snps can be found in the summary statistics
dim(trip.lumB.snp.BC.m) #4

trip.lumB.nomerge <- filter(trip.lumB.snp.BC.undup, !(trip.lumB.snp.BC.undup$SNP %in% trip.lumB.snp.BC.m$SNP)) #4
write.csv(trip.lumB.snp.BC.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumB.snp.overall.0605.csv",row.names = F)

write.csv(trip.lumB.nomerge,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumB.snp.overall.nofind.0605.csv",row.names = F)

# merge with proxy
trip.lumB.proxy <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/triple.lumB.proxy.csv")
trip.lumB.proxy.m <- merge(trip.lumB.proxy, trip.lumB.sum, by.x="Proxy", by.y="SNP")


# rs35054928: no proxy
# rs45631563
proxy <- trip.lumB.sum[which(trip.lumB.sum$SNP=="rs45631582"),] #r2=0.83

#rs11374964: no proxy
proxy2 <- trip.lumB.sum[which(trip.lumB.sum$SNP=="rs3824987"),]
#rs74911261: proxy can not be found


trip.lumB.snp.BC.all[nrow(trip.lumB.snp.BC.all)+1,] <- c('rs45631582',10,123343544,'rs45631563')
trip.lumB.snp.BC.all[nrow(trip.lumB.snp.BC.all)+1,] <- c('rs3824987',11,108346103,'rs10069690')


snp.list.all <- c()
for(i in 1:nrow(trip.lumB.snp.list.ccgwas)){
  snp.list <- c()
  snp <- subset(trip.lumB.snp.BC.all, CHR==trip.lumB.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumB.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumB.snp.list.ccgwas$BP[i]-500000)
  snp <- filter(snp, !duplicated(snp$SNP))
  colnames(snp)[1] <- c("SNP")
  snp.list <- rbind(trip.lumB.snp.list.ccgwas[i,c(1:3)],snp[,c(1:3)])
  snp.list$no.test <- paste0("test",i)
  snp.list$ccgwas.snp <- trip.lumB.snp.list.ccgwas[i,1]
  write.table(snp.list[,c(1:3)],paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumB/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
  snp.list.all <- rbind(snp.list.all,snp.list)
  
}
write.table(snp.list.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumB/test.all.txt")


####################################
##### triple and lumB-her2-neg
trip.lumB.her2.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='trip_lumB_her2_neg.out.results.gz')#2
trip.lumB.her2.snp.BC.all <- c()
for (i in 1:nrow(trip.lumB.her2.snp.list.ccgwas)){
  snp <- subset(snp.list.overallBC.m, CHR==trip.lumB.her2.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumB.her2.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumB.her2.snp.list.ccgwas$BP[i]-500000)
  snp$ccgwas.snp <- trip.lumB.her2.snp.list.ccgwas$SNP[i]
  snp <- snp[,c('SNP','CHR','BP','ccgwas.snp')]
  trip.lumB.her2.snp.BC.all <- rbind(trip.lumB.her2.snp.BC.all, snp)
}

dim(trip.lumB.her2.snp.BC.all) #3 independent snps have been found for 2 index SNPs
trip.lumB.her2.snp.BC.all <- as.data.frame(trip.lumB.her2.snp.BC.all)
trip.lumB.her2.snp.BC.undup <- filter(trip.lumB.her2.snp.BC.all, !duplicated(trip.lumB.her2.snp.BC.all$SNP)) #3
dim(trip.lumB.her2.snp.BC.undup)
trip.lumB.her2.sum <- fread ("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_trip_lumB_her2_neg_gcta_input0428.ma")

trip.lumB.her2.snp.BC.m <- merge(trip.lumB.her2.snp.BC.undup,trip.lumB.her2.sum,by.x="SNP", by.y="SNP") # examine if the snps can be found in the summary statistics
dim(trip.lumB.her2.snp.BC.m) #1

trip.lumB.her2.nomerge <- filter(trip.lumB.her2.snp.BC.undup, !(trip.lumB.her2.snp.BC.undup$SNP %in% trip.lumB.her2.snp.BC.m$SNP)) #1
dim(trip.lumB.her2.nomerge)
write.csv(trip.lumB.her2.snp.BC.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumB.her2.snp.overall.0605.csv",row.names = F)
write.csv(trip.lumB.her2.nomerge,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumB.her2.snp.overall.nofind.0605.csv",row.names = F)

# merge with proxy
trip.lumB.her2.proxy <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/trip.lumB.her2.proxy.csv")
trip.lumB.her2.proxy.m <- merge(trip.lumB.her2.proxy, trip.lumB.her2.sum, by.x="Proxy", by.y="SNP")

#rs11814448: proxy can not be found

#rs514192:proxy can not be found
proxy1 <- trip.lumB.her2.sum[which(trip.lumB.her2.sum$SNP=='rs578355'),]

trip.lumB.her2.snp.BC.all[nrow(trip.lumB.her2.snp.BC.all)+1,] <- c('rs578355',8,102483098,'rs514192')


snp.list.all <- c()
for(i in 1:nrow(trip.lumB.her2.snp.list.ccgwas)){
  snp.list <- c()
  snp <- subset(trip.lumB.her2.snp.BC.all, CHR==trip.lumB.her2.snp.list.ccgwas$CHR[i] & as.numeric(BP) < trip.lumB.her2.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > trip.lumB.her2.snp.list.ccgwas$BP[i]-500000)
  snp <- filter(snp, !duplicated(snp$SNP))
  colnames(snp)[1] <- c("SNP")
  snp.list <- rbind(trip.lumB.her2.snp.list.ccgwas[i,c(1:3)],snp[,c(1:3)])
  snp.list$no.test <- paste0("test",i)
  snp.list$ccgwas.snp <- trip.lumB.her2.snp.list.ccgwas[i,1]
  write.table(snp.list[,c(1:3)],paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumB_her2_neg/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
  snp.list.all <- rbind(snp.list.all,snp.list)
  
}
write.table(snp.list.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/trip_lumB_her2_neg/test.all.txt")

################################################################################
##### lumB and lumA
lumB.lumA.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='lumB_lumA.out.results.gz') #1
lumB.lumA.snp.BC.all <- c()
for (i in 1:nrow(lumB.lumA.snp.list.ccgwas)){
  snp <- subset(snp.list.overallBC.m, CHR==lumB.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < lumB.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > lumB.lumA.snp.list.ccgwas$BP[i]-500000)
  snp$ccgwas.snp <- lumB.lumA.snp.list.ccgwas$SNP[i]
  snp <- snp[,c('SNP','CHR','BP','ccgwas.snp')]
  lumB.lumA.snp.BC.all <- rbind(lumB.lumA.snp.BC.all, snp)
}
dim(lumB.lumA.snp.BC.all)#1
lumB.lumA.snp.BC.all <- as.data.frame(lumB.lumA.snp.BC.all)
lumB.lumA.snp.BC.undup <- filter(lumB.lumA.snp.BC.all, !duplicated(lumB.lumA.snp.BC.all$SNP)) #1


lumB.lumA.sum <- fread ("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_lumB_lumA_gcta_input0428.ma")
lumB.lumA.snp.BC.m <- merge(lumB.lumA.snp.BC.undup,lumB.lumA.sum,by.x="SNP", by.y="SNP") # examine if the snps can be found in the summary statistics
dim(lumB.lumA.snp.BC.m) #1

lumB.lumA.nomerge <- filter(lumB.lumA.snp.BC.undup, !(lumB.lumA.snp.BC.undup$SNP %in% lumB.lumA.snp.BC.m$SNP))
write.csv(lumB.lumA.snp.BC.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/lumB.lumA.snp.overall.0605.csv",row.names = F)
write.csv(lumB.lumA.nomerge,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/lumB.lumA.snp.overall.nofind.0605.csv",row.names = F)


snp.list.all <- c()
for(i in 1:nrow(lumB.lumA.snp.list.ccgwas)){
  snp.list <- c()
  snp <- subset(lumB.lumA.snp.BC.all, CHR==lumB.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < lumB.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > lumB.lumA.snp.list.ccgwas$BP[i]-500000)
  snp <- filter(snp, !duplicated(snp$SNP))
  colnames(snp)[1] <- c("SNP")
  snp.list <- rbind(lumB.lumA.snp.list.ccgwas[i,c(1:3)],snp[,c(1:3)])
  snp.list$no.test <- paste0("test",i)
  snp.list$ccgwas.snp <- lumB.lumA.snp.list.ccgwas[i,1]
  write.table(snp.list[,c(1:3)],paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/lumB_lumA/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
  snp.list.all <- rbind(snp.list.all,snp.list)
  
}
write.table(snp.list.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/lumB_lumA/test.all.txt")


###############################################################################
##### her2 and lumA
her2.lumA.snp.list.ccgwas <- subset(snp.list.ccgwas, snp.list.ccgwas$cc=='her2_enrich_lumA.out.results.gz')#1
her2.lumA.snp.BC.all <- c()
for (i in 1:nrow(her2.lumA.snp.list.ccgwas)){
  snp <- subset(snp.list.overallBC.m, CHR==her2.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < her2.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > her2.lumA.snp.list.ccgwas$BP[i]-500000)
  snp$ccgwas.snp <- her2.lumA.snp.list.ccgwas$SNP[i]
  snp <- snp[,c('SNP','CHR','BP','ccgwas.snp')]
  her2.lumA.snp.BC.all <- rbind(her2.lumA.snp.BC.all, snp)
}

dim(her2.lumA.snp.BC.all)#1
her2.lumA.snp.BC.all <- as.data.frame(her2.lumA.snp.BC.all)
her2.lumA.snp.BC.undup <- filter(her2.lumA.snp.BC.all, !duplicated(her2.lumA.snp.BC.all$SNP)) #1
dim(her2.lumA.snp.BC.undup)
her2.lumA.sum <- fread ("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/BCAC_her2_lumA_gcta_input0428.ma")
her2.lumA.snp.BC.m <- merge(her2.lumA.snp.BC.undup,her2.lumA.sum,by.x="SNP", by.y="SNP") # examine if the snps can be found in the summary statistics
dim(her2.lumA.snp.BC.m) #1

her2.lumA.nomerge <- filter(her2.lumA.snp.BC.undup, !(her2.lumA.snp.BC.undup$SNP %in% her2.lumA.snp.BC.m$SNP)) #0
write.csv(her2.lumA.snp.BC.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/her2.lumA.snp.overall.0605.csv",row.names = F)
write.csv(her2.lumA.nomerge,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/her2.lumA.snp.overall.nofind.0605.csv",row.names = F)


snp.list.all <- c()
for(i in 1:nrow(her2.lumA.snp.list.ccgwas)){
  snp.list <- c()
  snp <- subset(her2.lumA.snp.BC.all, CHR==her2.lumA.snp.list.ccgwas$CHR[i] & as.numeric(BP) < her2.lumA.snp.list.ccgwas$BP[i]+500000 & as.numeric(BP) > her2.lumA.snp.list.ccgwas$BP[i]-500000)
  snp <- filter(snp, !duplicated(snp$SNP))
  colnames(snp)[1] <- c("SNP")
  snp.list <- rbind(her2.lumA.snp.list.ccgwas[i,c(1:3)],snp[,c(1:3)])
  snp.list$no.test <- paste0("test",i)
  snp.list$ccgwas.snp <- her2.lumA.snp.list.ccgwas[i,1]
  write.table(snp.list[,c(1:3)],paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/her2_lumA/test",i,".txt",sep=""),row.names=F,col.names = F,quote = F)
  snp.list.all <- rbind(snp.list.all,snp.list)
  
}
write.table(snp.list.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/gcta/test_list_subtypes0605/her2_lumA/test.all.txt")


