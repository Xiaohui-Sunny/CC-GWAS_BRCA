# Fine-mapping by using susie

## debug
#https://github.com/stephenslab/susieR/issues/148

library(data.table)
library(tidyverse)
library(susieR)

## prepare the summary statistics for the lead SNPs(n=100) from the ccgwas summary statistics
#### chunk the data into genomic risk loci, where all variants within 0.5Mb either side the lead SNP are included in the chunk
#### Summary statistics of genetic association studies typically contain effect size (β^ coefficient from regression), p-value and minor allele frequencies. 
setwd('/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas')

files <- c('lumB_lumA.out.results.gz',  'her2_enrich_lumA.out.results.gz',
           'trip_lumA.out.results.gz', 
           'trip_lumB.out.results.gz', 'trip_lumB_her2_neg.out.results.gz')


index.snp <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/BCAC_subtype_ccgwas_sugg_LD_prune_0605.csv",head=T)
index.snp <- as.data.frame(index.snp)

index.snp <- subset(index.snp, index.snp$index=='1')


### extract results for SNPs within +/-500 kb of index SNP
for (i in 1:length(files)){
  snp.all <- c()
  dat <- data.table::fread(files[i], header=T)
  dat <- data.frame(dat)
  index.sub <- index.snp[index.snp$cc==files[i],]
  for (l in 1:nrow(index.sub)){
    snp <- subset(dat, CHR==index.sub$CHR[l] & BP < index.sub$BP[l]+500000 & BP > index.sub$BP[l]-500000)
    snp$index <- ifelse(snp$CHR== index.sub$CHR[l] & snp$BP == index.sub$BP[l], 1, 0)
    snp$indexsnp <- index.sub$SNP[l]
    snp.all <- rbind(snp.all, snp)
  }
  write.csv(snp.all, paste0(files[i], 'ccgwas_1mb_index_snp_0605.csv'), row.names = F)
  write.table(snp.all$SNP,paste0(files[i], 'ccgwas_1mb_index_snp_list_0605.txt'), row.names = F,quote = F)
}



## fine-mapping
dbsnp <- fread("/home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur.bim", head=F)
dim(dbsnp)#22665064
dbsnp$chr_bp <- paste0(dbsnp$V1,"_", dbsnp$V4)

###################################################################################################
setwd("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie")
### triple and lumA
trip.lumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/trip_lumA.out.results.gzccgwas_1mb_index_snp_0605.csv", head=T)
trip.lumA$chr_bp <-paste0(trip.lumA$CHR, '_',  trip.lumA$BP)
dim(trip.lumA) #48116    15
#### merge with dbSNP (1000 Genomes)
trip.lumA <- merge(trip.lumA, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumA) # 47763     21
trip.lumA.snp.list <- trip.lumA$V2
write.table(trip.lumA.snp.list, "trip_lumA_ccgwas_1mb_index_snp_list_0605.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumA <- fread("Dosage.trip.lumA.raw", head=T) # obtain from plink
dosage.trip.lumA <- dosage.trip.lumA[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumA))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumA)`, "_", 2))


trip.lumA <- merge(trip.lumA,colnames,by.x="V2",by.y="V1",fill=T) 
#### check the allele
trip.lumA$beta <- ifelse(trip.lumA$EA==trip.lumA$V2.y, trip.lumA$OLS_beta,-trip.lumA$OLS_beta)
trip.lumA$z.change <- trip.lumA$beta/trip.lumA$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumA.index <- trip.lumA$indexsnp
trip.lumA.index <- unique(trip.lumA.index)#47
dosage.trip.lumA <- as.data.frame(dosage.trip.lumA)

trip.lumA.all <- c()
for (i in 19:length(trip.lumA.index)){
  trip.lumA.sub <-c()
  trip.lumA.re <-c()
  trip.lumA.sub <- as.data.frame(subset(trip.lumA, trip.lumA$indexsnp==trip.lumA.index[i]))
  v = as.vector(unlist(trip.lumA.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.trip.lumA.sub <- dosage.trip.lumA[,col_loc]
  ld.matrix.trip.lumA.sub <- cor(dosage.trip.lumA.sub)
  res <- susie_rss(trip.lumA.sub$z.change, as.matrix(ld.matrix.trip.lumA.sub), L=10,estimate_prior_variance = T)
  trip.lumA.re <- as.data.frame(summary(res)$cs)
  trip.lumA.re$cc.index <- trip.lumA.index[i]
  for (j in 1:nrow(trip.lumA.re)){
    row_loc <- trip.lumA.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    trip.lumA.cs.variable <- trip.lumA.sub[row_loc,1]
    trip.lumA.re[j,7] <- paste (trip.lumA.cs.variable,collapse = ",")

    }
  trip.lumA.all <- rbind(trip.lumA.all, trip.lumA.re)
  write.csv(trip.lumA.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumA0606.cs.",trip.lumA.index[i],".csv"),row.names=F)
  }

write.csv(trip.lumA.all, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumA0606.all.cs.csv")
#rs10516070:169512860:A:C(the 2th in the trip.lumA.index list) did not find the credible sets
#rs12435035:(13th in the trip.lumA.index list) did not find the credible sets
#rs1264478:102424918:G:A:(18th in the trip.lumA.index list) did not find the credible sets




###########################################################################################################
### trip and lumB
trip.lumB <- fread("trip_lumB.out.results.gzccgwas_1mb_index_snp_0605.csv", head=T)
trip.lumB$chr_bp <-paste0(trip.lumB$CHR, '_',  trip.lumB$BP)
dim(trip.lumB) #8696
#### merge with dbSNP (1000 Genomes)
trip.lumB <- merge(trip.lumB, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumB)  #8633

#### generage the SNP list
trip.lumB.snp.list <- trip.lumB$V2
write.table(trip.lumB.snp.list, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumB_ccgwas_1mb_index_snp_list_0605.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumB <- fread("Dosage.trip.lumB.raw", head=T) # obtain from plink
dosage.trip.lumB <- dosage.trip.lumB[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumB))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumB)`, "_", 2))


trip.lumB <- merge(trip.lumB,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
trip.lumB$beta <- ifelse(trip.lumB$EA==trip.lumB$V2.y, trip.lumB$OLS_beta,-trip.lumB$OLS_beta)
trip.lumB$z.change <- trip.lumB$beta/trip.lumB$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumB.index <- trip.lumB$indexsnp
trip.lumB.index <- unique(trip.lumB.index) #9


dosage.trip.lumB <- as.data.frame(dosage.trip.lumB)
trip.lumB.all <- c()
for (i in 2:length(trip.lumB.index)){
  trip.lumB.sub <-c()
  trip.lumB.sub <- as.data.frame(subset(trip.lumB, trip.lumB$indexsnp==trip.lumB.index[i]))
  v = as.vector(unlist(trip.lumB.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.trip.lumB.sub <- dosage.trip.lumB[,col_loc]
  ld.matrix.trip.lumB.sub <- cor(dosage.trip.lumB.sub)
  res <- susie_rss(trip.lumB.sub$z.change, as.matrix(ld.matrix.trip.lumB.sub), L=10, estimate_prior_variance = T)
  trip.lumB.re <- as.data.frame(summary(res)$cs)
  for (j in 1:nrow(trip.lumB.re)){
    row_loc <- trip.lumB.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    trip.lumB.cs.variable <- trip.lumB.sub[row_loc,1]
    trip.lumB.re[j,6] <- paste (trip.lumB.cs.variable,collapse = ",")
  }
  trip.lumB.re$cc.index <- trip.lumB.index[i]
  trip.lumB.all <- rbind(trip.lumB.all, trip.lumB.re)
  write.csv(trip.lumB.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB0606.cs.",trip.lumB.index[i],".csv"),row.names=F)
}
write.csv(trip.lumB.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB0606.cs.all.csv")
#rs7589172( 1st SNP in the trip.lumB.index) find no credible sets



##################################################################################################
### trip and lumB-her2-negative
trip.lumB.her2.neg <- fread("trip_lumB_her2_neg.out.results.gzccgwas_1mb_index_snp_0503.csv", head=T)
trip.lumB.her2.neg$chr_bp <-paste0(trip.lumB.her2.neg$CHR, '_',  trip.lumB.her2.neg$BP)
dim(trip.lumB.her2.neg) #20158    15
#### merge with dbSNP (1000 Genomes)
trip.lumB.her2.neg <- merge(trip.lumB.her2.neg, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumB.her2.neg)  #20000    21
trip.lumB.her2.neg.snp.list <- trip.lumB.her2.neg$V2
write.table(trip.lumB.her2.neg.snp.list, "trip_lumB_her2_neg_ccgwas_1mb_index_snp_list_0503.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumB.her2.neg <- fread("Dosage.trip.lumB.her2.neg.raw", head=T) # obtain from plink
dosage.trip.lumB.her2.neg <- dosage.trip.lumB.her2.neg[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumB.her2.neg))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumB.her2.neg)`, "_", 2))

trip.lumB.her2.neg <- merge(trip.lumB.her2.neg,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
trip.lumB.her2.neg$beta <- ifelse(trip.lumB.her2.neg$EA==trip.lumB.her2.neg$V2.y, trip.lumB.her2.neg$OLS_beta,-trip.lumB.her2.neg$OLS_beta)
trip.lumB.her2.neg$z.change <- trip.lumB.her2.neg$beta/trip.lumB.her2.neg$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumB.her2.neg.index <- trip.lumB.her2.neg$indexsnp
trip.lumB.her2.neg.index <- unique(trip.lumB.her2.neg.index)#9

dosage.trip.lumB.her2.neg <- as.data.frame(dosage.trip.lumB.her2.neg)
trip.lumB.her2.neg.all <- c()
for (i in 5:length(trip.lumB.her2.neg.index)){
  trip.lumB.her2.neg.sub <-c()
  trip.lumB.her2.neg.re <-c()
  trip.lumB.her2.neg.sub <- as.data.frame(subset(trip.lumB.her2.neg, trip.lumB.her2.neg$indexsnp==trip.lumB.her2.neg.index[i]))
  v = as.vector(unlist(trip.lumB.her2.neg.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.trip.lumB.her2.neg.sub <- dosage.trip.lumB.her2.neg[,col_loc]
  ld.matrix.trip.lumB.her2.neg.sub <- cor(dosage.trip.lumB.her2.neg.sub)
  res <- susie_rss(trip.lumB.her2.neg.sub$z.change, as.matrix(ld.matrix.trip.lumB.her2.neg.sub), L=10,estimate_prior_variance = T)
  trip.lumB.her2.neg.re <- summary(res)$cs
   for (j in 1:nrow(trip.lumB.her2.neg.re)){
     row_loc <- trip.lumB.her2.neg.re[j,5]
     row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
     trip.lumB.her2.neg.cs.variable <- trip.lumB.her2.neg.sub[row_loc,1]
     trip.lumB.her2.neg.re[j,6] <- paste (trip.lumB.her2.neg.cs.variable,collapse = ",")
   }
  trip.lumB.her2.neg.re$cc.index <- trip.lumB.her2.neg.index[i]
  trip.lumB.her2.neg.all <- rbind(trip.lumB.her2.neg.all, trip.lumB.her2.neg.re)
  write.csv(trip.lumB.her2.neg.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB.her2.neg0503.cs.",trip.lumB.her2.neg.index[i],".csv"),row.names=F)
}

#rs1243184 (4th in the list) did not find the crediable sets
write.csv(trip.lumB.her2.neg.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB.her2.neg0503.cs.all.csv")
##################################################################################################
### her2 and lumA
her2.lumA <- fread("her2_enrich_lumA.out.results.gzccgwas_1mb_index_snp_0503.csv", head=T)
her2.lumA$chr_bp <-paste0(her2.lumA$CHR, '_',  her2.lumA$BP)
dim(her2.lumA) #5773   15
#### merge with dbSNP (1000 Genomes)
her2.lumA <- merge(her2.lumA, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(her2.lumA) #217816     21
her2.lumA.snp.list <- her2.lumA$V2
write.table(her2.lumA.snp.list, "her2_enrich_lumA_ccgwas_1mb_index_snp_list_0503.txt",row.names=T, col.names = T,quote=F)


dosage.her2.lumA <- fread("Dosage.her2.lumA.raw", head=T) # obtain from plink
dosage.her2.lumA <- dosage.her2.lumA[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.her2.lumA))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.her2.lumA)`, "_", 2))

her2.lumA <- merge(her2.lumA,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
her2.lumA$beta <- ifelse(her2.lumA$EA==her2.lumA$V2.y, her2.lumA$OLS_beta,-her2.lumA$OLS_beta)
her2.lumA$z.change <- her2.lumA$beta/her2.lumA$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
her2.lumA.index <- her2.lumA$indexsnp
her2.lumA.index <- unique(her2.lumA.index) #3


dosage.her2.lumA <- as.data.frame(dosage.her2.lumA)
her2.lumA.all <- c()
for (i in 1:length(her2.lumA.index)){
  her2.lumA.sub <-c()
  her2.lumA.sub <- as.data.frame(subset(her2.lumA, her2.lumA$indexsnp==her2.lumA.index[i]))
  v = as.vector(unlist(her2.lumA.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.her2.lumA.sub <- dosage.her2.lumA[,col_loc]
  ld.matrix.her2.lumA.sub <- cor(dosage.her2.lumA.sub)
  res <- susie_rss(her2.lumA.sub$z.change, as.matrix(ld.matrix.her2.lumA.sub), L=10,estimate_prior_variance = T)
  her2.lumA.re <- summary(res)$cs
  for (j in 1:nrow(her2.lumA.re)){
    row_loc <- her2.lumA.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    her2.lumA.cs.variable <- her2.lumA.sub[row_loc,1]
    her2.lumA.re[j,6] <- paste (her2.lumA.cs.variable,collapse = ",")
  }
  her2.lumA.re$cc.index <- her2.lumA.index[i]
  her2.lumA.all <- rbind(her2.lumA.all, her2.lumA.re)
  
  write.csv(her2.lumA.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/her2.lumA0503.cs.",her2.lumA.index[i],".csv"),row.names=F)
}
write.csv(her2.lumA.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/her2.lumA0503.cs.all.csv")

####################################################################################################################
### her2 and lumB
her2.lumB <- fread("her2_enrich_lumB.out.results.gzccgwas_1mb_index_snp_0503.csv", head=T)
her2.lumB$chr_bp <-paste0(her2.lumB$CHR, '_',  her2.lumB$BP)
dim(her2.lumB) #4512   15
#### merge with dbSNP (1000 Genomes)
her2.lumB <- merge(her2.lumB, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(her2.lumB) #4478   21
her2.lumB.snp.list <- her2.lumB$V2
write.table(her2.lumB.snp.list, "her2_enrich_lumB_ccgwas_1mb_index_snp_list_0503.txt",row.names=T, col.names = T,quote=F)

dosage.her2.lumB <- fread("Dosage.her2.lumB.raw", head=T) # obtain from plink
dosage.her2.lumB <- dosage.her2.lumB[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.her2.lumB))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.her2.lumB)`, "_", 2))

her2.lumB <- merge(her2.lumB,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
her2.lumB$beta <- ifelse(her2.lumB$EA==her2.lumB$V2.y, her2.lumB$OLS_beta,-her2.lumB$OLS_beta)
her2.lumB$z.change <- her2.lumB$beta/her2.lumB$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
her2.lumB.index <- her2.lumB$indexsnp
her2.lumB.index <- unique(her2.lumB.index) #2


dosage.her2.lumB <- as.data.frame(dosage.her2.lumB)
her2.lumB.all <- c()
for (i in 1:length(her2.lumB.index)){
  her2.lumB.sub <-c()
  her2.lumB.sub <- subset(her2.lumB, her2.lumB$indexsnp==her2.lumB.index[i])
  v = as.vector(unlist(her2.lumB.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.her2.lumB.sub <- dosage.her2.lumB[,col_loc]
  ld.matrix.her2.lumB.sub <- cor(dosage.her2.lumB.sub)
  res <- susie_rss(her2.lumB.sub$z.change, as.matrix(ld.matrix.her2.lumB.sub), L=10,estimate_prior_variance = T)
  her2.lumB.re <- summary(res)$cs
  for (j in 1:nrow(her2.lumB.re)){
    row_loc <- her2.lumB.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    her2.lumB.cs.variable <- her2.lumB.sub[row_loc,1]
    her2.lumB.re[j,6] <- paste (her2.lumB.cs.variable,collapse = ",")
  }
  her2.lumB.re$cc.index <- her2.lumB.index[i]
  her2.lumB.all <- rbind(her2.lumB.all, her2.lumB.re)
  
  write.csv(her2.lumB.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/her2.lumB0503.cs.",her2.lumB.index[i],".csv"),row.names=F)
  }
 
write.csv(her2.lumB.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/her2.lumB0503.cs.all.csv")




###################################################################################################################
### lumB and lumA
lumB.lumA <- fread("lumB_lumA.out.results.gzccgwas_1mb_index_snp_0503.csv", head=T)
lumB.lumA$chr_bp <-paste0(lumB.lumA$CHR, '_',  lumB.lumA$BP)
dim(lumB.lumA) #2039   15
#### merge with dbSNP (common)
lumB.lumA <- merge(lumB.lumA, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(lumB.lumA) #2017   21

####prepare index SNP for LD matrix
lumB.lumA.snp.list <- as.data.frame(lumB.lumA$V2)
write.table(lumB.lumA.snp.list, "lumB_lumA_ccgwas_1mb_index_snp_list_0503.txt",row.names=T, col.names = T,quote=F)

#### input the LD matrix from Plink
#ld.matrix.lumB.lumA <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/LDmatrix.lumB.lumA.ld", head=F)

#### perform fine-mapping
#res <- susie_rss(lumB.lumA$z.change, as.matrix(ld.matrix.lumB.lumA), L=10)
#Error:
#For large R or large XtX, consider installing the  Rfast package for better performance.
#XtX is not symmetric; forcing XtX to be symmetric by  replacing XtX with (XtX + t(XtX))/2
#Error in susie_suff_stat(XtX = R, Xty = z, n = 2, yty = 1, scaled_prior_variance = prior_variance,  : 
#                           The estimated prior variance is unreasonably large.
#                         Please check the input.

#Debug: change the allele to 1000 genomes
#### input the dosage of SNPs and than calculate the LD matrix by using R cor() function
dosage.lumB.lumA <- fread("Dosage.lumB.lumA.raw", head=T) # obtain from plink
dosage.lumB.lumA <- dosage.lumB.lumA[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.lumB.lumA))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.lumB.lumA)`, "_", 2))
lumB.lumA <- merge(lumB.lumA,colnames,by.x="V2",by.y="V1", fill=T)
#### check the allele
lumB.lumA$beta <- ifelse(lumB.lumA$EA==lumB.lumA$V2.y, lumB.lumA$OLS_beta,-lumB.lumA$OLS_beta)
lumB.lumA$z.change <- lumB.lumA$beta/lumB.lumA$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
lumB.lumA.index <- lumB.lumA$indexsnp
lumB.lumA.index <- unique(lumB.lumA.index) #1
dosage.lumB.lumA <- as.data.frame(dosage.lumB.lumA)

for (i in 1:length(lumB.lumA.index)){
  lumB.lumA.sub <-c()
  lumB.lumA.sub <- subset(lumB.lumA, lumB.lumA$indexsnp==lumB.lumA.index[1])
  v = as.vector(unlist(lumB.lumA.sub[,1]))
  rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
  col_loc <- rowcol[,1]
  dosage.lumB.lumA.sub <- dosage.lumB.lumA[,col_loc]
  ld.matrix.lumB.lumA.sub <- cor(dosage.lumB.lumA.sub)
  res <- susie_rss(lumB.lumA.sub$z.change, as.matrix(ld.matrix.lumB.lumA.sub), L=10,estimate_prior_variance = T)
  lumB.lumA.re <- summary(res)$cs
  for (j in 1:nrow(lumB.lumA.re)){
    row_loc <- lumB.lumA.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    lumB.lumA.cs.variable <- lumB.lumA.sub[row_loc,1]
    lumB.lumA.re[j,6] <- paste (lumB.lumA.cs.variable,collapse = ",")
  }
  write.csv(lumB.lumA.re,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/lumB.lumA0503.cs.",lumB.lumA.index[1],".csv"),row.names=F)
}


#res <- susie_rss(lumB.lumA$z.change, as.matrix(ld.matrix.lumB.lumA), L=10,estimate_prior_variance = T)
#Error:
#For large R or large XtX, consider installing the  Rfast package for better performance.
#XtX is not symmetric; forcing XtX to be symmetric by  replacing XtX with (XtX + t(XtX))/2
#Error in susie_suff_stat(XtX = R, Xty = z, n = 2, yty = 1, scaled_prior_variance = prior_variance,  : 
#                           The estimated prior variance is unreasonably large.
#                         Please check the input.

## change the allele to 1000 Genomes, it works
res <- susie_rss(lumB.lumA$z.change, as.matrix(ld.matrix.lumB.lumA), L=10,estimate_prior_variance = T)
lumB.lumA.re <- summary(res)$cs
write.csv (lumB.lumA.re, "/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/lumB.lumA.cs.0310.csv",row.names = F)


##################################################################################
# combine the SNP within 500kb of the index SNP
## Triple and LumB-her2-neg
### trip.lumB.cs.rs2912779 and trip.lumB.cs.rs1649154:123360913:T:C
trip.lumB <- fread("trip_lumB.out.results.gzccgwas_1mb_index_snp_0309.csv", head=T)
trip.lumB$chr_bp <-paste0(trip.lumB$CHR, '_',  trip.lumB$BP)
dim(trip.lumB) #26200
#### merge with dbSNP (1000 Genomes)
trip.lumB <- merge(trip.lumB, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumB)  #25983

#### generage the SNP list
trip.lumB.snp.list <- trip.lumB$V2
write.table(trip.lumB.snp.list, "trip_lumB_ccgwas_1mb_index_snp_list_0309.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumB <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/Dosage.trip.lumB.raw", head=T) # obtain from plink
dosage.trip.lumB <- dosage.trip.lumB[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumB))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumB)`, "_", 2))

trip.lumB <- merge(trip.lumB,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
trip.lumB$beta <- ifelse(trip.lumB$EA==trip.lumB$V2.y, trip.lumB$OLS_beta,-trip.lumB$OLS_beta)
trip.lumB$z.change <- trip.lumB$beta/trip.lumB$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumB.index <- trip.lumB$indexsnp
trip.lumB.index <- unique(trip.lumB.index) #10

dosage.trip.lumB <- as.data.frame(dosage.trip.lumB)
trip.lumB.sub <- as.data.frame(subset(trip.lumB, trip.lumB$indexsnp=='rs2912779' | trip.lumB$indexsnp=="rs1649154:123360913:T:C"))

trip.lumB.sub <- filter(trip.lumB.sub,!duplicated(trip.lumB.sub$SNP)) # remove the duplicated SNPs
v = as.vector(unlist(trip.lumB.sub[,1]))
rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
col_loc <- rowcol[,1]
dosage.trip.lumB.sub <- dosage.trip.lumB[,col_loc]
ld.matrix.trip.lumB.sub <- cor(dosage.trip.lumB.sub)
res <- susie_rss(trip.lumB.sub$z.change, as.matrix(ld.matrix.trip.lumB.sub), L=10,estimate_prior_variance = T)
trip.lumB.re <- as.data.frame(summary(res)$cs)

for (j in 1:nrow(trip.lumB.re)){
   row_loc <- trip.lumB.re[j,5]
   row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
   trip.lumB.cs.variable <- trip.lumB.sub[row_loc,1]
   trip.lumB.re[j,6] <- paste (trip.lumB.cs.variable,collapse = ",")
}

write.csv(trip.lumB.re,paste0("/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/trip.lumB.cs.","rs2912779 and rs1649154:123360913:T:C",".csv"),row.names=F)

#########################################################
## Triple and LumA

trip.lumA <- fread("trip_lumA.out.results.gzccgwas_1mb_index_snp_0309.csv", head=T)
trip.lumA$chr_bp <-paste0(trip.lumA$CHR, '_',  trip.lumA$BP)
dim(trip.lumA) #166211 15
#### merge with dbSNP (1000 Genomes)
trip.lumA <- merge(trip.lumA, dbsnp, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumA) #164937 18
trip.lumA.snp.list <- trip.lumA$V2

dosage.trip.lumA <- fread("/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/Dosage.trip.lumA.raw", head=T) # obtain from plink
dosage.trip.lumA <- dosage.trip.lumA[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumA))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumA)`, "_", 2))

trip.lumA <- merge(trip.lumA,colnames,by.x="V2",by.y="V1",fill=T) 
#### check the allele
trip.lumA$beta <- ifelse(trip.lumA$EA==trip.lumA$V2.y, trip.lumA$OLS_beta,-trip.lumA$OLS_beta)
trip.lumA$z.change <- trip.lumA$beta/trip.lumA$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
dosage.trip.lumA <- as.data.frame(dosage.trip.lumA)

### trip.lumA.cs.rs111458676:121287994:A:G and trip.lumA.cs.rs12123157:121237156:G:T
### rs78720385:17456392:C:T and rs61494113
### c19_pos17208709 and rs60960968:17342284:G:A
### rs147070975:69475386:C:T and rs12806763 and rs657315 and rs17136641:69328512:G:A and chr11_69281829_A_C and rs7927237 and rs139074399:69237593:A:G and rs72940459:69227074:C:T
### chr10_123391068_A_G and c10_pos123355714 and rs12268516 and rs45631630 and rs2912780 and rs2981429 and  rs72832370:123294692:C:T
### rs72709458:1283755:C:T and rs13167280
### rs12223381 and chr11_107965390_A_G
### rs10941679 and rs12520604:44476872:A:C and rs4866900:44445100:G:A
### rs2464264 and 12:115792584:C:T and 12:115723292:T:C and 12:115583816:G:A
### rs62355901:56053535:T:C and chr5_56045081_C_T
### rs2372966:217967719:T:C and chr2_217920769_G_T and rs72949860 and chr2_217912437_C_T
### chr4_175842786_A_G and chr4_175839725_C_T
### rs3803659:52586581:A:G and chr16_52599188_C_T
### rs4752537:123092921:A:G and rs35401699
trip.lumA.sub <- as.data.frame(subset(trip.lumA, trip.lumA$indexsnp=="rs4752537:123092921:A:G" | trip.lumA$indexsnp=="rs35401699"))
trip.lumA.sub <- filter(trip.lumA.sub,!duplicated(trip.lumA.sub$SNP))

v = as.vector(unlist(trip.lumA.sub[,1]))
rowcol <- do.call('rbind', lapply(v, function(x) which(colnames$V1==x,arr.ind = T)))
col_loc <- rowcol[,1]
dosage.trip.lumA.sub <- dosage.trip.lumA[,col_loc]
ld.matrix.trip.lumA.sub <- cor(dosage.trip.lumA.sub)
res <- susie_rss(trip.lumA.sub$z.change, as.matrix(ld.matrix.trip.lumA.sub), L=10,estimate_prior_variance = T)
trip.lumA.re <- as.data.frame(summary(res)$cs)
for (j in 1:nrow(trip.lumA.re)){
    row_loc <- trip.lumA.re[j,5]
    row_loc <- as.numeric(unlist(stringr::str_split(row_loc,",")))
    trip.lumA.cs.variable <- trip.lumA.sub[row_loc,1]
    trip.lumA.re[j,6] <- paste (trip.lumA.cs.variable,collapse = ",")
  }
write.csv(trip.lumA.re,paste0("/home/nfs/sunx3/bra_subtypes_ccgwas/result/susie/trip.lumA.cs.","rs4752537:123092921:A:G and rs35401699",".csv"),row.names=F)





