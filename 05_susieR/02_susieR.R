# Fine-mapping by using susie
library(data.table)
library(tidyverse)
library(susieR)

index.snp <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/BCAC_subtype_ccgwas_sugg_LD_prune.csv",head=T)
index.snp <- as.data.frame(index.snp)
index.snp <- subset(index.snp, index.snp$index=='1')

files <- c('trip_lumA.out.results.gz','trip_lumB.out.results.gz', 'trip_lumB_her2_neg.out.results.gz')
### extract results for SNPs within +/-500 kb of index SNP
for (i in 1:length(files)){
  snp.all <- c()
  dat <- data.table::fread(paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/",files[i]), header=T)
  dat <- data.frame(dat)
  index.sub <- index.snp[index.snp$cc==files[i],]
  for (l in 1:nrow(index.sub)){
    snp <- subset(dat, CHR==index.sub$CHR[l] & BP < index.sub$BP[l]+500000 & BP > index.sub$BP[l]-500000)
    snp$index <- ifelse(snp$CHR== index.sub$CHR[l] & snp$BP == index.sub$BP[l], 1, 0)
    snp$indexsnp <- index.sub$SNP[l]
    snp.all <- rbind(snp.all, snp)
  }
  write.csv(snp.all, paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/",files[i], "ccgwas_1mb_index_snp.csv"), row.names = F)
  write.table(snp.all$SNP,paste0("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/",files[i], "ccgwas_1mb_index_snp_list.txt"), row.names = F,quote = F)
}



## fine-mapping
### triple and lumA
trip.lumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumA.out.results.gzccgwas_1mb_index_snp.csv", head=T)
trip.lumA$chr_bp <-paste0(trip.lumA$CHR, '_',  trip.lumA$BP)
dim(trip.lumA) 
#### merge with 1000 Genomes
geno.1k <- fread("/home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur.bim", head=F)
dim(geno.1k)#22665064
geno.1k$chr_bp <- paste0(geno.1k$V1,"_", geno.1k$V4)

trip.lumA <- merge(trip.lumA, geno.1k, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumA) 
trip.lumA.snp.list <- trip.lumA$V2
write.table(trip.lumA.snp.list, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumA_ccgwas_1mb_index_snp_list.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumA <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/Dosage.trip.lumA.raw", head=T) # obtain from plink
dosage.trip.lumA <- dosage.trip.lumA[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumA))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumA)`, "_", 2))


trip.lumA <- merge(trip.lumA,colnames,by.x="V2",by.y="V1",fill=T) 
#### check the allele
trip.lumA$beta <- ifelse(trip.lumA$EA==trip.lumA$V2.y, trip.lumA$OLS_beta,-trip.lumA$OLS_beta)
trip.lumA$z.change <- trip.lumA$beta/trip.lumA$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumA.index <- trip.lumA$indexsnp
trip.lumA.index <- unique(trip.lumA.index)
dosage.trip.lumA <- as.data.frame(dosage.trip.lumA)

trip.lumA.all <- c()
for (i in 1:length(trip.lumA.index)){
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
  trip.lumA.re$cc.index <- trip.lumA.index[i]
  trip.lumA.all <- rbind(trip.lumA.all, trip.lumA.re)
  
  }

write.csv(trip.lumA.all, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumA.cs.all.csv",row.names=F)


###########################################################################################################
### trip and lumB
trip.lumB <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumB.out.results.gzccgwas_1mb_index_snp.csv", head=T)
trip.lumB$chr_bp <-paste0(trip.lumB$CHR, '_',  trip.lumB$BP)
dim(trip.lumB) 
trip.lumB <- merge(trip.lumB, geno.1k, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumB) 

#### generate the SNP list
trip.lumB.snp.list <- trip.lumB$V2
write.table(trip.lumB.snp.list, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumB_ccgwas_1mb_index_snp_list.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumB <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/Dosage.trip.lumB.raw", head=T) # obtain from plink
dosage.trip.lumB <- dosage.trip.lumB[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumB))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumB)`, "_", 2))


trip.lumB <- merge(trip.lumB, colnames, by.x="V2",by.y="V1",fill=T)
#### check the allele
trip.lumB$beta <- ifelse(trip.lumB$EA==trip.lumB$V2.y, trip.lumB$OLS_beta,-trip.lumB$OLS_beta)
trip.lumB$z.change <- trip.lumB$beta/trip.lumB$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumB.index <- trip.lumB$indexsnp
trip.lumB.index <- unique(trip.lumB.index) 


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
}
write.csv(trip.lumB.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB.cs.all.csv")


##################################################################################################
### trip and lumB-her2-negative
trip.lumB.her2.neg <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumB_her2_neg.out.results.gzccgwas_1mb_index_snp.csv", head=T)
trip.lumB.her2.neg$chr_bp <-paste0(trip.lumB.her2.neg$CHR, '_',  trip.lumB.her2.neg$BP)
dim(trip.lumB.her2.neg) 

trip.lumB.her2.neg <- merge(trip.lumB.her2.neg, geno.1k, by.x="chr_bp", by.y="chr_bp")
dim(trip.lumB.her2.neg)  
trip.lumB.her2.neg.snp.list <- trip.lumB.her2.neg$V2
write.table(trip.lumB.her2.neg.snp.list, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/trip_lumB_her2_neg_ccgwas_1mb_index_snp_list.txt",row.names=T, col.names = T,quote=F)

dosage.trip.lumB.her2.neg <- fread("/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/susie/Dosage.trip.lumB.her2.neg.raw", head=T) # obtain from plink
dosage.trip.lumB.her2.neg <- dosage.trip.lumB.her2.neg[,-c(1:6)]
colnames <- as.data.frame(colnames(dosage.trip.lumB.her2.neg))
colnames <- as.data.frame(stringr::str_split_fixed(colnames$`colnames(dosage.trip.lumB.her2.neg)`, "_", 2))

trip.lumB.her2.neg <- merge(trip.lumB.her2.neg,colnames,by.x="V2",by.y="V1",fill=T)
#### check the allele
trip.lumB.her2.neg$beta <- ifelse(trip.lumB.her2.neg$EA==trip.lumB.her2.neg$V2.y, trip.lumB.her2.neg$OLS_beta,-trip.lumB.her2.neg$OLS_beta)
trip.lumB.her2.neg$z.change <- trip.lumB.her2.neg$beta/trip.lumB.her2.neg$OLS_se

##### generate the LD matrix for each index SNP and its nearby SNPs
trip.lumB.her2.neg.index <- trip.lumB.her2.neg$indexsnp
trip.lumB.her2.neg.index <- unique(trip.lumB.her2.neg.index)

dosage.trip.lumB.her2.neg <- as.data.frame(dosage.trip.lumB.her2.neg)
trip.lumB.her2.neg.all <- c()
for (i in 2:length(trip.lumB.her2.neg.index)){
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
}

write.csv(trip.lumB.her2.neg.all,"/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/susie/trip.lumB.her2.neg.cs.all.csv")
