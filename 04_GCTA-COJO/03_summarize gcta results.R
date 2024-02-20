#summary for gcta(1 step)
library(tidyverse)

#her2 and lumA # 1

temp <- as.data.frame(NA)
files <- c()

for (i in 1:1){
  temp <- read.table (paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/her2_lumA/0605test",i,".jma.cojo",sep=""), head=T)
  temp$test <- paste0("test",i)
  files <- rbind(files, temp)
}
write.csv(files, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/gcta_result_her2_lumA_0605.csv", row.names = F)


#lumB and lumA #1
temp <- read.table (paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/lumB_lumA/0605test",1,".jma.cojo",sep=""), head=T)
temp$test <- paste0("test",1)
write.csv(temp, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/gcta_result_lumB_lumA_0605.csv", row.names = F)

#Triple and lumA #21
temp <- as.data.frame(NA)
files <- c()

for (i in 1:21){
  temp <- read.table (paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/trip_lumA/0605test",i,".jma.cojo",sep=""), head=T)
  temp$test <- paste0("test",i)
  files <- rbind(files, temp)
}

write.csv(files, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/gcta_result_trip_lumA_0605.csv", row.names = F)


#Triple and lumB #4
temp <- as.data.frame(NA)
files <- c()

for (i in 1:4){
  temp <- read.table (paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/trip_lumB/0605test",i,".jma.cojo",sep=""), head=T)
  temp$test <- paste0("test",i)
  files <- rbind(files, temp)
}
write.csv(files, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/gcta_result_trip_lumB_0605.csv", row.names = F)



#Triple and trip_lumB_her2_neg #2
temp <- as.data.frame(NA)
files <- c()

for (i in 1:2){
  temp <- read.table (paste("/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/trip_lumB_her2_neg/0605test",i,".jma.cojo",sep=""), head=T)
  temp$test <- paste0("test",i)
  files <- rbind(files, temp)
}

write.csv(files, "/home/nfs/sunx3/project/bra_subtypes_ccgwas/result/gcta/gcta_result_trip_lumB_her2_neg_0605.csv", row.names = F)

