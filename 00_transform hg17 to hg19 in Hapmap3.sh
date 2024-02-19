#transform hg17 to hg19 in Hapmap3 by liftover

#introduction:liftover
#https://www.biostars.org/p/306616/
#https://genome.sph.umich.edu/wiki/LiftOver

#download liftover
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver

#download annotation:
#hg38Tohg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
#hg19Tohg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#hg18Tohg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz


# transform hapmap3 map, ped to bed
/home/nfs/sunx3/software/plink2.0/plink \
  --file /home/nfs/sunx3/data/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode \
  --make-bed \
  --out /home/nfs/sunx3/data/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode

# the input of liftover should be bed (chr, start, end(start+1)), and without colnames
awk '{print "chr"$1,"\t",$4,"\t",$4+1,"\t",$2}' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3.hg18.bed

# transform by using liftover
#$ chmod +x ./filePath  make liftover avaliable in linux
chmod +x /home/nfs/sunx3/software/liftover/liftOver 
# transform  hg18.bed to hg19.bed
/home/nfs/sunx3/software/liftover/liftOver \
  /home/nfs/sunx3/data/hapmap3/hapmap3.hg18.bed \
  /home/nfs/sunx3/software/liftover/hg18ToHg19.over.chain.gz \
  /home/nfs/sunx3/data/hapmap3/hapmap3.hg19.bed \
  unmap
# change the name hg19.bed to map (end(start))
awk '{print $1,"\t",$4,"\t",0,"\t",$2}' hapmap3.hg19.bed > hapmap3_r1_b37_fwd_consensus.qc.poly.recode.map

perl -p -i -e 's/chr//g' hapmap3_r1_b37_fwd_consensus.qc.poly.recode.map

mv hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped hapmap3_r1_b37_fwd_consensus.qc.poly.recode.ped
/home/nfs/sunx3/software/plink2.0/plink \
   --file hapmap3_r1_b37_fwd_consensus.qc.poly.recode \
   --make-bed \
   --out hapmap3_r1_b37_fwd_consensus.qc.poly.recode

