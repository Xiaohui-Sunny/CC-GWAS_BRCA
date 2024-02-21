#MAGMA gene-based anlaysis
## Introduction: 
### https://gwaslab.com/2021/12/13/magma/
### https://www.jianshu.com/p/15c478da414d
### https://www.cnblogs.com/jialinliu/p/13854577.html
### https://ctg.cncr.nl/software/magma
### Note: FUMA can also perform MAGMA analysis (never try)

## Analyze (whole-genome wide gene-based analysis)
### Step1: annotation: assign a SNP to a gene if the SNP's location falls inside the region of a gene(transcription start to the stop)
#### triple and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumA.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome0608

#### triple and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumB.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome0608

#### triple and lumB-her2-neg
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumB.her2.neg.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome0608

#### Her2 and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/her2.lumA.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome0608

#### Her2 and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/her2.lumB.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome0608

#### lumA and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/lumA.lumB.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome0608

#### Triple and her2
/home/nfs/sunx3/software/magma_v1.10/magma \
  --annotate \
  --snp-loc /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.her2.snp.loc0608.txt \
  --gene-loc /home/nfs/sunx3/data/genelocation_protein_gene/NCBI37.3/NCBI37.3.gene.loc \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome0608


### Step 2: gene-based analysis (summary data)
#### triple and lumA (N=16499+60715) case number of triple and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumA.magma0608.txt   N=77214 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome.gene.based0608

#### triple and lumB (N=16499+12414) case number of triple and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumB.magma0608.txt   N=28913 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome.gene.based0608

#### triple and lumB-her2-neg (N=16499+16943 )
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.lumB.her2.neg.magma0608.txt   N=33442 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome.gene.based0608

#### Her2 and lumA (N=5859+60715) 
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/her2.lumA.magma0608.txt   N=66574 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome.gene.based0608

#### Her2 and lumB (N=5859+12414)
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/her2.lumB.magma0608.txt   N=18273 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome.gene.based0608

#### lumA and lumB (N=60715+12414)
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/lumA.lumB.magma0608.txt   N=73129 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome.gene.based0608

#### triple and her2 (N=16499+5859)
/home/nfs/sunx3/software/magma_v1.10/magma \
  --bfile /home/nfs/sunx3/data/1000genomes/g1000_eur/g1000_eur \
  --pval /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/magma/trip.her2.magma0608.txt   N=22358 \
  --gene-annot /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome0608.genes.annot \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome.gene.based0608




### Step 3: gene-set anlyasis (use MSigDB Canonical all pathways)
#### Note: the output from magma used Entrez ID (NCBI gene ID) to represent the gene
#### Note: the full description for MSigDB pathway can be found https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=BP

#### all gene sets
#### triple and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome.gene.set

#### triple and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome.gene.set

#### triple and lumB-her2-neg
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome.gene.set

#### Her2 and lumA 
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome.gene.set

#### Her2 and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome.gene.set

####lumA and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome.gene.based.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/msigdb.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome.gene.set


#### GO gene sets
#### triple and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome0504.go.set

#### triple and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome0504.go.set

#### triple and lumB-her2-neg
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome0504.go.set

#### Her2 and lumA 
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome0504.go.set

#### Her2 and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome0504.go.set

####lumA and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome0504.go.set

####Triple and Her2
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c5.go.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome0504.go.set

### canonical pathways 
#### triple and lumA
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumA.wholegenome0504.curated.set

#### triple and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.wholegenome0504.curated.set

#### triple and lumB-her2-neg
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.lumB.her2.neg.wholegenome0504.curated.set

#### Her2 and lumA 
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumA.wholegenome0504.curated.set

#### Her2 and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/her2.lumB.wholegenome0504.curated.set

####lumA and lumB
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/lumA.lumB.wholegenome0504.curated.set

####Triple and Her2
/home/nfs/sunx3/software/magma_v1.10/magma \
  --gene-results /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome.gene.based0504.genes.raw \
  --set-annot /home/nfs/sunx3/data/msigdb/c2.all.v7.5.1.entrez.gmt \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/magma/trip.her2.wholegenome0504.curated.set