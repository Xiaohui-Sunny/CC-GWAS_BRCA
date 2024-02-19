#Annotation for 100 index SNP from ccGWAS
cd software
cd annovar

#prepare input file
#.avinput: five columns: chromosome, start position, end position, reference nucleoties, the observed nucleoties. Save text file to .avinput
#example:
#14	29976101	29976101	G	T	rs79609746
#9	17882518	17882518	T	C	rs7873122
#9	115904820	115904820	A	G	rs9632905
#input file saved as "G:\My work\BRA subtypes_ccGWAS\results\tempo"

#Annotation analysis
/home/nfs/sunx3/software/annovar/annotate_variation.pl \
 -geneanno \
 -dbtype refGene \
 -buildver hg19 /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ccgwas/BCAC_subtype_ccgwas_sugg_index_0428.avinput humandb
