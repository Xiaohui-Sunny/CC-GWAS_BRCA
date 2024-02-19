#ldsc analysis
cd ldsc
source activate ldsc
conda activate ldsc
# conda deactive  


# bugs
##environment
# https://github.com/bulik/ldsc/issues/112
# https://github.com/bulik/ldsc/blob/master/environment.yml
# https://github.com/perslab/CELLECT/issues/36

## long time in munge
# https://github.com/bulik/ldsc/issues/281 add the code --chunksize 500000

###########################################################################
# 2020.04.20 
# If I used the previous code, there is an error in the python ldsc.py
# Error:
#ValueError: Improperly formatted sumstats file: ('No columns to parse from file',)
# To fix this bug, panda was updated to 0.20 (pip install pandas==0.20). However, munge_sumstats.py will take a long time to generate the summary data. Thus I added --chunksize 500000.

# environment.yml(04.20.2022)

(ldsc) sunx3@omega5:~/software/ldsc$ conda list
# packages in environment at /home/nfs/sunx3/anaconda3/envs/ldsc:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main
_openmp_mutex             4.5                       1_gnu
argparse                  1.3.0                    pypi_0    pypi
bitarray                  0.8.1                    pypi_0    pypi
ca-certificates           2022.3.29            h06a4308_0
certifi                   2020.6.20          pyhd3eb1b0_3
libedit                   3.1                  heed3624_0
libffi                    3.2.1             hf484d3e_1007
libgcc-ng                 9.3.0               h5101ec6_17
libgomp                   9.3.0               h5101ec6_17
libstdcxx-ng              9.3.0               hd4cf53a_17
ncurses                   6.0                  h9df7e31_2
nose                      1.3.0                    pypi_0    pypi
numpy                     1.16.0                   pypi_0    pypi
openssl                   1.0.2u               h7b6447c_0
pandas                    0.20.0                   pypi_0    pypi    # This is upated
pip                       19.3.1                   py27_0
pybedtools                0.7.0                    pypi_0    pypi
pysam                     0.19.0                   pypi_0    pypi
python                    2.7.13              heccc3f1_16
python-dateutil           2.8.2                    pypi_0    pypi
pytz                      2022.1                   pypi_0    pypi
readline                  7.0                  ha6073c6_4
scipy                     0.18.1                   pypi_0    pypi
setuptools                44.0.0                   py27_0
six                       1.16.0                   pypi_0    pypi
sqlite                    3.23.1               he433501_0
tk                        8.6.11               h1ccaba5_0
wheel                     0.37.1             pyhd3eb1b0_0
zlib                      1.2.12               h7f8727e_1




##munge data

## according to the orginal BC subtype GWAS paper, the total sample size is 106278 cases + 91477 controls= 197775
# change the sample size based on Xiang's code
# the default is to filter the SNP with info >0.9, if change the r2 >0.3 need to change the python.py
# BCAC_triple
python munge_sumstats.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Trip_ldsc_input_hapmap3_filtered.txt \
   --N 197775 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.3 and rerun (04.20)
python munge_sumstats_r2_0.3.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Trip_ldsc_input_hapmap3_filtered2.txt \
   --N 197775 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0419_2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.8, effective sample size and rerun(04.28)
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Trip_ldsc_input_neff.txt \
   --N 107976 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
# 91477 + 16499 =107976


# BCAC_LuminaB
python munge_sumstats.py \
  --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_ldsc_input_hapmap3_filtered2.txt \
  --N 197775 \
  --chunksize 500000 \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0419_2 \
  --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
## change the info to >0.3 and rerun (04.20)
python munge_sumstats_r2_0.3.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_ldsc_input_hapmap3_filtered2.txt \
   --N 197775 \
   --chunksize 98136 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0419_2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.8, effective sample size and rerun(04.28)
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_ldsc_input_neff.txt \
   --N 103891 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
# 91477 + 12414 =103891

# BCAC_LuminaA
python munge_sumstats.py \
  --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaA_ldsc_input.txt \
  --N 197775 \
  --chunksize 500000 \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA \
  --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
## change the info to >0.3 and rerun (04.20)
python munge_sumstats_r2_0.3.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaA_ldsc_input_hapmap3_filtered2.txt \
   --N 197775 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0419_2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.8, effective sample size and rerun(04.28)
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaA_ldsc_input_neff.txt \
   --N 152192 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
# 91477 + 60715 =152192



#BCAC_Her2
python munge_sumstats.py \
  --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Her2_enrich_ldsc_input.txt \
  --N 197775 \
  --chunksize 500000 \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20419 \
  --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
## change the info to >0.3 and rerun (04.20)
python munge_sumstats_r2_0.3.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Her2_enrich_ldsc_input_hapmap3_filtered2.txt \
   --N 197775 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20419_2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.8, effective sample size and rerun(04.28)
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Her2_enrich_ldsc_input_neff.txt \
   --N 97336 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
# 91477 + 5859 =152192


#BCAC_LumB_her2
python munge_sumstats.py \
  --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_Her2_ldsc_input_hapmap3_filtered2.txt \
  --N 197775 \
  --chunksize 500000 \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20419_2 \
  --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
## change the info to >0.3 and rerun (04.20)
python munge_sumstats_r2_0.3.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_Her2_ldsc_input_hapmap3_filtered2.txt \
   --N 197775 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20419_2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

## change the info to >0.8, effective sample size and rerun(04.28)
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_Her2_ldsc_input_neff.txt \
   --N 108420 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist
# 91477 + 16943 =108420

##ld regression
#LumA and LumB
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_0428

#LumA and Her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Her2_0428

#LumA and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB-her2_0428

#LumA and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Triple_0428

#LumB and Her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_Her2_0428

#LumB and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_LumB-her2_0428

#LumB and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_Triple_0428

#Her2 and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2_LumB-her2_0428

#Her2 and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2_Triple_0428

#LumB_her2 and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB-her2_Triple_0428

## cross-trait ldsc, convert to liability
#heribility: continuous trait
#liability: binatary trait
#For breast cancer, we need to convert heribility obatiend from ldsc to liability 
#LumA and LumB
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_liability \
  --pop-prev 0.00327,0.0007 \
  --samp-prev 0.275,0.04
  # --pop-prev for lumA is 0.006403*0.73*0.7=0.003271933 (population prevalence of BC in USA as a representative figure= 0.006403 (IARC), 73% of ER+ and/or PR+, Her2- BC, assume 70% are grade 1 and 2)
  # --pop-prev for lumB is 0.006403*0.11=0.00070433 (11% of ER+ and/or PR+, Her2+ BC) # check ACS Breast Cancer Facts & Figures 2019-2020
  # --samp-prev for lumA is jk= 45436/(45436+6659+9336+3026+9067+91477)=0.275368
  # --samp-prev for lumB is median(Neff_lumB)/sum(Neff_subtypes,N0)= 6659/(45436+6659+9336+3026+9067+91477)=0.04035733

## change effective sample size and rerun(0428)
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_liability \
  --pop-prev 0.00327,0.0007 \
  --samp-prev 0.2977583,0.0608807
  # --pop-prev for lumA is 0.006403*0.73*0.7=0.003271933 (population prevalence of BC in USA as a representative figure= 0.006403 (IARC), 73% of ER+ and/or PR+, Her2- BC, assume 70% are grade 1 and 2)
  # --pop-prev for lumB is 0.006403*0.11=0.00070433 (11% of ER+ and/or PR+, Her2+ BC) # check ACS Breast Cancer Facts & Figures 2019-2020
  # --samp-prev for lumA is jk= 60715/(60715+12414+16943+5859+16499+91477)=0.2977583
  # --samp-prev for lumB is median(Neff_lumB)/sum(Neff_subtypes,N0)= 12414/(60715+12414+16943+5859+16499+91477)=0.0608807

#LumA and LumB_her2
 python ldsc.py \
   --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0419_2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20419_2.sumstats.gz \
   --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_her2_liability \
   --pop-prev 0.00327,0.0014 \
   --samp-prev 0.275,0.057
  # --pop-prev for lumB_her2_neg is 0.006403*0.73*0.3=0.001402257 (73% of ER+ and/or PR+, Her2- BC, assume 30% are grade 3)
  # --samp-prev for lumB_her2_neg is median(Neff_lumB_her2)/sum(Neff_subtypes,N0)= 9336/(45436+6659+9336+3026+9067+91477)=0.05658148

## change effective sample size and rerun(0428)
 python ldsc.py \
   --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her20428.sumstats.gz \
   --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_her2_liability \
   --pop-prev 0.00327,0.0014 \
   --samp-prev 0.2977583,0.0830918
  # --pop-prev for lumB_her2_neg is 0.006403*0.73*0.3=0.001402257 (73% of ER+ and/or PR+, Her2- BC, assume 30% are grade 3)
  # --samp-prev for lumB_her2_neg is median(Neff_lumB_her2)/sum(Neff_subtypes,N0)= 16943/(60715+12414+16943+5859+16499+91477)=0.0830918


#LumA and her2
 python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0419_2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20419_2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_her2_liability \
  --pop-prev 0.00327,0.00026 \
  --samp-prev 0.275,0.018
  # --pop-prev for her2_enrich is 0.006403*0.04=0.00025612 (4% of ER- and PR-, Her2+ BC)
  # --samp-prev for her2_enrich is median(Neff_her2_enrich)/sum(Neff_subtypes,N0)= 3026/(45436+6659+9336+3026+9067+91477)= 0.01833928
 
## change effective sample size and rerun(0428)
 python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her20428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_her2_liability \
  --pop-prev 0.00327,0.00026 \
  --samp-prev 0.2977583,0.02873369
 # --pop-prev for her2_enrich is 0.006403*0.04=0.00025612 (4% of ER- and PR-, Her2+ BC)
  # --samp-prev for her2_enrich is median(Neff_her2_enrich)/sum(Neff_subtypes,N0)= 5859/(60715+12414+16943+5859+16499+91477)= 0.02873369
 
#LumA and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0419_2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0419_2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Trip_liability \
  --pop-prev 0.00327,0.00077 \
  --samp-prev 0.275,0.055
  # --pop-prev for Triple_neg is 0.006403*0.12=0.00076836 (12% of ER- and PR-, Her2- BC)
  # --samp-prev for Triple_neg is median(Neff_Triple_neg)/sum(Neff_subtypes,N0)=16499/(45436+6659+9336+3026+9067+91477)= 0.05495118

## change effective sample size and rerun(0428)
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA0428.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip0428.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Trip_liability \
  --pop-prev 0.00327,0.00077 \
  --samp-prev 0.2977583,0.08091434
  # --pop-prev for Triple_neg is 0.006403*0.12=0.00076836 (12% of ER- and PR-, Her2- BC)
  # --samp-prev for Triple_neg is median(Neff_Triple_neg)/sum(Neff_subtypes,N0)=16499/(60715+12414+16943+5859+16499+91477)= 0.08091434






#LumB and LumB_her2_neg
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB_her2.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out LumB.2_LumB_her2.2_liability \
  --pop-prev 0.0007,0.0014 \
  --samp-prev 0.04,0.057

#LumB and Her2_enriched
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/Her2.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out LumB.2_her2.2_liability \
  --pop-prev 0.0007,0.00026 \
  --samp-prev 0.04,0.018

#LumB and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/Trip.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out LumB.2_Trip.2_liability \
  --pop-prev 0.0007,0.00077 \
  --samp-prev 0.04,0.055

#LumB_her2_neg and Her2
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB_her2.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/Her2.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out LumB_her2.2_Her2.2_liability \
  --pop-prev 0.0014,0.00026 \
  --samp-prev 0.057,0.018
  
#LumB2_her2_neg and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB_her2.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/Trip.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out LumB_her2.2_Trip.2_liability \
  --pop-prev 0.0014,0.00077 \
  --samp-prev 0.057,0.055

#Her2_enrich and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/bra_subtypes_ccgwas/data/Her2.2.sumstats.gz,/home/nfs/sunx3/bra_subtypes_ccgwas/data/Trip.2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out Her2.2_Trip.2_liability \
  --pop-prev 0.00026,0.00077 \
  --samp-prev 0.018,0.055






























##Partitioned heritability
#Lumb_her2
python ldsc.py \
	--h2 /home/nfs/sunx3/bra_subtypes_ccgwas/data/LumB_her2.sumstats.gz \
	--ref-ld-chr /home/nfs/sunx3/software/ldsc/baseline/baseline. \
	--w-ld-chr /home/nfs/sunx3/software/ldsc/weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr /home/nfs/sunx3/software/ldsc/1000G_frq/1000G.mac5eur. \
  --out LumB-her2-baseline

