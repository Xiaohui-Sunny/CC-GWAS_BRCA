# ldsc analysis
cd ldsc
source activate ldsc
conda activate ldsc
# conda deactive  

##environment
# https://github.com/bulik/ldsc/issues/112
# https://github.com/bulik/ldsc/blob/master/environment.yml
# https://github.com/perslab/CELLECT/issues/36

## long time in munge
# https://github.com/bulik/ldsc/issues/281 add the code --chunksize 500000

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




## munge data
### BCAC_triple
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Trip_ldsc_input_neff.txt \
   --N 107976 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

### BCAC_LuminaB
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_ldsc_input_neff.txt \
   --N 103891 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

### BCAC_LuminaA
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaA_ldsc_input_neff.txt \
   --N 152192 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA8 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

### BCAC_Her2
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_Her2_enrich_ldsc_input_neff.txt \
   --N 97336 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist

### BCAC_LumB_her2
python munge_sumstats_r2_0.8.py \
   --sumstats /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/BCAC_LuminaB_Her2_ldsc_input_neff.txt \
   --N 108420 \
   --chunksize 500000 \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2 \
   --merge-alleles /home/nfs/sunx3/software/ldsc/example/w_hm3.snplist


## ld regression
### LumA and LumB
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB

### LumA and Her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Her2

### LumA and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB-her2

### LumA and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Triple

### LumB and Her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_Her2

### LumB and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_LumB-her2

### LumB and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumB_Triple

### Her2 and LumB_her2
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2_LumB-her2

### Her2 and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2_Triple

### LumB_her2 and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB-her2_Triple

## cross-trait ldsc, convert to liability
### LumA and LumB
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_liability \
  --pop-prev 0.00327,0.0007 \
  --samp-prev 0.2977583,0.0608807


### LumA and LumB_her2
 python ldsc.py \
   --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumB_her2.sumstats.gz \
   --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
   --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_LumB_her2_liability \
   --pop-prev 0.00327,0.0014 \
   --samp-prev 0.2977583,0.0830918

### LumA and her2
 python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Her2.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_her2_liability \
  --pop-prev 0.00327,0.00026 \
  --samp-prev 0.2977583,0.02873369

 
### LumA and Triple
python ldsc.py \
  --rg /home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/LumA.sumstats.gz,/home/nfs/sunx3/project/bra_subtypes_ccgwas/data/ldsc/Trip.sumstats.gz \
  --ref-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --w-ld-chr /home/nfs/sunx3/software/ldsc/example/eur_w_ld_chr/ \
  --out /home/nfs/sunx3/project/bra_subtypes_ccgwas/result/ldsc/LumA_Trip_liability \
  --pop-prev 0.00327,0.00077 \
  --samp-prev 0.2977583,0.08091434
 
