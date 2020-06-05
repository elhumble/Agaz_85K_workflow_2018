# Script to calculate error rate from typed SNPs

library(dplyr)

#~~ Extract typed SNPs

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3 ",
              "--extract data/out/agaz/plink/typed_snps.txt ",
              "--out data/out/agaz/plink/agaz_plate1_3_typed ",
              "--make-bed --allow-extra-chr --debug ",
              "--no-fid --no-pheno --no-sex --no-parents"))

# calculate error rate using all typed SNPs

system(paste0("~/software/plink --bfile data/out/agaz/plink/agaz_plate1_3_typed ",
              "--genome --out data/out/agaz/plink/agaz_plate1_3_typed ",
              "--allow-extra-chr --debug --no-fid --no-pheno --no-sex --no-parents"))

error <- fread("data/out/agaz/plink/agaz_plate1_3_typed.genome", header = T) %>%
  filter(Z2 > 0.9) %>%
  select(c(FID1, FID2, Z2))

error

1 - (error %>%
       summarise(mean(Z2)))

# Z2 score: the proportion of SNPs at which two individuals (replicates) share both alleles IBD
