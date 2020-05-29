#~~ Summary statistics for converted Agaz SNPs
#~~ Call rate, error rate

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(ggthemr)
library(gridExtra)
source("theme_emily.R")
library(wesanderson)


#ggthemr(palette = "pale", layout = "clean", 
#        line_weight = 0.7, text_size = 20, type = "outer")

#~~ Filter plink file for polymorphic and typed SNPs

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3 ", 
              "--extract data/out/agaz/plink/poly_snps.txt ",
              "--missing --freq ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly ",
              "-recode --allow-extra-chr --debug --nonfounders ",
              "--no-fid --no-pheno --no-sex --no-parents"))


system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3 ",
              "--extract data/out/agaz/plink/typed_snps.txt ",
              "--out data/out/agaz/plink/agaz_plate1_3_typed ",
              "--make-bed --allow-extra-chr --debug --nonfounders ",
              "--no-fid --no-pheno --no-sex --no-parents"))


#~~ Individual missingness

imiss <- read.table("data/out/agaz/plink/agaz_plate1_3_poly.imiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100) %>%
  separate(IID, c("well", "Sample_name", "array", "plate"))

head(imiss)
hist(imiss$percent_miss)

#~~ SNP missingness

lmiss <- read.table("data/out/agaz/plink/agaz_plate1_3_poly.lmiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100)

head(lmiss)
hist(lmiss$percent_miss)


# Error rate using all typed SNPs

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



