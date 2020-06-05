#~~ Summary statistics for polymoprhic Agaz SNPs

library(data.table)
library(dplyr)
library(ggplot2)
source("scripts/theme_emily.R")
library(wesanderson)

#~~ Get missingness and MAF

# relatives (offspring) in calculation and remove dups
system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly ", 
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--missing --freq ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly ",
              "--allow-extra-chr --debug"))

# no relatives (offspring removed) in calculation and remove dups
system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ", 
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--missing --freq ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--allow-extra-chr --debug"))

#~~ Compare MAFs

relatives_frq <- fread("data/out/agaz/plink/agaz_plate1_3_poly.frq") %>%
  select(SNP, MAF)
colnames(relatives_frq) <- c("SNP", "MAF_relatives")

no_relatives_frq <- fread("data/out/agaz/plink/agaz_plate1_3_poly_parents.frq") %>%
  select(SNP, MAF)
colnames(no_relatives_frq) <- c("SNP", "MAF_no_relatives")


all_frq <- left_join(relatives_frq, no_relatives_frq, by = "SNP")

pal <- wes_palette("Darjeeling2")

maf_compare <- ggplot(all_frq) +
  geom_point(aes(MAF_relatives, MAF_no_relatives), 
             size = 1, alpha = 1/10) +
  theme_emily() +
  xlab("MAF relatives") + ylab("MAF no relatives")
maf_compare

ggsave("figs/MAF_comparison.png", maf_compare, width = 10, height = 10, units = "cm")


#~~ Individual missingness

imiss <- read.table("data/out/agaz/plink/agaz_plate1_3_poly_parents.imiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100) %>%
  separate(IID, c("well", "Sample_name", "array", "plate"))

head(imiss)
hist(imiss$percent_miss)

#~~ SNP call rates

# All polymorphic SNPs

lmiss <- read.table("data/out/agaz/plink/agaz_plate1_3_poly_parents.lmiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100)

head(lmiss)
hist(lmiss$percent_miss)
summary(100 - lmiss$percent_miss)


# Call rates for OTVs

otv <- fread("data/out/agaz/plink/otv_snps.txt", header = F)

otv <- inner_join(otv, lmiss, by = c("V1" = "SNP"))
summary(100 - otv$percent_miss) # call rate for rescued OTV SNPs

# Call rates without OTVs

no_otv <- anti_join(lmiss, otv, by = c("SNP" = "V1"))
summary(100 - no_otv$percent_miss) # call rate for non-OTV SNPs

sum(otv$V1 %in% lmiss$SNP)


