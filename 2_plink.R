# Run plink on plink genotype files
library(data.table)
library(dplyr)
library(tidyr)

# Run plink missingness

system("plink --file data/BGI_out/2.Genotyping/raw --missing --out data/BGI_out/2.Genotyping/raw_miss --allow-extra-chr --debug --no-fid --no-pheno --no-sex --no-parents")

imiss <- read.table("data/BGI_out/2.Genotyping/raw_miss.imiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100)
head(imiss)
hist(imiss$percent_miss)

# read in concentration data

lmiss <- read.table("data/BGI_out/2.Genotyping/raw_miss.lmiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100)

hist(lmiss$percent_miss)

# split ID column

imiss <- imiss %>%
  separate(IID, c("A", "B"), by = "-")

metadata <- read.csv("/Users/emilyhumble/Dropbox/phd/projects/paper9_snp_chip/SAMPLES/data/BGI/Plate1_list_June2018.csv",
                     header = T)

df <- left_join(imiss, metadata, by = c("B" = "Sample_name"))

plot(df$percent_miss~df$GelBand)
plot(df$percent_miss~ df$C_ng_microl)
