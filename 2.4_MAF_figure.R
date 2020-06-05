library(data.table)
library(dplyr)
source("scripts/theme_emily.R")
library(gridExtra)
library(cowplot)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(patchwork)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Estimated MAF            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Prepare dataframe of estimated MAFs
# This is using code from chip proposal workflow

annotations <- fread("data/processed/SNP_info_df_27_01_2018.txt") %>%
  unite("cust_id", Contig, Position)

proposal <- fread("data/processed/proposal_2018-02-01.txt" ,header = T) %>%
  filter(Tile == 1)

est_maf <- left_join(proposal, annotations, by = "cust_id")

# Read in MAFs for transcriptome data

bowtie <- read.csv("data/processed/bowtie_MAFs.csv")
bwa <- read.csv("data/processed/bwa_MAFs.csv")
swap <- read.csv("data/processed/swap454SNPs.csv")
newbler <- read.csv("data/processed/NewblerSNPs.csv")

transcriptome_mafs <- bowtie %>%
  full_join(bwa, by = c("Contig_Name", "SNP_Position")) %>%
  full_join(swap, by = c("Contig_Name", "SNP_Position")) %>%
  full_join(newbler, by = c("Contig_Name", "SNP_Position")) %>%
  unite(cust_id, c("Contig_Name", "SNP_Position"), sep = "_")

transcriptome_mafs <- est_maf %>%
  left_join(transcriptome_mafs, by = "cust_id") %>%
  filter(!grepl("Con", cust_id)) %>%
  filter(!grepl("^JIH", cust_id)) %>%
  filter(!grepl("^Canis", organism)) %>%
  select(cust_id, MAF.x, MAF.y, MAF.x.x, MAF.y.y) %>%
  dplyr::mutate(MAF = MAF.x) %>%
  dplyr::mutate(MAF = ifelse(is.na(MAF), MAF.y,
                             MAF)) %>%
  dplyr::mutate(MAF = ifelse(is.na(MAF), MAF.y.y,
                             MAF)) %>%
  dplyr::mutate(MAF = ifelse(is.na(MAF), MAF.x.x,
                             MAF)) %>%
  filter(!is.na(MAF)) %>%
  select(cust_id, MAF)


#~~ Combine transcriptome MAFs

est_maf <- left_join(est_maf, transcriptome_mafs, by = "cust_id") %>%
  mutate(MAF = MAF.x) %>%
  mutate(MAF = ifelse(is.na(MAF), MAF.y,
                      MAF)) %>%
  select(-c(MAF.x, MAF.y)) %>%
  mutate(type = ifelse(grepl("^Con", cust_id), "RAD", "trans")) %>%
  select(cust_id, MAF)

colnames(est_maf) <- c("cust_id", "est_MAF")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Observed MAF             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

obs_maf <- read.table("data/out/agaz/plink/agaz_plate1_3_poly_parents.frq", header = T)

hist(obs_maf$MAF)

obs_maf <- obs_maf %>%
  mutate(type = case_when(grepl("Contig", CHR) ~ "RAD",
                          TRUE ~ "Transcriptome"))
mean(obs_maf$MAF)

obs_maf %>%
  filter(type == "Transcriptome") %>%
  summarise(mean = mean(MAF),
            sd = sd(MAF))

obs_maf %>%
  filter(type == "RAD") %>%
  summarise(mean = mean(MAF),
            sd = sd(MAF))

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Plot           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

mafs <- obs_maf %>%
  left_join(est_maf, by = c("SNP" = "cust_id"))

#~~ Correlations

pal <- wes_palette("Darjeeling2")

mafs <- mafs %>%
  mutate(type = case_when(grepl("Con", SNP) ~ "RAD",
                          TRUE ~ "Transcriptome"),
         plot = case_when(grepl("Transcriptome", type) ~ "D",
                          TRUE ~ "C"))

corr <- ggplot(mafs, aes(est_MAF, MAF, col = plot)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method=lm, se=T, col = "grey20") +
  theme_emily() +
  facet_wrap(~plot, nrow = 2) +
  ylab("Empirical MAF") + xlab(expression(italic("in silico")~textstyle(MAF))) +
  theme(legend.position = "none",
        strip.text.y = element_blank(),
        strip.text = element_text(hjust = 0,
                                  face = "plain", 
                                  size = "16")) +
  scale_color_manual(values = pal[c(2,3)])

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}

corr_eqn(filter(mafs, type == "RAD")$est_MAF,
         filter(mafs, type == "RAD")$MAF)

corr_eqn(filter(mafs, type == "Transcriptome" & !is.na(est_MAF))$est_MAF,
         filter(mafs, type == "Transcriptome" & !is.na(est_MAF))$MAF)

#~~ Distributions

pal <- c("#D69C4E","#046C9A","#eed7b8","#cce1ea")

distB <- mafs %>%
  select(SNP, MAF, type, est_MAF) %>%
  gather("obs_est", "MAF", -c(SNP, type)) %>%
  unite("plot", c("type", "obs_est"), remove = F) %>%
  mutate(title = case_when(type == "RAD" ~ "A",
                           TRUE ~ "B")) %>%
  mutate(alpha = case_when(plot == "RAD_MAF" | plot == "Transcriptome_MAF" ~ 0.1,
                           TRUE ~ 0.7)) %>%
  ggplot(aes(x=MAF, fill = plot)) +
  #geom_histogram(binwidth = 0.01) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.01) +
  labs(y = "Density", x = "Minor Allele Frequency") +
  theme_emily() +
  facet_wrap(~title, nrow = 2) +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(hjust = 0,
                                  face = "plain", 
                                  size = "16"),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.9)) +
  scale_fill_manual(breaks=c("RAD_MAF", "Transcriptome_MAF", "RAD_est_MAF", "Transcriptome_est_MAF"),
                    labels=c("Empirical MAF", "Empirical MAF", expression(italic("in silico")~textstyle(MAF)), expression(italic("in silico")~textstyle(MAF))),
                    values = pal[c(2,1,4,3)])
distB



fig_2 <- distB + corr
ggsave("figs/Figure_2.tiff", fig_2, height = 15, width = 18, units = "cm")
       


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Using all original MAF info            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


orig_trans_mafs <- bowtie %>%
  full_join(bwa, by = c("Contig_Name", "SNP_Position")) %>%
  full_join(swap, by = c("Contig_Name", "SNP_Position")) %>%
  full_join(newbler, by = c("Contig_Name", "SNP_Position")) %>%
  unite(cust_id, c("Contig_Name", "SNP_Position"), sep = "_") %>%
  select(cust_id, MAF.x) %>%
  mutate(type = "Transcriptome_discovery")

colnames(orig_trans_mafs) <- c("SNP", "MAF", "type")

# this file comes from ArcGaz_biallelic_sub.recode.vcf filtered only for depth and not maf (server pipeline hard drive)
# no other filters were applied
# plink frq applied
# all on eddie
# module load igmm/apps/vcftools/0.1.13
# vcftools --vcf ArcGaz_biallelic_sub.recode.vcf --out ArcGaz_biallelic_sub_qual_miss --min-meanDP 5 --max-meanDP 18 --max-missing 0.6 --recode --recode-INFO-all
# grep -v '#' ArcGaz_biallelic_sub_qual_miss.recode.vcf | wc -l
# vcftools --vcf data/ArcGaz_biallelic_sub_qual_miss.recode.vcf --plink --out data/ArcGaz_biallelic_sub_qual_miss
# module load roslin/plink/1.90p
# plink --file ArcGaz_biallelic_sub_qual_miss --freq --out ArcGaz_biallelic_sub_qual_miss --allow-extra-chr --debug --nonfounders

# in processed:
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_nov_19/ArcGaz_biallelic_sub_qual_miss.frq ./

orig_rad_mafs <- fread("data/processed/ArcGaz_biallelic_sub_qual_miss.frq") %>%
  select(SNP, MAF) %>%
  mutate(type = "RAD_discovery") %>%
  filter(MAF > 0)
 
colnames(orig_rad_mafs) <- c("SNP", "MAF", "type")

typed_mafs <- mafs %>%
  select(SNP, MAF) %>%
  mutate(type = case_when(grepl("Contig", SNP) ~ "typed_RAD",
                          TRUE ~ "typed_Transcriptome"))

df <- rbind(orig_trans_mafs, orig_rad_mafs, typed_mafs) %>%
  mutate(facet = case_when(grepl("Transcriptome", type) ~ "transcriptome",
                   TRUE ~ "RAD"))

png(file="figs/MAF_distribution.png", units = "in", res = 300, height=5, width=9)
ggplot(df, aes(x=MAF, fill = type)) +
  geom_histogram(binwidth = 0.01) +
  labs(y = "Number of SNPs", x = "Minor Allele Frequency") +
  theme_emily() +
  #scale_y_log10() +
  facet_wrap(~facet, scales = "free")
dev.off()
#legend.position = c(0.8, 0.8)

