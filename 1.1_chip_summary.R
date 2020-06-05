#~~ SNP chip summary stats for Agaz samples only
#~~ 85,359 SNPs

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(readxl)
source("scripts/theme_emily.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Initial data wrangling     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in all SNP categories

data_path <- "data/out/agaz/"   # path to the data
ps_files <- dir(path = "data/out/agaz/",
                pattern = "\\.ps$")

data <- tibble(filename = ps_files) %>% # create dataframe of file names
  mutate(file_contents = map(filename,      # read in files into new col
                             ~ fread(file.path(data_path, .))))

data <- unnest(data) %>%
  filter(filename != "Recommended.ps") %>%
  filter(filename != "OffTargetVariant.ps") 

total <- nrow(data)

#~~ Get rescued OTVs

data_path <- "data/out/agaz/OTV"   # path to the data
ps_files <- dir(path = "data/out/agaz/OTV/",
                pattern = "\\.ps$")
ps_files <- ps_files[!grepl("OTV", ps_files)]

otv <- tibble(filename = ps_files) %>% # create dataframe of file names
  mutate(file_contents = map(filename,      # read in files into new col
                             ~ fread(file.path(data_path, .))))
otv <- unnest(otv) %>%
  filter(filename != "Recommended.ps")


#~~ Join

data <- rbind(data, otv) 

#~~ Summaries

data %>%
  group_by(filename) %>%
  summarise(number = n(), percentage = (number / total)*100)


#~~ Add custom IDs using annotation file

annot <- read.csv("data/library_files/Axiom_Agaz_90K_Annotation.r1.csv", skip = 18) %>%
  select(Probe.Set.ID, cust_id)

data <- left_join(data, annot,
                  by = c("probeset_id" = "Probe.Set.ID")) %>%
  mutate(type = case_when(grepl("Contig", cust_id) ~ "RAD",
                          grepl("^Ag", cust_id) | grepl("^4", cust_id) ~ "Transcriptome",
                          grepl("JIH", cust_id) ~ "JIH",
                          TRUE ~ "Canine"))

#~~ Write list of otvs only for downstream comparison

otv %>%
  left_join(annot,                  
            by = c("probeset_id" = "Probe.Set.ID")) %>%
  mutate(type = case_when(grepl("Contig", cust_id) ~ "RAD",
                          grepl("^Ag", cust_id) | grepl("^4", cust_id) ~ "Transcriptome",
                          grepl("JIH", cust_id) ~ "JIH",
                          TRUE ~ "Canine")) %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/agaz/plink/otv_snps.txt", quote = F,
              row.names = F, col.names = F)

# Add affy info and priority categories from design phase

pri <- fread("data/processed/axiomdesign_Seq_FurSeal_Prioritised_Consideration_01Feb2018") %>%
  select(Snpid, SNP_PRIORITY, forwardPconvert, reversePconvert)

# Summary of tiled SNPs by RAD / Trans / pre-val

pri %>%
  left_join(data, by = c("Snpid" = "cust_id")) %>%
  group_by(type) %>%
  summarise(n = n())
  
# NA = RAD
# Transcriptome includes 6 MHC

# Add priority categories and convert scores to main dataframe

data <- data %>%
  left_join(pri, by = c("cust_id" = "Snpid")) %>%
  mutate(forwardPconvert = as.numeric(forwardPconvert),
         reversePconvert = as.numeric(reversePconvert))

# Split into polymorphic and all typed SNPs

poly <- data %>%
  filter(filename == "PolyHighResolution.ps" | filename == "NoMinorHom.ps")

typed <- data %>%
  filter(filename == "PolyHighResolution.ps" | filename == "NoMinorHom.ps" | filename == "MonoHighResolution.ps")


#~~~~~~~~~~~~~~~~~#
#     Numbers     #
#~~~~~~~~~~~~~~~~~#

# Proportion polymorphic over total tiled
nrow(poly) / 85359 * 100

# Proportion succesfully typed over total tiled
nrow(typed) / 85359 * 100

# Proportion polymorphic over total typed
nrow(poly) / nrow(typed) * 100


nrow(filter(poly, type == "RAD")) # n poly RAD SNPs
nrow(filter(poly, type == "Transcriptome")) # n poly Transcriptome SNPs
nrow(filter(poly, type == "Canine" | type == "JIH")) # n poly pre-validated SNPs
nrow(filter(poly, type == "Canine")) # n of poly Canine BeadChip SNPs
nrow(filter(poly, type == "JIH")) # n of poly pre-val trans


# Write lists of polymorphic and typed SNPs for downstream analysis:

poly %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/agaz/plink/poly_snps.txt", quote = F,
              row.names = F, col.names = F)

typed %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/agaz/plink/typed_snps.txt", quote = F,
              row.names = F, col.names = F)

# poly hi res

poly %>%
  filter(filename == "PolyHighResolution.ps") %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/agaz/plink/poly_hi_res_snps.txt", quote = F,
              row.names = F, col.names = F)

# Average design score for printed SNPs?

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Exploring array performance by SNP priority category  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Break down of all tiled SNPs

data %>%
  group_by(SNP_PRIORITY) %>%
  summarise(n())

a <- data %>%
  group_by(SNP_PRIORITY) %>%
  summarise(total = n()) %>%
  select(total)

# Break down of all typed SNPs

typed %>%
  group_by(SNP_PRIORITY) %>%
  summarise(total = n())

b <- typed %>%
  group_by(SNP_PRIORITY) %>%
  summarise(total = n()) %>%
  select(total)

# Break down of all polymorphic SNPs

poly %>%
  group_by(SNP_PRIORITY) %>%
  summarise(total = n())


#~~ Proportion typed / tiled across priorities 1:7 (conversion rate)

(b$total / a$total) * 100
sum(b$total) / sum(a$total) # Overall


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Exploring array performance by SNP type         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Conversion rate = proportion of SNPs yielding high quality genotypes

# Break down of all tiled SNPs

data %>%
  group_by(type) %>%
  summarise(n())

a <- data %>%
  group_by(type) %>%
  summarise(total = n()) %>%
  select(total)

# Break down of all typed SNPs

typed %>%
  group_by(type) %>%
  summarise(n())

b <- typed %>%
  group_by(type) %>%
  summarise(total = n()) %>%
  select(total)

# Break down of all polymorphic SNPs

poly %>%
  group_by(type) %>%
  summarise(n())

#~~ Proportion typed / tiled across SNP categories

b$total / a$total


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Conversion rate of pre-validated success rate     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Conversion rate = proportion of SNPs yielding high quality genotypes

# Read in list of 40 validated RAD SNP IDs from Humble et al 2016

val_rad <- read_excel("data/processed/validated_RAD_SNPs.xlsx") %>%
  filter(SequencingOutcome == 1) %>%
  select(Contig, Location) %>%
  unite("cust_id", c("Contig", "Location"))

# Conversion rate of all pre-validated SNPs

nrow(typed %>%
       filter(type == "Canine" | type == "JIH" | cust_id %in% val_rad$cust_id)) /
  nrow(data %>%
         filter(type == "Canine" | type == "JIH" | cust_id %in% val_rad$cust_id)) * 100

# Conversion rate of JIH and Canine

nrow(typed %>%
       filter(type == "Canine" | type == "JIH")) /
  nrow(data %>%
         filter(type == "Canine" | type == "JIH")) * 100

# Conversion rate of Canine

nrow(typed %>%
       filter(type == "Canine" )) /
  nrow(data %>%
         filter(type == "Canine" )) * 100 # 58 / 173

# Conversion rate of JIH

nrow(typed %>%
       filter(type == "JIH")) /
  nrow(data %>%
         filter(type == "JIH")) * 100

# Conversion rate of JIH and 40 RAD pre-val

nrow(typed %>%
       filter(type == "JIH" | cust_id %in% val_rad$cust_id)) /
  nrow(data %>%
         filter(type == "JIH" | cust_id %in% val_rad$cust_id)) * 100

# Conversion rate of 40 pre-val RAD SNPs

length(val_rad$cust_id[val_rad$cust_id %in% typed$cust_id]) / 
  length(val_rad$cust_id[val_rad$cust_id %in% data$cust_id]) * 100


#~~ Conversion rate of unvalidated SNPs

# All minus validated SNPs including 40 pre-val RAD

nrow(typed %>%
       filter(!cust_id %in% val_rad$cust_id) %>%
       filter(type == "RAD" | type == "Transcriptome")) /
  nrow(data %>%
         filter(!cust_id %in% val_rad$cust_id) %>%
         filter(type == "RAD" | type == "Transcriptome")) * 100


#~~ Canine SNPs less likely to convert than unvalidated SNPs?

# Test that the probabilities for the two groups are equal

canine_unval <- matrix(c(nrow(typed %>% filter(type == "Canine")), # typed canine
                         
                         nrow(data %>%  filter(type == "Canine" )) - 
                           nrow(typed %>% filter(type == "Canine")), # untyped canine
                         
                         nrow(typed %>%
                                filter(!cust_id %in% val_rad$cust_id) %>%
                                filter(type == "RAD" | type == "Transcriptome")), # typed unval
                         nrow(data %>%
                                filter(!cust_id %in% val_rad$cust_id) %>%
                                filter(type == "RAD" | type == "Transcriptome")) -
                           nrow(typed %>%
                                  filter(!cust_id %in% val_rad$cust_id) %>%
                                  filter(type == "RAD" | type == "Transcriptome"))), # untyped unval
                       ncol = 2, byrow = TRUE)

colnames(canine_unval) <- c("Typed", "Untyped")
rownames(canine_unval) <- c("Canine", "Unvalidated") 
canine_unval

# 2-sample test for equality of proportions without continuity correction

prop.test(x = canine_unval, correct = FALSE)

# Chi-squared test (same compution, different function)

chisq.test(canine_unval, correct = F)
chi <- chisq.test(canine_unval, correct = F)
chi$expected # (>5?)


#~~ Remaining validated SNPs still less likely to convert than unvalidated SNPs?

val_unval <- matrix(c(nrow(typed %>%
                             filter(type == "JIH" | cust_id %in% val_rad$cust_id)), # typed val - canine
                      
                      nrow(data %>%
                             filter(type == "JIH" | cust_id %in% val_rad$cust_id)) -
                        nrow(typed %>%
                               filter(type == "JIH" | cust_id %in% val_rad$cust_id)), # untyped val - untyped canine
                      
                      nrow(typed %>%
                             filter(!cust_id %in% val_rad$cust_id) %>%
                             filter(type == "RAD" | type == "Transcriptome")), # typed unval
                      
                      nrow(data %>%
                             filter(!cust_id %in% val_rad$cust_id) %>%
                             filter(type == "RAD" | type == "Transcriptome")) -
                        nrow(typed %>%
                               filter(!cust_id %in% val_rad$cust_id) %>%
                               filter(type == "RAD" | type == "Transcriptome"))), # untyped unval
                    ncol = 2, byrow = TRUE)

colnames(val_unval) <- c("Typed", "Untyped")
rownames(val_unval) <- c("Remaining validated", "Unvalidated") 
val_unval


# 2-sample test for equality of proportions without continuity correction

prop.test(x = val_unval, correct = FALSE)

# Chi-squared test (same compution, different function)

chisq.test(val_unval, correct = F)
chi <- chisq.test(val_unval, correct = F)
chi$expected # (>5?)



#~~ Classification of pre-validated SNPs that did not convert

val <- data %>%
  filter(type == "Canine" | type == "JIH" | cust_id %in% val_rad$cust_id) %>%
  filter(filename == "Other.ps" | filename == "CallRateBelowThreshold.ps" | filename == "OffTargetVariant.ps")

ggplot(val, aes(filename, fill = type)) +
  geom_bar()


#~~ Number of annotated SNPs

transcriptome <- filter(typed, type == "Transcriptome") %>%
  mutate(cust_id = gsub("_v1", ".v1", cust_id)) %>%
  separate(cust_id, c("Contig", "Position"), sep = "_") %>%
  mutate(Contig = gsub(".v1", "_v1", Contig))

trans_annot <- fread("data/processed/Transcript_annotations.txt") %>% # This file is saved on HD: server data
  gather(category, annotation, -Contig_Name) %>%
  filter(!is.na(annotation)) %>%
  left_join(transcriptome, by = c("Contig_Name" = "Contig")) %>%
  filter(!is.na(filename))

trans_annot %>%
  group_by(category) %>%
  summarise(n())

# Meinolf's 6 transcripts are annotated as both immune and MHC.
# Most MHC SNPs are also immune. But not all.
# Not all immune SNPs are MHC

# 5/6 MHC SNPs typed. 3 polymorphic

# Thank you to Meinolf Ottensmann for designing the six MHC SNPs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          RAD SNP stats           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Scaffold lengths

Agaz_v1.2_lengths <- read.table("data/PBJelly_lengths.txt")

# Combine with typed RAD SNP information

RAD_SNPs <- filter(typed, grepl("Con", cust_id)) %>%
  separate(cust_id, c("Contig", "Location")) %>%
  left_join(Agaz_v1.2_lengths, by = c("Contig" = "V1"))

# no of Scaffolds with SNPs
length(unique(RAD_SNPs$Contig))

# no of SNPs / Scaffold length

snps_per_contig <- RAD_SNPs %>%
  dplyr::group_by(V2) %>%
  dplyr::summarise(n = n()) %>%
  mutate(LengthKb = V2 / 1000)

cor(snps_per_contig$LengthKb, snps_per_contig$n) # pearson

png(file="figs/RAD_SNPs_per_Scaffold.png", units = "in", res = 300, height=7, width=8)

ggplot(snps_per_contig, aes(x=LengthKb, y = n)) +
  geom_point(alpha=0.9, col = "grey30") +
  labs(y = "Number of SNPs", x = "Scaffold Length (kb)") +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme_emily()

dev.off()


#~~ SNP spacing

RAD_SNPs <- RAD_SNPs %>%
  mutate(Location = as.numeric(Location))

# split df into list of Contigs
dist_snp <- split(RAD_SNPs , f = RAD_SNPs$Contig)

# arrange by SNP order within each contig
dist_snp <- lapply(dist_snp, function(x) 
  arrange(x, Location)) 

dist_snp <- lapply(dist_snp, function(x) 
  mutate(x, distleft = (x$Location - lag(x$Location))))

# Get distance of first SNP on contig

dist_snp <- lapply(dist_snp, function(x)
  mutate(x, distleft = ifelse(is.na(distleft),Location,distleft)))

head(dist_snp[[1]])
library(plyr)
dist_snp <- ldply(dist_snp)


nrow(filter(dist_snp, distleft < 10000))
nrow(filter(dist_snp, distleft < 100))
nrow(filter(dist_snp, distleft > 10000)) # 10 kb
nrow(filter(dist_snp, distleft > 100000)) # 100 kb

summary(dist_snp$distleft/1000)

# Plot distances

dist_snp_plot <- dist_snp %>%
  mutate(distleft = ifelse(distleft > 300000, 300000, distleft)) %>%
  mutate(distleftKb = distleft / 1000)

png(file="figs/typed_RAD_SNP_distances.png", units = "in", res = 300, height=5, width=6)

ggplot(dist_snp_plot, aes(x=distleftKb)) +
  geom_histogram(alpha=0.8, fill = "grey30") +
  labs(y = "Frequency", x = "SNP distance (kb)") +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels=c(seq(0,200, by=100), "300+")) +
  theme_emily()

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          RAD SNP POLY stats           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Combine with polu RAD SNP information

RAD_SNPs <- filter(poly, grepl("Con", cust_id)) %>%
  separate(cust_id, c("Contig", "Location")) %>%
  left_join(Agaz_v1.2_lengths, by = c("Contig" = "V1"))

# no of Scaffolds with SNPs
length(unique(RAD_SNPs$Contig))

# no of SNPs / Scaffold length

snps_per_contig <- RAD_SNPs %>%
  dplyr::group_by(V2) %>%
  dplyr::summarise(n = n()) %>%
  mutate(LengthKb = V2 / 1000)

cor(snps_per_contig$LengthKb, snps_per_contig$n) # pearson

png(file="figs/RAD_SNPs_per_Scaffold.png", units = "in", res = 300, height=7, width=8)

ggplot(snps_per_contig, aes(x=LengthKb, y = n)) +
  geom_point(alpha=0.9, col = "grey30") +
  labs(y = "Number of SNPs", x = "Scaffold Length (kb)") +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme_emily()

dev.off()

length(unique(RAD_SNPs$Contig))


#~~ SNP spacing

RAD_SNPs <- RAD_SNPs %>%
  mutate(Location = as.numeric(Location))

# split df into list of Contigs
dist_snp <- split(RAD_SNPs , f = RAD_SNPs$Contig)

# arrange by SNP order within each contig
dist_snp <- lapply(dist_snp, function(x) 
  arrange(x, Location)) 

dist_snp <- lapply(dist_snp, function(x) 
  mutate(x, distleft = (x$Location - lag(x$Location))))

# Get distance of first SNP on contig

dist_snp <- lapply(dist_snp, function(x)
  mutate(x, distleft = ifelse(is.na(distleft),Location,distleft)))

head(dist_snp[[1]])
library(plyr)
dist_snp <- ldply(dist_snp)


nrow(filter(dist_snp, distleft < 10000))
nrow(filter(dist_snp, distleft < 100))
nrow(filter(dist_snp, distleft > 10000)) # 10 kb
nrow(filter(dist_snp, distleft > 100000)) # 100 kb

summary(dist_snp$distleft/1000)

# Plot distances

dist_snp_plot <- dist_snp %>%
  mutate(distleft = ifelse(distleft > 300000, 300000, distleft)) %>%
  mutate(distleftKb = distleft / 1000)

png(file="figs/poly_RAD_SNP_distances.png", units = "in", res = 300, height=7, width=8)

ggplot(dist_snp_plot, aes(x=distleftKb)) +
  geom_histogram(alpha=0.8, fill = "grey30") +
  labs(y = "Frequency", x = "SNP distance (kb)") +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(labels=c(seq(0,200, by=100), "300+")) +
  theme_emily()

dev.off()

