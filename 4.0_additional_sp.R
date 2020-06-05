#~~ Additional species analysis

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(adegenet)
library(ggthemr)
library(wesanderson)

#~~ Summary statistics

data_path <- "data/out/all/"   # path to the data
ps_files <- dir(path = "data/out/all/",
                pattern = "\\.ps$")

data <- tibble(filename = ps_files) %>% # create dataframe of file names
  mutate(file_contents = map(filename,      # read in files into new col
                             ~ fread(file.path(data_path, .))))

data <- unnest(data) %>%
  filter(filename != "Recommended.ps") %>%
  filter(filename != "OffTargetVariant.ps") 

total <- nrow(data)

#~~ Get rescued OTVs

data_path <- "data/out/all/OTV"   # path to the data
ps_files <- dir(path = "data/out/all/OTV/",
                pattern = "\\.ps$")
ps_files <- ps_files[!grepl("OTV", ps_files)]

otv <- tibble(filename = ps_files) %>% # create dataframe of file names
  mutate(file_contents = map(filename,      # read in files into new col
                             ~ fread(file.path(data_path, .))))
otv <- unnest(otv) %>%
  filter(filename != "Recommended.ps")


#~~ Join

data <- rbind(data, otv) 

data %>%
  group_by(filename) %>%
  summarise(number = n(), percentage = (number / total)*100)


#~~ Add extra info from annotation file

annot <- read.csv("data/library_files/Axiom_Agaz_90K_Annotation.r1.csv", skip = 18) %>%
  select(Probe.Set.ID, cust_id)

data <- left_join(data, annot,
                  by = c("probeset_id" = "Probe.Set.ID")) %>%
  mutate(type = case_when(grepl("Contig", cust_id) ~ "RAD",
                          grepl("^Ag", cust_id) | grepl("^4", cust_id) ~ "Transcriptome",
                          grepl("JIH", cust_id) ~ "JIH",
                          TRUE ~ "Canine"))

#~~ Add affy info and priority categories

pri <- fread("data/processed/axiomdesign_Seq_FurSeal_Prioritised_Consideration_01Feb2018") %>%
  select(Snpid, SNP_PRIORITY, forwardPconvert, reversePconvert)

data <- data %>%
  left_join(pri, by = c("cust_id" = "Snpid")) %>%
  mutate(forwardPconvert = as.numeric(forwardPconvert),
         reversePconvert = as.numeric(reversePconvert))

#~~ Filter for good and polymorphic SNPs

poly <- data %>%
  filter(filename == "PolyHighResolution.ps" | filename == "NoMinorHom.ps")

typed <- data %>%
  filter(filename == "PolyHighResolution.ps" | filename == "NoMinorHom.ps" | filename == "MonoHighResolution.ps")

nrow(poly) / 85359 * 100
nrow(typed) / 85359 * 100


#~~ Filter plink file for polymorphic SNPs

poly %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/all/plink/poly_snps.txt", quote = F,
              row.names = F, col.names = F)

system("~/software/plink --file data/out/all/plink/agaz_plate1_3 --extract data/out/all/plink/poly_snps.txt --freq --make-bed --recodeA --out data/out/all/plink/agaz_plate1_3_poly --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")

typed %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>% # edit JIH SNP names
  select(cust_id) %>%
  write.table("data/out/all/plink/typed_snps.txt", quote = F,
              row.names = F, col.names = F)

system("~/software/plink --file data/out/all/plink/agaz_plate1_3 --extract data/out/all/plink/typed_snps.txt --freq --make-bed --out data/out/all/plink/agaz_plate1_3_typed --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")


#~~ Individual call rate from full PLINK file

system("~/software/plink --bfile data/out/all/plink/agaz_plate1_3_typed --missing --out data/out/all/plink/agaz_plate1_3_typed --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")

imiss <- read.table("data/out/all/plink/agaz_plate1_3_typed.imiss", header = T) %>%
  mutate(percent_miss = F_MISS * 100,
         N_CALL = N_GENO - N_MISS)

imiss %>%
  filter(grepl("RR", IID)) %>%
  summarise(mean=mean(100 - percent_miss),
            meanN = mean(N_CALL),
            minN = min(N_CALL),
            maxN = max(N_CALL))

imiss %>%
  filter(grepl("CA", IID)) %>%
  summarise(mean=mean(100 - percent_miss),
            meanN = mean(N_CALL),
            minN = min(N_CALL),
            maxN = max(N_CALL))

imiss %>%
  filter(!grepl("RR|CA", IID)) %>%
  summarise(mean=mean(100 - percent_miss),
            meanN = mean(N_CALL),
            minN = min(N_CALL),
            maxN = max(N_CALL))

#~~ Alternative call rate:

# Call rate

call_rate <- fread("data/out/all/AxiomGT1.report.txt", header = T, skip = 364)

call_rate <- call_rate %>%
  mutate(sp = case_when(grepl("CA", cel_files) ~ "CA",
                        grepl("RR", cel_files) ~ "RR",
                        TRUE ~ "AG"))

ggplot(call_rate, aes(sp, call_rate)) +
  geom_boxplot()

call_rate %>%
  filter(grepl("RR", sp)) %>%
  summarise(mean=mean(call_rate))

call_rate %>%
  filter(grepl("CA", sp)) %>%
  summarise(mean=mean(call_rate))

# Basically the same

# Time to most recent common ancestor

# CA: GSL vs AG
# RR: Stellar vs AG 5.2 my
# GS: Grey vs AG 23 my 5.2 my

#~~ Prepare ID files for structure

imiss %>%
  filter(grepl("RR", FID)) %>%
  select(FID) %>%
  mutate(IID = FID) %>%
  write.table("data/out/all/plink/SSL.txt", quote = F, row.names = F)

imiss %>%
  filter(grepl("CA", FID)) %>%
  select(FID) %>%
  mutate(IID = FID) %>%
  write.table("data/out/all/plink/GSL.txt", quote = F, row.names = F)


imiss %>%
  filter(!grepl("GS|RR|CA", FID)) %>%
  select(FID) %>%
  mutate(IID = FID) %>%
  write.table("data/out/all/plink/AFS.txt", quote = F, row.names = F)

#~~ Get SNP freqs for each species group

system("~/software/plink --bfile data/out/all/plink/agaz_plate1_3_typed --keep data/out/all/plink/SSL.txt --freq --out data/out/all/plink/SSL_typed --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")
system("~/software/plink --bfile data/out/all/plink/agaz_plate1_3_typed --keep data/out/all/plink/GSL.txt --freq --out data/out/all/plink/GSL_typed --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")
system("~/software/plink --bfile data/out/all/plink/agaz_plate1_3_typed --keep data/out/all/plink/AFS.txt --freq --out data/out/all/plink/AFS_typed --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")


#~~ Count number of polymorphic SNPs

ssl_freq <- fread("data/out/all/plink/SSL_poly.frq")
gsl_freq <- fread("data/out/all/plink/GSL_poly.frq")
afs_freq <- fread("data/out/all/plink/AFS_poly.frq")

(nrow(ssl_freq %>%
        filter(MAF > 0.1)) / nrow(ssl_freq)) * 100

(nrow(ssl_freq %>%
        filter(MAF > 0.1)) / 74130 ) * 100

(nrow(gsl_freq %>%
        filter(MAF > 0.1)) / nrow(gsl_freq)) * 100

(nrow(gsl_freq %>%
        filter(MAF > 0.1)) / 73922) * 100

#(nrow(gs_freq %>%
#       filter(MAF > 0.1)) / nrow(gs_freq)) * 100

(nrow(afs_freq %>%
        filter(MAF > 0.1)) / nrow(afs_freq)) * 100


#~~ Filter individuals from plink file

set.seed(5)
agaz_struc <- imiss %>% # Randomly select 4 adult AGF
  filter(!grepl("RR|GS|CA|AGP", FID)) %>%
  select(FID) %>%
  mutate(IID = FID) %>%
  sample_n(4)

imiss %>%
  filter(grepl("RR|CA", FID)) %>% # GS
  select(FID) %>%
  mutate(IID = FID) %>%
  rbind(agaz_struc) %>%
  write.table("data/out/all/struc_inds.txt", quote = F, row.names = F)

#~~ Get intersect of polymorphic SNPs

gsl_poly <- gsl_freq %>%
  filter(MAF != 0) %>%
  select(SNP)

ssl_poly <- ssl_freq %>%
  filter(MAF != 0) %>%
  select(SNP)

#gs_poly <- gs_freq %>%
#  filter(MAF != 0) %>%
#  select(SNP)

afs_poly <- afs_freq %>%
  filter(MAF != 0) %>%
  select(SNP)

intersect(gsl_poly, ssl_poly, afs_poly) %>%
  write.table("data/out/all/intersect_snps.txt", quote = F,
              row.names = F, col.names = F)

#~~ Filter structure file for inds and intersect 

system("~/software/plink --file data/out/all/agaz_plate1-3 --extract data/out/all/intersect_snps.txt --keep data/out/all/struc_inds.txt --recodeA --out data/out/all/struc_inds --geno 0.1 --maf 0.01 --hwe 0.001 --allow-extra-chr --make-bed --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")
system("~/software/plink --bfile data/out/all/struc_inds --freq --out data/out/all/struc_inds --allow-extra-chr --debug --nonfounders --no-fid --no-pheno --no-sex --no-parents")

#~~ Any more filtering?

freq <- fread("data/out/all/struc_inds.frq")
gl <- read.PLINK("data/out/all/struc_inds.raw")

# assign population colours

library(RColorBrewer)

cbPalette <- c(brewer.pal(6, "Dark2"))
cbPalette <- wes_palette("Moonrise2")

col <- ifelse(grepl("RR", gl@pop), "SSL",
              ifelse(grepl("CA", gl@pop), "GSL",
                     ifelse(grepl("GS", gl@pop), "GS","AFS")))

unique(col)

# PCA

pca1 <- glPca(gl, returnDotProd=T) # naxes 10
4


# ggplot

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]
ind_names <- gl@ind.names

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  mutate(pop = col) 

# eig

eig <- data.frame(pca1$eig)
eig$percentage = (eig[, 1]/sum(eig$pca1.eig))*100
sum(eig$percentage)
sum(eig$percentage[1:2])

eig$percentage <- round(eig$percentage, digits = 1)
eig$percentage[1]
eig$percentage[2]
eig$percentage[3]

#~~ Plot PCA

png(file="figs/PCA_1_2.png", units = "in", res = 300, height=7, width=7)

ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette,
                     name = "Species") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "right") +
  ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain")) +
  scale_x_reverse()

dev.off()

# PCA 1 & 3

png(file="figs/PCA_1_3.png", units = "in", res = 300, height=7, width=7)

ggplot(ggplot_pca, aes(pc1, pc3)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette,
                     name = "Species") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "right") +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain")) +
  scale_x_reverse()

dev.off()
