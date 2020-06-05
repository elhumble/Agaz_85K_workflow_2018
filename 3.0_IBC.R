library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
source("scripts/theme_emily.R")
library(wesanderson)
library(inbreedR)
library(patchwork)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Inbreeding analyis      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Prepare plink file

# Filter for MAF, geno and hwe
# Remove duplicate individuals
# Remove X linked SNPs
# Write to raw out for sMLH

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--exclude data/processed/x_linked_variants.txt ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 --recodeAD --make-bed ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding ",
              "--allow-extra-chr --debug"))


#~~ Function to get sMLH values from raw plink file

get_sMLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  #row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  sMLH <- as.data.frame(sMLH(x))
  sMLH$ANIMAL <- ids
  sMLH$NAS <- NAs
  colnames(sMLH) <- c("sMLH", "ANIMAL", "NAs")
  sMLH <- sMLH
  
}

# Run sMLH on PLINK raw

raw_file <- "data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding.raw"
sMLH <- get_sMLH_from_plinkraw(raw_file)

ggplot(sMLH, aes(sMLH)) +
  geom_histogram()

ggplot(sMLH, aes(NAs)) +
  geom_histogram()

pal <- wes_palette("Darjeeling2")

sMLH_dist <- ggplot(sMLH, aes(sMLH)) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = pal[2], fill = pal[2]) +
  geom_density(alpha=0.6, col = pal[2], fill = pal[2]) + #3262AB
  labs(y = "Density", x = "sMLH") +
  theme_emily() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = pal[2]) +
  ggtitle("A")

sMLH_dist

#~~ Get Fhats and ROH

# recode bim files for GCTA

recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}


bim_file <- "data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding.bim"
recode_bim_1chr(bim_file)

# Run GCTA

plink_files <- gsub(".bim", "", bim_file)

for (i in 1:length(plink_files)){
  system(paste0("~/software/gcta_1.92.4beta_mac/bin/gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 4"))
}

# Load fhats

fhats <- fread("data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding.ibc", header = T)

ibcs <- fhats %>%
  left_join(sMLH, by = c("IID" = "ANIMAL"))

# Run PLINK homozyg
# Same filters as above

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--exclude data/processed/x_linked_variants.txt ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 ",
              "--homozyg --homozyg-window-snp 20 --homozyg-snp 20 ",
              "--homozyg-kb 1000 --homozyg-gap 1000 --homozyg-density 100 ",
              "--homozyg-window-missing 5 --homozyg-window-het 1 --homozyg-het 1 ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding ",
              "--allow-extra-chr --debug"))

# ROH length dist

roh_dist <- fread("data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding.hom")
summary(roh_dist$KB)
hist(roh_dist$KB)

# Median of total length of short ROH (<5000)
roh_dist %>%
  group_by(IID) %>%
  filter(KB < 5000) %>%
  summarise(sum = sum(KB) / 1000) %>%
  summarise(median = median(sum),
            mean = mean(sum))

# Median of total length of long ROH (>5000)
roh_dist %>%
  group_by(IID) %>%
  filter(KB >= 5000) %>%
  summarise(sum = sum(KB) / 1000) %>%
  summarise(median = median(sum),
            mean = mean(sum))

# Median of total length of all ROH
roh_dist %>%
  group_by(IID) %>%
  summarise(sum = sum(KB) / 1000) %>%
  summarise(median = median(sum),
            mean = mean(sum))

# Max length in all individuals

roh_dist %>%
  group_by(IID) %>%
  summarise(max = max(KB)) %>%
  summarise(min = min(max) / 1000)


roh_lengths_A <- filter(roh_dist, KB < 5000) %>%
  ggplot(aes(KB)) + geom_histogram(alpha=0.9, col = "grey33", fill = "grey33") +
  labs(y = "Density", x = "ROH length (kb)") +
  theme_emily() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = pal[2]) +
  ggtitle("(A) ROH < 5 Mb")

roh_lengths_B <- filter(roh_dist, KB >= 5000) %>%
  ggplot(aes(KB)) + geom_histogram(alpha=0.9, col = "grey33", fill = "grey33") +
  labs(y = "Density", x = "ROH length (kb)") +
  theme_emily() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = pal[2]) +
  ggtitle(expression("(B) ROH" >= "5 Mb"))

#~~ Individual values

roh <- fread("data/out/agaz/plink/agaz_plate1_3_poly_parents_inbreeding.hom.indiv")
hist(roh$KB/2300000)

# Summary number of ROH segments per individual
summary(roh$NSEG)

ibcs <- ibcs %>%
  left_join(roh, by = "IID") %>%
  dplyr::mutate(froh = KB/2300000)

summary(ibcs$froh)
hist(ibcs$froh)

roh_plot <- ggplot(ibcs, aes(froh, fill = as.factor(Run))) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "grey33", fill = "grey33") +
 # geom_density(alpha=0.6, col = "grey33", fill = "grey33") + #3262AB
  labs(y = "Density", x = expression(italic(F["ROH"]))) +
  theme_emily() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = pal[2]) +
  ggtitle("A")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Correlation plots     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ibcs %>%
  select(Fhat1, Fhat2, Fhat3, sMLH, froh) %>%
  pairs()


corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}


#~~ Correlation coefficients

corr_eqn(ibcs$sMLH, ibcs$Fhat3)
corr_eqn(ibcs$sMLH, ibcs$froh)
corr_eqn(ibcs$Fhat3, ibcs$froh)


#~~ Plots

labels3 <- data.frame(x = 1.026, y = 0.03, label = corr_eqn(ibcs$sMLH, ibcs$Fhat3))

fhat3_sMLH <- ggplot(ibcs, aes(x=sMLH, y = Fhat3)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") + #3262AB
  labs(x = "sMLH", y = expression(italic(hat(F)["III"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  # geom_text(aes(x = -0.03, y = 0.15, label = lm_eqn(m1)), parse = TRUE,
  #           size = 5, fontface = "plain") +
  geom_text(data = labels3, aes(x = x, y = y, label = label), parse = TRUE,
            size = 4, fontface = "plain") +
  #ggtitle('(C)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain")) +
  theme_emily() +
  ggtitle("B")


labels4 <- data.frame(x = 0.035, y = 0.028, label = corr_eqn(ibcs$froh, ibcs$Fhat3))

fhat3_roh <- ggplot(ibcs, aes(x=froh, y = Fhat3)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") + #3262AB
  labs(x = expression(italic(F["ROH"])), y = expression(italic(hat(F)["III"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  # geom_text(aes(x = -0.03, y = 0.15, label = lm_eqn(m1)), parse = TRUE,
  #           size = 5, fontface = "plain") +
  geom_text(data = labels4, aes(x = x, y = y, label = label), parse = TRUE,
            size = 4, fontface = "plain") +
  # ggtitle('(C)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain")) +
  theme_emily() +
  ggtitle("C")

labels5 <- data.frame(x = 0.078, y = 1.035, label = corr_eqn(ibcs$froh, ibcs$sMLH))

sMLH_roh <- ggplot(ibcs, aes(x=froh, y = sMLH)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") + #3262AB
  labs(x = expression(italic(F["ROH"])), y = "sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  # geom_text(aes(x = -0.03, y = 0.15, label = lm_eqn(m1)), parse = TRUE,
  #           size = 5, fontface = "plain") +
  geom_text(data = labels5, aes(x = x, y = y, label = label), parse = TRUE,
            size = 4, fontface = "plain") +
  # ggtitle('(C)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain")) +
  theme_emily() +
  ggtitle("D")


# Correlation figure

fig_3 <- roh_plot + fhat3_sMLH + fhat3_roh + sMLH_roh
fig_4 <- roh_lengths_A + roh_lengths_B

ggsave("figs/Figure_3.tiff", fig_3, width = 19, height = 16, units = "cm")
ggsave("figs/Figure_4.tiff", fig_4, width = 18, height = 10, units = "cm")

 
#~~ Other distributions

ggplot(ibcs, aes(froh)) +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = pal[2], fill = pal[2]) +
  geom_density(alpha=0.6, col = pal[2], fill = pal[2]) + #3262AB
  labs(y = "Density", x = "sMLH") +
  theme_emily() +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = pal[2])


#~~~~~~~~~~~~~~~~~~~~~#
#         g2          #
#~~~~~~~~~~~~~~~~~~~~~#

get_g2_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "numeric")
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  
  g2 <- g2_snps(x, nperm = 1000, nboot = 1000, CI = 0.95)
  g2
  g2 <- g2
  
}

# run once
# g2 <- get_g2_from_plinkraw(raw_file)
# saveRDS(g2, file = "data/out/agaz/plink/g2.rds")

# Read output from RDS 

g2 <- readRDS("data/out/agaz/plink/g2.rds")
plot(g2, col = "grey")


g2_plot <- data.frame(g2$g2_boot)
lcl <- g2$CI_boot[1]
ucl <- g2$CI_boot[2]
g2_boot_summary <- data.frame(lcl, ucl)

ibcs_vars_summary <- data.frame(c(g2$g2,
                                  var(ibcs$sMLH, na.rm = T),
                                  var(ibcs$Fhat1),
                                  var(ibcs$Fhat2),
                                  var(ibcs$Fhat3),
                                  var(ibcs$froh)),
                                c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3", "froh"))


colnames(ibcs_vars_summary) <- c("val", "var")

# g2 bootstrapping distribution showing empirical g2 with CIs

cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3", "#D95F02", "#E7298A")

pal <- c("black", wes_palette("Darjeeling1"))

#g2_CI_plot <- 
ggplot(g2_plot, aes(g2.g2_boot)) + 
  geom_histogram(colour = "grey45", fill = "grey45") +
  geom_errorbarh(aes(xmin = g2_boot_summary$lcl , xmax = g2_boot_summary$ucl , y = 165),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  geom_vline(data = ibcs_vars_summary, aes(xintercept = val, colour = var), size = 0.8, 
             linetype = c("dashed", "solid", "solid", "solid", "solid", "solid"), show.legend = T) +
  scale_colour_manual(values = pal[c(2,3,4,5,1,6)], name = "",
                      breaks = c("g2","Fhat1", "Fhat2", "Fhat3", "sMLH", "froh"),
                      labels = c(expression(italic(g[2])), 
                                 expression("var"(italic(hat(F)["I"]))), 
                                 expression("var"(italic(hat(F)["II"]))), 
                                 expression("var"(italic(hat(F)["III"]))), 
                                 expression("var"("sMLH")),
                                 expression(italic("F"["ROH"])))) +
  labs(y = "Counts", x = expression(italic(g[2]))) +
  theme_emily()
