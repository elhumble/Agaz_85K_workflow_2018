library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
source("scripts/theme_emily.R")
library(wesanderson)
library(stringr)
library(sequoia)
library(patchwork)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Relatedness            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Prune for LD

# Run PLINK --indep

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--exclude data/processed/x_linked_variants.txt ",
              "--indep 50 5 2 --out data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--allow-extra-chr --debug"))

system("wc -l data/out/agaz/plink/agaz_plate1_3_poly_parents.prune.in")


#~~ Filter for LD, high MAF 0.3, 90% geno for informative loci

# Run PLINK genome and write out for sequoia

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--remove data/out/agaz/plink/dups_to_rm_plink.txt ",
              "--exclude data/processed/x_linked_variants.txt ",
              "--extract data/out/agaz/plink/agaz_plate1_3_poly_parents.prune.in ",
              "--geno 0.1 --maf 0.3 --hwe 0.001 ",
              "--genome ",
              "--out data/out/agaz/plink/sequoia ",
              "--recodeA --allow-extra-chr --debug"))

#~~ Read in genome output

gen <- fread("data/out/agaz/plink/sequoia.genome", header = T)

summary(gen$PI_HAT)


#~~ Specify inference criteria

# https://academic.oup.com/bioinformatics/article/26/22/2867/228512

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1-(1/(2^(3/2))) ~ "Second-degree",
                              kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1 -(1/2^(5/2)) ~ "Third-degree",
                              kinship < 1/2^(9/2) & Z0 > 1-(1/2^(5/2)) ~ "Unrelated",
                              TRUE ~ "Unknown"))


#~~ Combine with ngsrelate output

ngsrel <- fread("data/processed/vcf_maf.res")

# individual order is the same
gen$R1 <- ngsrel$R1 
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING


#~~ Plotting

# Z scores

pal <- c(wes_palette("Darjeeling1"))

ggplot(gen, aes(Z0, Z1)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=1)") +
  scale_colour_manual(values = pal) +
  theme_emily() +
  
  ggplot(gen, aes(Z0, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=2)") +
  ylim(c(0,1)) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  
  ggplot(gen, aes(Z1, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=1)") + ylab("Pr(IBD=2)") +
  ylim(c(0,1)) +
  scale_colour_manual(values = pal) +
  theme_emily()


# R1, R0, KING

pal <- c(wes_palette("Darjeeling1")[2],
         wes_palette("GrandBudapest1")[4],
         wes_palette("Darjeeling2")[4],
         wes_palette("Darjeeling2")[3],
         wes_palette("Darjeeling2")[2])


a <- ggplot(gen, aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.position = "none")
  
b <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.9)) +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  ggtitle("B")



# PI_HAT
 

c <- ggplot(gen, aes(x=PI_HAT)) + 
  geom_histogram(col = pal[5], alpha = 0.9, fill = pal[5], binwidth = 0.01) + 
  scale_y_continuous(trans='log10') +
  xlab("Pairwise relatedness") + ylab("Count") +
  theme_emily() +
  ggtitle("A")



# Figure 5

c + b

ggsave("figs/Figure_5.tiff", c + b, width = 21, height = 10, units = "cm")  
ggsave("figs/Figure_5.png", c + b, width = 21, height = 10, units = "cm")  

# Figure S3


d <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.position = "none") +
  ggtitle("A")

e <- ggplot(gen, aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  ggtitle("B")


d + e

ggsave("figs/Figure_S3.tiff", d + e, width = 21, height = 10, units = "cm")  
ggsave("figs/Figure_S3.png", d + e, width = 21, height = 10, units = "cm")  

# Remove pairs of individuals close to boundaries

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria_2 = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) + 0.01 & kinship < 1/2^(5/2) + 0.01 & Z0 > 0.365 + 0.01 & Z0 < 1-(1/(2^(3/2))) + 0.01 ~ "Second-degree",
                              kinship >= 1/2^(9/2) + 0.01 & kinship < 1/2^(7/2) + 0.01 & Z0 > 1-(1/2^(3/2)) + 0.01 & Z0 < 1 -(1/2^(5/2)) + 0.01 ~ "Third-degree",
                              kinship < 1/2^(9/2) + 0.01 & Z0 > 1-(1/2^(5/2)) + 0.01 ~ "Unrelated",
                              TRUE ~ "Unknown"))

f <- ggplot(filter(gen, criteria_2 != "Unknown"), aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria_2, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.position = "none") +
  ggtitle("A")

g <- ggplot(filter(gen, criteria_2 != "Unknown"), aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = factor(criteria_2, 
                                 levels = c("Parent-offspring",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  ggtitle("B")


f + g

ggsave("figs/Figure_S4.tiff", f + g, width = 21, height = 10, units = "cm")  
ggsave("figs/Figure_S4.png", f + g, width = 21, height = 10, units = "cm")  

#~~~~~~~~~~~~~~~~~~~~~~~#
#        Numbers        #
#~~~~~~~~~~~~~~~~~~~~~~~#

# Proportion of pairs with relatedness UN
nrow(filter(gen, criteria == "Unrelated")) / nrow(gen)
nrow(filter(gen, criteria == "Unrelated"))


# Proportion of third degree relatives i.e. cousins
nrow(filter(gen, criteria == "Third")) / nrow(gen)
nrow(filter(gen, criteria == "Third"))
nrow(filter(gen, criteria_2 == "Third"))

# Proportion of second degree relatives i.e. half sib / avu / grandparent-grandchild
nrow(filter(gen, criteria == "Second")) / nrow(gen)
nrow(filter(gen, criteria == "Second"))
nrow(filter(gen, criteria_2 == "Second"))

# Proportion of PO pairs
nrow(filter(gen, criteria == "PO")) / nrow(gen)
nrow(filter(gen, criteria == "PO"))

# Proportion of uncategorized
nrow(filter(gen, criteria == "OTHER")) / nrow(gen)
nrow(filter(gen, criteria == "OTHER"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Add pedigree information      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in mum pup pair ped file

ped <- fread("data/out/agaz/plink/new_fam_for_plink.txt") %>%
  select(V1, V4, V3) %>%
  `colnames<-`(c("IID", "MOTHER", "FATHER"))

#~~ Join with genome ibd out to compare

# Flip mum and pup in alphabetical order

ped_flip <- ped %>%
  mutate(counter = 1:nrow(.))

ids <- ped_flip[,c(1,2)]
ids <- data.frame(t(apply(ids, 1, str_sort))) %>%
  mutate(counter = 1:nrow(.))

ped_flip <- left_join(ped_flip, ids, by = "counter") %>%
  unite(pair, X1, X2) %>%
  mutate(type = case_when(!is.na(MOTHER) ~ "MO",
                          TRUE ~ "UN"))

#~~ Flip ID1 and ID2 in alphabetical order for PLINK genome output

gen_flip <- gen %>%
  mutate(counter = 1:nrow(.))

ids <- gen_flip[,c(1,3)]
ids <- data.frame(t(apply(ids, 1, str_sort))) %>%
  mutate(counter = 1:nrow(.))

gen_flip <- left_join(gen_flip, ids, by = "counter") %>%
  unite(pair, X1, X2)


#~~ Join known ped with PLINK relationships

check <- gen_flip %>%
  left_join(ped_flip, by = "pair") %>%
  mutate(type = case_when(is.na(type) ~ "UN",
                          TRUE ~ type))

#~~ incorrect relationships based on PI_HAT

filter(check, PI_HAT < 0.4 & !is.na(MOTHER))  # shouldn't be mum-pup pairs (4)
filter(check, PI_HAT > 0.4 & is.na(MOTHER))  # should be mum-pup pairs (3)

#~~ incorrect relationships based on freq

filter(check, criteria != "PO" & !is.na(MOTHER))  # shouldn't be mum-pup pairs (4)
filter(check, criteria == "PO" & is.na(MOTHER))  # should be mum-pup pairs (3)

# 53 - 5 + 4 = 52 mum pup pairs


#~~ Load Sequoia out

load("data/processed/SeqOUT.Rdata")
SummarySeq(SeqOUT)

# Compare assigned parents to field / microsat pedigree

PC.par <- PedCompare(Ped1 = ped[, c("IID", "MOTHER", "FATHER")],
                     Ped2 = SeqOUT$PedigreePar)
PC.par$Counts["TT",,]

PC.par$P1only # parents in provided pedigree that were not assigned in sequoia i.e. shouldn't mum-pup pairs
PC.par$P2only # errors (correct parent not present) should be mother offspring pairs

PC.par$Mismatch # errors (correct parent present)

# half sibs assigned my sequoia

load("data/processed/MaybeRel.Rdata")
Maybe$MaybeRel


#~~ Join sequoia pedigree with PLINK genome out

sequoia_pedigree <- SeqOUT$Pedigree

sequoia_ped_flip <- sequoia_pedigree %>%
  mutate(counter = 1:nrow(.))

ids <- sequoia_ped_flip[c(1,2)]
ids <- data.frame(t(apply(ids, 1, str_sort))) %>%
  mutate(counter = 1:nrow(.))

sequoia_ped_flip <- left_join(sequoia_ped_flip, ids, by = "counter") %>%
  unite(pair, X1, X2) %>%  
  mutate(seq_type = case_when(!is.na(dam) ~ "MO",
                              TRUE ~ "UN"))

df <- check %>%
  left_join(sequoia_ped_flip, by = "pair") %>%
  mutate(seq_type = case_when(is.na(seq_type) ~ "UN",
                              TRUE ~ seq_type))

#~~ Highlight MO pairs assigned by sequoia

ggplot(df, aes(Z0, Z1)) + 
  geom_point(aes(colour = factor(seq_type)), alpha = 0.6, size = 4) +
  scale_color_manual(values = pal[c(2,1)]) +
  theme_emily()

ggplot(df, aes(Z1, Z2)) + 
  geom_point(aes(colour = factor(seq_type)), alpha = 0.6, size = 3) +
  scale_color_manual(values = pal[c(2,1)]) +
  theme_emily()


#~~ Highlight half sibs assigned by sequoia

M0001 <- filter(sequoia_pedigree, sire == "M0001")
M0001 <- as.data.frame(t(combn(M0001$id, 2)))

M0002 <- filter(sequoia_pedigree, sire == "M0002")
M0002 <- as.data.frame(t(combn(M0002$id, 2)))

flip_rows <- function(df) {
  
  df_flip <- df %>%
    mutate(counter = 1:nrow(.))
  
  ids <- df_flip[c(1,2)]
  ids <- data.frame(t(apply(ids, 1, str_sort))) %>%
    mutate(counter = 1:nrow(.))
  
  df_flip <- left_join(df_flip, ids, by = "counter") %>%
    unite(pair, X1, X2)
  
  df_flip <- df_flip
  
}

M0001_flip <- flip_rows(M0001) %>%
  mutate(HS = T)

M0002_flip <- flip_rows(M0002) %>%
  mutate(HS = T)

HS <- rbind(M0001_flip, M0002_flip)

df_HS <- df %>%
  left_join(HS, by = "pair") %>%
  mutate(seq_type = case_when(!is.na(HS) ~ "HS",
                              !is.na(dam) ~ "MO",
                      TRUE ~ "UN"))

ch <- wes_palette("FantasticFox1")
em <- c(pal[c(2,3)],ch[3])


ggplot(df_HS, aes(Z0, Z1)) +
  geom_jitter(width = 0.1, aes(colour = factor(seq_type)), alpha = 0.6, size = 3) +
  scale_color_manual(breaks = c("HS", "MO", "UN"),
                     labels = c("Half-sibs", "Mother-offsping", "Unknown"),
                     values = em[c(3,1,2)]) +
  theme_emily()


ggplot(df_HS, aes(Z0, Z1)) +
  geom_jitter(width = 0.1, aes(colour = factor(type)), alpha = 0.6, size = 3) +
  theme_emily()


#~~~~~~~~~~~~~~~~~~~~#
#       Plot         #
#~~~~~~~~~~~~~~~~~~~~#

a <- ggplot(gen, aes(x=PI_HAT)) + 
  geom_histogram(col = pal[2], alpha = 0.9, fill = pal[2], binwidth = 0.01) + 
  scale_y_continuous(trans='log10') +
  xlab("Pairwise relatedness") + ylab("Count") +
  theme_emily() +
  ggtitle("A")

b <- ggplot(df_HS, aes(Z0, Z1)) +
  geom_jitter(width = 0.05, aes(colour = factor(seq_type)), alpha = 0.6, size = 3) +
  scale_color_manual(breaks = c("HS", "MO", "UN"),
                     labels = c("Half-sibs", "Mother-offsping", "Unknown"),
                     values = em[c(3,1,2)],
                     name = "") +
  theme_emily() +
  theme(legend.position = c(0.8, 0.9)) +
  ggtitle("B")

# without unknowns in M-O

df_HS <- df_HS %>%
  mutate(legend = case_when(Z1 > 0.75 & seq_type == "UN" ~ "MO", TRUE ~ seq_type))

b <- ggplot(df_HS, aes(Z0, Z1)) +
  geom_jitter(width = 0.05, aes(colour = factor(legend)), alpha = 0.6, size = 3) +
  scale_color_manual(breaks = c("HS", "MO", "UN"),
                     labels = c("Half-sibs", "Mother-offsping", "Unknown"),
                     values = em[c(3,1,2)],
                     name = "") +
  theme_emily() +
  theme(legend.position = c(0.8, 0.9)) +
  ggtitle("B")


png(file="figs/Figure_4.png", units = "in", res = 300, height=4, width=9)
grid.arrange(a, b, ncol = 2)
dev.off()

png(file="figs/Figure_4.png", units = "in", res = 300, width = 9, height = 4)
plot_grid(a, b, ncol = 2, align="v")
dev.off()


#~~ Write list of avuncular relationships

avu <- df_HS %>%
  filter(Z1 >= 0.1 & Z1 <= 0.75 & 
         Z0 >= 0.25 & Z0 <= 0.8) %>%
  select(pair, Z0, Z1, Z2, PI_HAT) %>%
  separate(pair, c("X1", "X2"))


write.csv(avu[c(1:6)], "data/out/agaz/plink/cryptic_relatives.csv",
          quote = F, row.names = F)

#~~ Age difference

LifeHist <- fread("data/processed/lifehist.txt")

age_diff <- avu %>%
  left_join(LifeHist, by = c("X1" = "V1")) %>%
  left_join(LifeHist, by = c("X2" = "V1")) %>%
  mutate(agediff = abs(V3.x - V3.y))

hist(age_diff$agediff)

png("figs/age_diff.png", units = "in", res = 300, width = 6, height = 5)

ggplot(age_diff, aes(x=agediff)) +
  geom_histogram(alpha=0.6, binwidth = 2.5) +
  labs(y = "Count", x = "Age difference") +
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_blank()) +
  theme_emily()

dev.off()

