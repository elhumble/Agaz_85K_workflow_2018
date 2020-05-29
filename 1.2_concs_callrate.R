library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
source("theme_emily.R")
library(wesanderson)

# Concentration vs call rate for all samples

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    DNA Concentration vs call rate    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in call rate data from Axiom report

call_rate <- fread("data/out/AxiomGT1.report.txt", header = T, skip = 363) %>%
  separate(cel_files, c("Well", "Sample"), sep = "-") %>%
  separate(Sample, c("Sample", "Chip"), sep = "_\\(A")

# Add code for species
call_rate <- call_rate %>%
  mutate(sp = case_when(grepl("CA", Sample) ~ "CA",
                        grepl("RR", Sample) ~ "RR",
                        grepl("GS", Sample) ~ "GS",
                        TRUE ~ "AG"))


nrow(filter(call_rate, sp == "AG")) # 283 - 8

ggplot(call_rate, aes(sp, call_rate)) +
  geom_boxplot()

# Read in concentrations from sample prep

concs <- read.csv("../SAMPLES/data/BGI/Batch_1_May_2018.csv") %>%
  select(Sample_name, C_ng_microl)

concs <- call_rate %>%
  left_join(concs, by = c("Sample" = "Sample_name")) 

concs <- concs[-c(2:5, 286:289, 291:294),]


hist(concs$C_ng_microl)
nrow(filter(concs, C_ng_microl < 50))
nrow(filter(concs, C_ng_microl < 40))
nrow(filter(concs, C_ng_microl < 50 & C_ng_microl >= 40))
nrow(filter(concs, C_ng_microl >= 50))

pal <- wes_palette("Darjeeling2")

png(file="figs/Call_rate_concentration.png", units = "in", res = 300, height=5, width=7)

ggplot(concs, aes(C_ng_microl, call_rate, fill = factor(sp))) +
  geom_point(size = 3, shape = 21, alpha = 0.8) +
  theme_emily() +
  geom_hline(yintercept = 97, col = "grey30", linetype = "dashed") +
  geom_vline(xintercept = 50, col = "grey30", linetype = "dashed") +
  scale_fill_manual(breaks = c("AG", "RR", "CA"),
                     labels = c("Antarctic fur seal", "Stellar sea lion", "Galapagos sea lion"),
                     values = pal[c(2,3,4)],
                     name = "Species") +
  xlab(expression(paste("DNA Concentration (ng/",mu,"l)"))) + ylab("Call rate")

dev.off()


# Linear regression of conc ~ call rate

lm_call_conc <- lm(concs$C_ng_microl~concs$call_rate)
summary_lm <- summary(lm_call_conc)$coefficients
summary(lm_call_conc)

plot(lm_call_conc, 1)

# slope
summary_lm[2,1]
# t
summary_lm[2,3]
# df
lm_call_conc$df.residual
# P
summary_lm[2,4]


# Remove outliers

hi <- concs %>%
  filter(call_rate > 97)

lm_call_conc_hi <- lm(hi$C_ng_microl~hi$call_rate)
summary_lm_hi <- summary(lm_call_conc_hi)$coefficients
summary(lm_call_conc_hi)


# Fur seals only

afs_concs <- concs %>%
  filter(sp == "AG")

lm_call_conc_afs <- lm(afs_concs$C_ng_microl~afs_concs$call_rate)
summary_lm_afs <- summary(lm_call_conc_afs)$coefficients
summary(lm_call_conc_afs)

plot(lm_call_conc_afs, 1)

