# Script to execute Best Practices Steps with APT

library(dplyr)
library(readxl)
library(data.table)
library(tidyr)
library(purrr)

#~~ Create list of CEL files

plate1 <- paste("data/raw/", (list.files("data/raw", pattern = "*plate1.CEL")), sep = "")
plate2 <- paste("data/raw/", (list.files("data/raw", pattern = "*plate2.CEL")), sep = "")
plate3 <- paste("data/raw/", (list.files("data/raw", pattern = "*plate3.CEL")), sep = "")

cel_list <- c(plate1, plate2, plate3)

write.table(cel_list, "data/processed/cel.list", quote = F, col.names = "cel_files", row.names = F)

# Analyse plates in batches corresponding to when they were run.

#~~ Best Practices step 2 : Sample Dish QC values - Quality metric for each sample

system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-geno-qc --cdf-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.cdf --qca-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.qca --qcc-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.qcc --cel-files data/processed/cel.list --out-file data/out/seal.DQC.txt --log-file data/out/apt-geno-qc.log")

# Explore output

dqc <- fread("data/out/seal.DQC.txt")
plot(dqc$axiom_dishqc_DQC)

#~~ Best Practices step 3 : filter based on DQC (<0.82)

# Which individuals are below DQC threshold?
# Mostly additional species

dqc %>%
  filter(axiom_dishqc_DQC < 0.82) %>%
  select(cel_files)

# Filter

dqc_filt <- dqc %>%
  mutate(cel_files = paste("data/raw/", cel_files, sep = "")) %>%
  filter(axiom_dishqc_DQC > 0.82)

write.table(dqc_filt[1], "data/processed/cel.list2", quote = F, col.names = "cel_files", row.names = F)


#~~ Best Practices step 4 : Generate sample QC call rates

system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir data/out/ --cel-files data/processed/cel.list2 --log-file test.log")

# Do the same for whole dataset

system("mkdir data/out/call_rates_all")
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir data/out/call_rates_all --cel-files data/processed/cel.list --log-file test.log")

#~~ Best Practices Step 5 : Filter based on call rates (~97%) and generate new cel file

call_rate <- fread("data/out/AxiomGT1.report.txt", skip = 365) # hard coded
hist(call_rate$call_rate)

call_rate %>%
  filter(call_rate < 97) %>%
  select(cel_files)

# Most of these individuals are other species.

# Remove samples below 97% call rate plus additional species

call_rate_filt_sp <- call_rate %>%
  separate(cel_files, c("well", "sampleID"), sep = "-") %>%
  filter(call_rate > 97) %>%
  filter(grepl("^A", sampleID)) %>%
  unite(cel_files, c("well", "sampleID"), sep = "-") %>%
  mutate(cel_files = paste("data/raw/", cel_files, sep = ""))

write.table(call_rate_filt_sp[1], "data/processed/cel.list3", quote = F, col.names = "cel_files", row.names = F)


#~~ Best Practices Step 6 : Plate QC (~98.5%)

qc_plate <- function(dqc_filt, call_rate) {
  
  a <- select(dqc_filt, cel_files, axiom_dishqc_DQC) %>%
    mutate(cel_files = gsub("data/raw/", "", cel_files))
  b <- select(call_rate, cel_files, call_rate)
  
  df <- left_join(a, b) %>%
    mutate(cel_files = gsub("\\(Axiom\\_AGaz90k\\)\\_", "", cel_files)) %>%
    separate(cel_files, c("sample", "plate"), sep = "_plate") %>%
    separate(plate, c("plate", "file")) %>%
    select(-file)
  
  df %>%
    filter(axiom_dishqc_DQC > 0.82 & call_rate > 97) %>%
    group_by(plate) %>%
    summarise(n() / 96)
  
}


# Could remove very bad plates. Plate QC is lower when additional species are present.
# Proceed for now.

#~~ Best Practices Step 7 : Genotype passing samples using AxiomGT1.Step2

system("mkdir data/out/agaz")
system("mkdir data/out/all")

# Plate including additional species (all)

system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir data/out/all --cel-files data/processed/cel.list2 --log-file data/out/all/apt-genotype-axiom2.log --write-models --summaries")

# Plate without additional species (agaz)
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir data/out/agaz --cel-files data/processed/cel.list3 --log-file data/out/agaz/apt-genotype-axiom2.log --write-models --summaries")

#~~ Best Practices Step 8 : Get metrics and classifications

# Get SNP metrics

# all
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-metrics --posterior-file data/out/all/AxiomGT1.snp-posteriors.txt --call-file data/out/all/AxiomGT1.calls.txt --metrics-file data/out/all/metrics.txt --log-file data/out/all/ps-metrics.log")
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-classification --species-type diploid --metrics-file data/out/all/metrics.txt --output-dir data/out/all --log-file data/out/all/ps-classification")

# agaz
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-metrics --posterior-file data/out/agaz/AxiomGT1.snp-posteriors.txt --call-file data/out/agaz/AxiomGT1.calls.txt --metrics-file data/out/agaz/metrics.txt --log-file data/out/agaz/ps-metrics.log")
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-classification --species-type diploid --metrics-file data/out/agaz/metrics.txt --output-dir data/out/agaz/ --log-file data/out/agaz/ps-classification")

# Get numbers (85359 SNPs on array)

system("wc -l data/out/all/*.ps")
system("wc -l data/out/agaz/*.ps")

# Hybridised:
# No minor hom
# Poly High Res
# Mono High Res

# Polymorphic
# No minor hom
# Poly High Res


# Could run OTV caller on snps classified into otv category by ps-classification

install.packages("~/software/SNPolisher_2.0.tar.gz", repos=NULL, type = "source")
library(SNPolisher)
system("mkdir data/out/all/OTV")
system("mkdir data/out/agaz/OTV")

#~~ All

SNPolisher::OTV_Caller(summaryFile = "data/out/all/AxiomGT1.summary.txt", 
                       posteriorFile = "data/out/all/AxiomGT1.snp-posteriors.txt",
                       callFile = "data/out/all/AxiomGT1.calls.txt", 
                       confidenceFile = "data/out/all/AxiomGT1.confidences.txt",
                       pidFile = "data/out/all/OffTargetVariant.ps",
                       output.dir = "data/out/all/OTV", OTV.only = T)


system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-metrics --posterior-file data/out/all/OTV/OTV.snp-posteriors.txt --call-file data/out/all/OTV/OTV.calls.txt --metrics-file data/out/all/OTV/ps-metrics.txt --log-file data/out/all/OTV/ps-metrics.log")
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-classification --species-type diploid --metrics-file data/out/all/OTV/ps-metrics.txt --output-dir data/out/all/OTV --log-file data/out/all/OTV/ps-classification")

#~~ Agaz


SNPolisher::OTV_Caller(summaryFile = "data/out/agaz/AxiomGT1.summary.txt", 
                       posteriorFile = "data/out/agaz/AxiomGT1.snp-posteriors.txt",
                       callFile = "data/out/agaz/AxiomGT1.calls.txt", 
                       confidenceFile = "data/out/agaz/AxiomGT1.confidences.txt",
                       pidFile = "data/out/agaz/OffTargetVariant.ps",
                       output.dir = "data/out/agaz/OTV", OTV.only = T)


system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-metrics --posterior-file data/out/agaz/OTV/OTV.snp-posteriors.txt --call-file data/out/agaz/OTV/OTV.calls.txt --metrics-file data/out/agaz/OTV/ps-metrics.txt --log-file data/out/agaz/OTV/ps-metrics.log")
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/ps-classification --species-type diploid --metrics-file data/out/agaz/OTV/ps-metrics.txt --output-dir data/out/agaz/OTV --log-file data/out/agaz/OTV/ps-classification")


#~~ Generate PLINK files

system("mkdir data/out/all/plink")
system("mkdir data/out/agaz/plink")

# all
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file data/out/all/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file data/out/all/plink/agaz_plate1_3 --log-file data/out/all/plink/plink.log")

# agaz
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file data/out/agaz/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file data/out/agaz/plink/agaz_plate1_3 --log-file data/out/agaz/plink/plink.log")


#~~ Format map files

# Join by Marker ID and Probe.Set.ID
# Fix _v1.1
# FIx JIH_SNP_ID

fix_map <- function(map, annot) {
  
  temp <- left_join(map, annot, by = c("Marker ID" = "Probe.Set.ID")) %>%
    mutate(cust_id_orig = cust_id) %>%
    mutate(cust_id = gsub("_v1.1", "v1.1", cust_id)) %>%
    mutate(cust_id = gsub("JIH_AG_SNP", "JIH.AG.SNP", cust_id)) %>%
    mutate(cust_id = gsub("_rs", "_", cust_id)) %>%
    separate(cust_id, c("Chr", "Position"), sep = "_", remove = F) %>%
    mutate(Chr = gsub("v1.1", "_v1.1", Chr),
           Chr = gsub("JIH.AG.SNP", "JIH_AG_SNP", Chr),
           cust_id_orig = gsub("_rs", "_", cust_id_orig))
  
  # re-write map file : triple check. This MUST be same order
  
  # chr, id, dist, pos
  # note that ids have changed from _ to . etc.
  
  new_map <- data.frame(temp$Chr, temp$cust_id_orig, temp$`Genetic distance`, temp$Position)
  colnames(new_map) <- c("Chr", "ID", "Distance", "Position")
  
  # recode NAs to 0s for SNPs with no known pos
  
  new_map <- new_map %>%
    mutate(Position = as.numeric(as.character(Position)),
           Position = replace_na(Position, 0))
  new_map <- new_map
  
}

annot <- read.csv("data/library_files/Axiom_Agaz_90K_Annotation.r1.csv", skip = 18) %>%
  select(Probe.Set.ID, cust_id)

#~~ All

map_all <- fread("data/out/all/plink/agaz_plate1_3.map")

write.table(fix_map(map_all, annot), "data/out/all/plink/agaz_plate1_3.map", 
            quote = F, row.names = F, col.names = F)

#~~ Agaz

map_agaz <- fread("data/out/agaz/plink/agaz_plate1_3.map")

write.table(fix_map(map_agaz, annot), "data/out/agaz/plink/agaz_plate1_3.map", 
            quote = F, row.names = F, col.names = F)

#~~ Recode IDs in ped file

# Bash commands:

# Agaz

# awk 'gsub("\\_\\(Axiom\\_AGaz90k\\)\\_plate.+\.CEL", "", $1)1' agaz_plate1_3.ped > test
# awk 'gsub("^.+\\-", "", $1)1' test > agaz_plate1_3.ped
# rm test

# All: Not run

# awk 'gsub("\\_\\(Axiom\\_AGaz90k\\)\\_plate.+\.CEL", "", $1)1' data/out/all/plink/agaz_plate1_3.ped > test
# awk 'gsub("^.+\\-", "", $1)1' test > data/out/all/plink/agaz_plate1_3.ped
# rm test

# Manually edited dups

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Pedigree info        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ PED info

seals <- read.csv("data/SSBseals_wrangled_Dec_2017.csv") %>%
  select(PUP, MOTHER, MOTHER_ID1) %>%
  mutate(PUP = as.character(PUP),
         MOTHER = as.character(MOTHER),
         MOTHER_ID1 = as.character(MOTHER_ID1))

AT <- data.frame(PUP = "ATP16015", 
                 MOTHER = "ATF16015",
                 MOTHER_ID1 = NA)

seals <- rbind(seals, AT) %>%
  mutate(MOTHER = case_when(MOTHER_ID1 == "AGF99049" ~ "AGF99049",
                            TRUE ~ MOTHER)) %>%
  select(-MOTHER_ID1)


ped_file <- read.table("data/out/agaz/plink/agaz_plate1_3.ped", skip  = 1)

ped <- as.data.frame(ped_file[,1])
colnames(ped) <- "IID"

ped_update <- ped

ped <- ped %>%
  left_join(seals, by = c("IID" = "PUP")) %>%
  mutate(MOTHER = as.character(MOTHER)) %>%
  mutate(a = MOTHER %in% IID) %>%
  mutate(MOTHER = case_when(a == T ~ MOTHER,
                            TRUE ~ "NA")) %>%
  mutate(MOTHER = na_if(MOTHER, "NA")) %>%
  mutate(FATHER = NA) %>%
  select(-a)

sum(!is.na(ped$MOTHER))
sum(ped$MOTHER %in% ped$IID)

#~~ Susie's code

ped$MOTHER <- as.character(ped$MOTHER)
ped$FATHER <- as.character(ped$FATHER)

famped <- NULL

for(i in ped$IID){
  ped1 <- ped[which(ped$IID == i),]
  {
    ped2 <- ped1[,1:3]
    ped2 <- rbind(data.frame(IID = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    famped <- rbind(famped, ped2)
    famped <- famped[which(famped[,1] != 0),]
    rm(ped2)
    rm(ped1)
  }
}

famped <- famped %>%
  dplyr::mutate(FID = IID,
                Phenotype = -9,
                Sex = 2) %>%
  dplyr::mutate(FATHER = ifelse(is.na(FATHER),NA,FATHER),
                MOTHER = ifelse(is.na(MOTHER),NA,MOTHER))

famped <- famped %>%
  dplyr::distinct(IID, .keep_all = T)

# Create fam file

# fam

head(famped)

new_famped <- ped_update %>%
  left_join(famped, by = "IID")
head(new_famped)

write.table(famped[c(4,1,2,3,6,5)], "data/out/agaz/plink/agaz_plate1_3.fam", col.names = F,
            row.names = F, quote = F)
