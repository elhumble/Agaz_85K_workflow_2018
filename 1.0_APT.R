# Script to execute Best Practices Steps with APT

library(tidyverse)
library(readxl)
library(data.table)

#~~ Create list of CEL files (Run in batches = plates?) 1.

cel_list <- paste("data/raw/", (list.files("data/raw", pattern = "*.CEL")), sep = "")
write.table(cel_list, "data/cel.list", quote = F, col.names = "cel_files", row.names = F)

#~~ Best Practices step 2 : Sample Dish QC values

system("apt-geno-qc --cdf-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.cdf --qca-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.qca --qcc-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K.r1.qcc --cel-files data/cel.list --out-file data/out/seal.CQC.txt")

# Explore out

dqc <- fread("data/out/seal.CQC.txt")
plot(dqc$axiom_dishqc_DQC)

#~~ Best Practices step 3 : filter based on DQC

# filter samples below 0.82 threshold

# Remove from cel list (not necessary as all over threshold) = cel_list2.txt

#~~ Best Practices step 4 : Generate sample QC call rates

system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-probeset-genotype --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1 --cel-files data/cel.list --out-dir data/out --xml-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.test.xml --no-gender-force")

# alt
# system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir test_out/ --cel-files data/cel.list --log-file test.log")



#~~ Best Practices Step 5 : Filter based on call rates


call_rate <- read.table("data/out/AxiomGT1.report.txt", header = T)
hist(call_rate$call_rate)


# again not necessary -- all above call rate


#~~ Best Practices Step 7 : Genotype passing samples (and plates - not run) using AxiomGT1.Step2

system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-probeset-genotype --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1 --cel-files data/cel.list --out-dir data/out --xml-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.test.xml --no-gender-force --write-models --summaries")

# alt
# system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-genotype-axiom --analysis-files-path data/library_files/Axiom_Agaz_90K_Analysis.r1/ --arg-file data/library_files/Axiom_Agaz_90K_Analysis.r1/Axiom_Agaz_90K_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml --out-dir test_out/ --cel-files data/cel.list --log-file test.log")


#~~ Best Practices Step 8 : Get metrics and classifications

# Get SNP metrics
system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/ps-metrics --posterior-file data/out/AxiomGT1.snp-posteriors.txt --call-file data/out/AxiomGT1.calls.txt --metrics-file data/out/metrics.txt")
       
system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/ps-classification --species-type diploid --metrics-file data/out/metrics.txt --output-dir data/out/")

# Get numbers

system("wc -l data/out/ps_classification/*")

  
could run otv caller on snps classified into otv category by ps-classification
  

#~~ Generate PLINK files

system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file data/out/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file data/out/plate_1")

# alt

system("~/programs/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file test_out/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file test_out/plate_1")


# Fix map file (this bit should be in another script)

library(data.table)
library(dplyr)
library(tidyr)

map <- fread("test_out/plate_1.map")

# Read in annotation file

annot <- read.csv("data/library_files/Axiom_Agaz_90K_Annotation.r1.csv", skip = 18) %>%
  select(Probe.Set.ID, cust_id)

# Join by Marker ID and Probe.Set.ID
# _v1.1
# JIH_SNP_ID


temp <- left_join(map, annot, by = c("Marker ID" = "Probe.Set.ID")) %>%
  mutate(cust_id_orig = cust_id) %>%
  mutate(cust_id = gsub("_v1.1", "v1.1", cust_id)) %>%
  mutate(cust_id = gsub("JIH_AG_SNP", "JIH.AG.SNP", cust_id)) %>%
  mutate(cust_id = gsub("_rs", "_", cust_id)) %>%
  separate(cust_id, c("Chr", "Position"), sep = "_", remove = F)

# re-write map file : triple check. This MUST be same order

# chr, id, dist, pos

new_map <- data.frame(temp$Chr, temp$cust_id_orig, temp$`Genetic distance`, temp$Position)
colnames(new_map) <- c("Chr", "ID", "Distance", "Position")
write.table(new_map, "test_out/plate_1_chr.map", 
            quote = F, row.names = F, col.names = F)

# filter list of RAD SNPs

rad <- write.table(new_map[1:2], "test_out/rad_snps.txt",
                   quote = F, row.names = F, col.names = F)

#~~ Filter plink file and run LD decay


