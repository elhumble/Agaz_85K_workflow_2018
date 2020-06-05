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
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file data/out/all/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file data/out/all/plink/agaz_plate1_3_apt --log-file data/out/all/plink/agaz_plate1_3_apt.log")

# agaz
system("~/software/apt-2.10.0-x86_64-apple-yosemite/bin/apt-format-result --calls-file data/out/agaz/AxiomGT1.calls.txt --annotation-file data/library_files/Axiom_Agaz_90K.r1.20180809.annot.db --export-plink-file data/out/agaz/plink/agaz_plate1_3_apt --log-file data/out/agaz/plink/agaz_plate1_3_apt.log")

