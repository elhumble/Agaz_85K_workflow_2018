library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(data.table)

#~ Prep flanking sequences for mapping to X chromosome

affy_meta <- fread("data/processed/axiomdesign_Seq_FurSeal_Prioritised_Consideration_01Feb2018") %>%
  select(Snpid, Seventyonemer) %>%
  mutate(Seventyonemer = gsub("\\[", "", Seventyonemer)) %>%
  mutate(Seventyonemer = gsub("\\/.+]", "", Seventyonemer))

names <- data.frame(paste(">", affy_meta$Snpid, sep = ""))
fasta <- cbind(names, affy_meta$Seventyonemer)
fasta <- as.vector(t(fasta))
fasta <- as.data.frame(fasta)

write.table(fasta, "data/processed/flanks_submitted_87608.fasta", quote = F, col.names = F, row.names = F)

#~~ Eddie

# scp data/processed/flanks_submitted_87608.fasta ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/

# Running alignment on EDDIE (bwa and blast)

#~~~~~~~~~

#~~ BWA output

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/flanks_submitted_87608.sam data/processed/flanks_submitted_87608.sam 

# read in sam file
sam_file <- "data/processed/flanks_submitted_87608.sam"

# check where alignment info starts
num_to_skip <- system(paste0("grep -n @PG ", sam_file), intern = TRUE) %>% 
  str_split(":") %>% 
  .[[1]] %>% 
  .[1] %>% 
  as.numeric()

# read data from sam file
probe_sam <- read_delim(sam_file, skip = num_to_skip, 
                            delim = "\t", 
                            col_names = c("snp", "flag", "chr", "pos", "mapq",
                                          "cigar", "rnext", "pnext", "tlen",
                                          "seq", "qual", "edit_dist", "mismatch",
                                          "alignment_score", "subopt_alignment_score",
                                          "alt_hits")) %>% 
  filter(flag %in% c(0,16))

# filter snps mapping to dog X

probe_sam_X <- probe_sam %>%
  filter(chr == "X")

# write list of X-linked SNPs for discarding with PLINK

write.table(probe_sam_X$snp, "data/processed/x_linked_variants.txt",
            quote = F, row.names = F, col.names = F)



#~~ Compare to BLAST output

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/agaz_chip/probes2canislupus data/processed/probes2canislupus 

blast <- fread("data/processed/probes2canislupus")

blast_X <- blast %>%
  filter(V2 == "X")
