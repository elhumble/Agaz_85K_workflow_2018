# APT outputs minimal PLINK files which need updating

library(dplyr)
library(data.table)
library(tidyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Update map files        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Join by Marker ID and Probe.Set.ID
# Fix _v1.1
# Fix JIH_SNP_ID

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

#~~ Rewrite all

map_all <- fread("data/out/all/plink/agaz_plate1_3_apt.map")

write.table(fix_map(map_all, annot), "data/out/all/plink/agaz_plate1_3.map", 
            quote = F, row.names = F, col.names = F)

#~~ Rewrite agaz

map_agaz <- fread("data/out/agaz/plink/agaz_plate1_3_apt.map")

write.table(fix_map(map_agaz, annot), "data/out/agaz/plink/agaz_plate1_3.map", 
            quote = F, row.names = F, col.names = F)


#~~~~~~~~~~~~~~~~~~~~~~~#
#      Update IDs       #
#~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Recode IDs in ped file

# Bash commands:

# Agaz

# awk 'gsub("\\_\\(Axiom\\_AGaz90k\\)\\_plate.+\.CEL", "", $1)1' agaz_plate1_3_apt.ped > test
# awk 'gsub("^.+\\-", "", $1)1' test > agaz_plate1_3.ped
# rm test
# Manually edited dups: _A and _B

# All: Not run

# awk 'gsub("\\_\\(Axiom\\_AGaz90k\\)\\_plate.+\.CEL", "", $1)1' data/out/all/plink/agaz_plate1_3_apt.ped > test
# awk 'gsub("^.+\\-", "", $1)1' test > data/out/all/plink/agaz_plate1_3.ped
# rm test


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Update pedigree info        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ PED info

# seals <- read.csv("data/SSBseals_wrangled_Dec_2017.csv") %>%
#   select(PUP, MOTHER, MOTHER_ID1) %>%
#   mutate(PUP = as.character(PUP),
#          MOTHER = as.character(MOTHER),
#          MOTHER_ID1 = as.character(MOTHER_ID1))
# 
# AT <- data.frame(PUP = "ATP16015",
#                  MOTHER = "ATF16015",
#                  MOTHER_ID1 = NA)
# 
# seals <- rbind(seals, AT) %>%
#   mutate(MOTHER = case_when(MOTHER_ID1 == "AGF99049" ~ "AGF99049",
#                             TRUE ~ MOTHER)) %>%
#   select(-MOTHER_ID1)
# 
# 
# ped_file <- read.table("data/out/agaz/plink/agaz_plate1_3.ped", skip  = 1)
# 
# ped <- as.data.frame(ped_file[,1])
# colnames(ped) <- "IID"
# 
# ped_update <- ped
# 
# ped <- ped %>%
#   left_join(seals, by = c("IID" = "PUP")) %>%
#   mutate(MOTHER = as.character(MOTHER)) %>%
#   mutate(a = MOTHER %in% IID) %>%
#   mutate(MOTHER = case_when(a == T ~ MOTHER,
#                             TRUE ~ "NA")) %>%
#   mutate(MOTHER = na_if(MOTHER, "NA")) %>%
#   mutate(FATHER = NA) %>%
#   select(-a)




#~~ Plate data from Anneke

seals <- read.csv("data/pup_dna_July2019.csv") %>%
  select(PUP, MOTHER)

# Add 2016 animals

AT <- data.frame(PUP = "ATP16015", 
                 MOTHER = "ATF16015")

seals <- rbind(seals, AT)

#~~ Read list of genotyped IDs from ped file

ped_file <- read.table("data/out/agaz/plink/agaz_plate1_3.ped")
ped <- as.data.frame(ped_file[,1])
colnames(ped) <- "IID"
IDs <- ped

#~~ NA mother ID if she has not been genotyped in plates 1-3

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
  dplyr::distinct(IID, .keep_all = T) %>%
  mutate(MOTHER = case_when(IID == "AGP01096" ~ "AGF00052", # lose this mother when running distinct()
                            TRUE ~ MOTHER)) %>%
  mutate(FATHER = replace_na(FATHER, 0)) %>%
  mutate(MOTHER = replace_na(MOTHER, 0))

# fam info

head(famped)

new_famped <- IDs %>%
  left_join(famped, by = "IID") # correct order

# correct MOTHER ID AGF99027 with AGF99049

new_famped <- new_famped %>%
  mutate(MOTHER = case_when(MOTHER == "AGF99027" ~ "AGF99049",
            TRUE ~ MOTHER))

nrow(filter(new_famped, MOTHER != "0"))


write.table(new_famped[c(4,1,2,3)], "data/out/agaz/plink/new_fam_for_plink.txt", col.names = F,
            row.names = F, quote = F)

# write out file with sex info 

sex <- new_famped[c(4,1)] %>%
  mutate(sex = 2)

write.table(sex, "data/out/agaz/plink/new_sex_for_plink.txt", col.names = F,
            row.names = F, quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Filter for polymorphic SNPs      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ Extract polymorphic SNPs

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3 ", 
              "--extract data/out/agaz/plink/poly_snps.txt ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly ",
              "--recode --allow-extra-chr --debug ",
              "--no-fid --no-pheno --no-sex --no-parents"))


#~~ update parents info using file from pedigree

system(paste0("~/software/plink --file data/out/agaz/plink/agaz_plate1_3_poly ", 
              "--update-parents data/out/agaz/plink/new_fam_for_plink.txt ",
              "--update-sex data/out/agaz/plink/new_sex_for_plink.txt ",
              "--out data/out/agaz/plink/agaz_plate1_3_poly_parents ",
              "--recode --allow-extra-chr --debug"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Write duplicate removal file      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# file for duplicate removal AGF14014 has highest genotyping rate

data.frame(FID = c("AGF14014_A", "AGF14014_B"),
                  IID = c("AGF14014_A", "AGF14014_B")) %>%
  write.table("data/out/agaz/plink/dups_to_rm_plink.txt", 
              quote = F, row.names = F, col.names = F)
