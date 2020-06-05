library(sequoia)
library(data.table)

Geno <- GenoConvert(InFile = "data/out/agaz/plink/sequoia.raw")
LifeHist <- read.table("data/processed/lifehist.txt", header = F) # won't run with data.table

# Added dummy year to ATF female. When NA mother doesn't get assigned

# duplicate check and parentage assignment

ParOUT <- sequoia(GenoM = Geno, LifeHistData = LifeHist,
                  MaxSibIter = 0, Err = 0.004,
                  quiet = F, Plot = T)


stats <- SnpStats(Geno, ParOUT$PedigreePar)
MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1-stats[,"AF"])

# Polish dataset

# Remove one duplicate from each dup pair
Geno2 <- Geno[!rownames(Geno) %in% ParOUT$DupGenotype$ID2, ]

# Drop SNPs with high error rate and / or low MAF
Geno2 <- Geno2[, -which(stats[,"Err.hat"]>0.05 | MAF < 0.1)]

# Drop low call rate samples
Indiv.Mis <- apply(Geno2, 1, function(x) sum(x == -9)) / ncol(Geno2)

Geno2 <- Geno2[Indiv.Mis < 0.2, ]

# Run full pedigree reconstruction

SeqOUT <- sequoia(GenoM = Geno2,
                  LifeHistData = LifeHist,
                  MaxSibIter = 20,
                  Err = 0.004)

# Inspect assigned parents

SummarySeq(SeqOUT)

# Get maybe rel

Maybe <- GetMaybeRel(GenoM = Geno2,
                     Pedigree = SeqOUT$PedigreePar, ParSib="sib")


# Save file

save(SeqOUT, file = "data/processed/SeqOUT.Rdata")
save(Maybe, file = "data/processed/MaybeRel.Rdata")


