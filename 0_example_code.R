library("SNPolisher")


# THIS RUNS::::: <><><><><><><><><><><><><><><><><><><><><><>

Ps_Visualization(pidFile = "test_out/PolyHighResolution.ps",
                 output.File = "test_out/PolyHighResolution.pdf",
                 output.dir = "plots/", 
                 summaryFile = "test_out/AxiomGT1.summary.txt", 
                 callFile = "test_out/AxiomGT1.calls.txt", 
                 posteriorFile = "test_out/AxiomGT1.snp-posteriors.txt", 
                 temp.dir = "test_out/temp", 
                 plot.ref = FALSE, max.num.SNP.draw = 6, col.AA = "green", col.AB = "blue", col.BB = "grey")

#:::: <><><><><><><><><><><><><><><><><><><><><><>


### EXAMPLE 1: HUMAN

# Ps_Visualization
Ps_Visualization(pidFile = "test_out/PolyHighResolution.ps", 
                 output.File = "test_out/PolyHighResolution.pdf", 
                 output.dir = "test_out/", 
                 callFile = "test_out/AxiomGT1.calls.txt", 
                 summaryFile = "test_out/AxiomGT1.summary.txt", 
                 posteriorFile = "test_out/AxiomGT1.snp-posteriors.txt", 
                 #multiallele.posteriorFile = "data/example_data/example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 #multiallele.priorFile = "data/example_data/example_1_human/AxiomGT1.priors.multi.txt", 
                 #specialSNPsFile = "data/example_data/example_1_human/example1.specialSNPs", 
                 reportFile = "test_out/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, temp.dir = "test_out/Temp")

# no reference plots
Ps_Visualization(pidFile = "data/example_data/example_1_human/PolyHighResolution.ps", 
                 output.File = "PolyHighResolution_no_ref.pdf", 
                 output.dir = "data/example_data/out", 
                 callFile = "data/example_data/example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "data/example_data/example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "data/example_data/example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "data/example_data/example_1_human/AxiomGT1.snp-posteriors.multi.txt",
                 multiallele.priorFile = "data/example_data/example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "data/example_data/example_1_human/example1.specialSNPs", 
                 reportFile = "data/example_data/example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, temp.dir = "data/example_data/example_1_human/Temp", 
                 plot.ref = FALSE, keep.temp.dir = TRUE)

# changing color of BB cluster, highlighting samples, changing plot labels, using temporary directory
Ps_Visualization(pidFile = "example_1_human/PolyHighResolution.ps", 
                 output.File = "PolyHighResolution_colors.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, temp.dir = "example_1_human/Temp", 
                 plot.ref = FALSE, use.temp.dir = TRUE, col.BB = "darkgreen", 
                 sampleFile = "example_1_human/samples.txt",labelsFile = "example_1_human/labels.txt")

# plotting confidences
Ps_Visualization(pidFile = "example_1_human/PolyHighResolution.ps", 
                 output.File = "PolyHighResolution_confidences.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 confidenceFile = "example_1_human/AxiomGT1.confidences.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, temp.dir = "example_1_human/Temp",
                 plot.ref = FALSE, keep.temp.dir = TRUE, plot.calls = FALSE, plot.confs = TRUE)


# updated confidences cut-offs

confidences <- read.table("example_1_human/AxiomGT1.confidences.txt", 
                          header=T, stringsAsFactors = F)

summary(unlist(confidences[,-1]))
quantile(unlist(confidences[,-1]),c(0.75,0.8,0.85,0.9,0.95))

Ps_Visualization(pidFile = "example_1_human/PolyHighResolution.ps", 
                 output.File = "PolyHighResolution_confidences_thresholds.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 confidenceFile = "example_1_human/AxiomGT1.confidences.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE, 
                 keep.temp.dir = TRUE, plot.calls = FALSE, plot.confs = TRUE, 
                 confs.thresh=c(0,0.00001,0.00002,0.00003,1))

# updated confidences with calls
Ps_Visualization(pidFile = "example_1_human/PolyHighResolution.ps", 
                 output.File = "PolyHighResolution_confidences_calls.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 confidenceFile = "example_1_human/AxiomGT1.confidences.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 temp.dir = "example_1_human/Temp", 
                 plot.ref = FALSE, keep.temp.dir = TRUE, 
                 plot.intensity = FALSE, plot.calls = TRUE, plot.confs = TRUE, 
                 confs.thresh=c(0,0.00001,0.00002,0.00003,1))




# Ps_Vis_Density
Ps_Vis_Density(pidFile = "example_1_human/PolyHighResolution.ps", 
               output.File = "PolyHighResolution_density.pdf", 
               output.dir = "example_1_human/plots", 
               callFile = "example_1_human/AxiomGT1.calls.txt", 
               summaryFile = "example_1_human/AxiomGT1.summary.txt", 
               posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
               multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
               multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
               specialSNPsFile = "example_1_human/example1.specialSNPs", 
               reportFile = "example_1_human/AxiomGT1.report.txt", 
               plot.prior = TRUE, 
               plot.posterior = TRUE, 
               temp.dir = "example_1_human/Temp", 
               plot.ref = FALSE, keep.temp.dir = TRUE)

# no priors or posteriors
Ps_Vis_Density(pidFile = "example_1_human/PolyHighResolution.ps", 
               output.File = "PolyHighResolution_density_no_posteriors.pdf", 
               output.dir = "example_1_human/plots", 
               callFile = "example_1_human/AxiomGT1.calls.txt", 
               summaryFile = "example_1_human/AxiomGT1.summary.txt", 
               posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
               multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
               multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
               specialSNPsFile = "example_1_human/example1.specialSNPs", 
               reportFile = "example_1_human/AxiomGT1.report.txt", 
               plot.prior = FALSE, plot.posterior = FALSE, temp.dir = "example_1_human/Temp", 
               plot.ref = FALSE, keep.temp.dir = TRUE, vertical.line=TRUE)

# densities in green with 5 levels
Ps_Vis_Density(pidFile = "example_1_human/PolyHighResolution.ps", 
               output.File = "PolyHighResolution_density_green.pdf", 
               output.dir = "example_1_human/plots", 
               callFile = "example_1_human/AxiomGT1.calls.txt", 
               summaryFile = "example_1_human/AxiomGT1.summary.txt", 
               posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
               multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
               multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
               specialSNPsFile = "example_1_human/example1.specialSNPs", 
               reportFile = "example_1_human/AxiomGT1.report.txt", 
               plot.prior = FALSE, plot.posterior = FALSE, 
               temp.dir = "example_1_human/Temp", 
               plot.ref = FALSE, keep.temp.dir = TRUE, vertical.line=TRUE, 
               density.col = "green", density.levels = 5)





# plotting all of the categories from Ps_Classification

Ps_Visualization(pidFile = "example_1_human/NoMinorHom.ps", 
                 output.File = "NoMinorHom.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)

Ps_Visualization(pidFile = "example_1_human/MonoHighResolution.ps", 
                 output.File = "MonoHighResolution.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)

Ps_Visualization(pidFile = "example_1_human/OffTargetVariant.ps", 
                 output.File = "OffTargetVariant.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)

Ps_Visualization(pidFile = "example_1_human/CallRateBelowThreshold.ps", 
                 output.File = "CallRateBelowThreshold.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)

Ps_Visualization(pidFile = "example_1_human/Other.ps", 
                 output.File = "Other.pdf", 
                 output.dir = "example_1_human/plots", 
                 callFile = "example_1_human/AxiomGT1.calls.txt", 
                 summaryFile = "example_1_human/AxiomGT1.summary.txt", 
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt", 
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt", 
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt", 
                 specialSNPsFile = "example_1_human/example1.specialSNPs", 
                 reportFile = "example_1_human/AxiomGT1.report.txt", 
                 plot.prior = TRUE, plot.posterior = TRUE, 
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)




# OTV_Caller
# making original plots of OTV SNPs
Ps_Visualization(pidFile = "example_1_human/OffTargetVariant.ps",
                 output.File = "OffTargetVariant.pdf",
                 output.dir = "example_1_human/plots",
                 callFile = "example_1_human/AxiomGT1.calls.txt",
                 summaryFile = "example_1_human/AxiomGT1.summary.txt",
                 posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.txt",
                 multiallele.posteriorFile = "example_1_human/AxiomGT1.snp-posteriors.multi.txt",
                 multiallele.priorFile = "example_1_human/AxiomGT1.priors.multi.txt",
                 specialSNPsFile = "example_1_human/example1.specialSNPs",
                 reportFile = "example_1_human/AxiomGT1.report.txt",
                 plot.prior = TRUE, plot.posterior = TRUE,
                 temp.dir = "example_1_human/Temp", plot.ref = FALSE)

# running OTV_Caller
OTV_Caller(summaryFile="example_1_human/AxiomGT1.summary.txt", 
           posteriorFile="example_1_human/AxiomGT1.snp-posteriors.txt", 
           callFile="example_1_human/AxiomGT1.calls.txt", 
           confidenceFile="example_1_human/AxiomGT1.confidences.txt", 
           pidFile="example_1_human/OffTargetVariant.ps", 
           output.dir="example_1_human/OTV", 
           OTV.only=TRUE)

# Ps_Visualization on OTV results
Ps_Visualization(pidFile="example_1_human/OTV/OTV.keep.ps", 
                 output.File="example_1_human/OTV/OTV_results.pdf", 
                 summaryFile="example_1_human/AxiomGT1.summary.txt", 
                 callFile="example_1_human/OTV/OTV.calls.txt", 
                 posteriorFile="example_1_human/OTV/OTV.snp-posteriors.txt", 
                 temp.dir="example_1_human/OTV/Temp", 
                 keep.temp.dir=TRUE, plot.ref=FALSE)

# Ps_Visualization with original data on OTV results
Ps_Visualization(pidFile="example_1_human/OTV/OTV.keep.ps", 
                 output.File="example_1_human/OTV/Cluster_OTV_original.pdf", 
                 summaryFile="example_1_human/AxiomGT1.summary.txt", 
                 callFile="example_1_human/AxiomGT1.calls.txt", 
                 posteriorFile="example_1_human/AxiomGT1.snp-posteriors.txt", 
                 temp.dir="example_1_human/OTV/Temp", 
                 keep.temp.dir=TRUE, plot.ref=FALSE)





### EXAMPLE 2: WHEAT


library("SNPolisher")

axiom_dir <- "C:/Users/user_name/Desktop/example_2_wheat/"
output_dir <- "C:/Users/user_name/Desktop/example_2_wheat/Output/"
otv_dir <- "C:/Users/user_name/Desktop/example_2_wheat/OTV/"
adjust_dir <- "C:/Users/user_name/Desktop/example_2_wheat/Call_Adjust/"


# using the paste function
paste("Lorem ipsum dolor sit amet,", "consectetur adipiscing elit.")
paste(axiom_dir,"AxiomGT1.calls.txt")
paste(axiom_dir,"AxiomGT1.calls.txt", sep="")


# using the paste function with a slash
new_dir <- "C:/Users/user_name/Desktop/example_2_wheat"
paste(new_dir,"AxiomGT1.calls.txt", sep="")
paste(new_dir,"AxiomGT1.calls.txt", sep="/")



# Ps_Visualization
Ps_Visualization("PolyHighResolution.ps",
                 output.File = paste(output_dir, "PolyHighResolution.pdf",sep=""),
                 summaryFile = paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), 
                 plot.ref=FALSE, max.num.SNP.draw=6)

# THIS RUNS::::: <><><><><><><><><><><><><><><><><><><><><><>

Ps_Visualization(pidFile = "test_out/PolyHighResolution.ps",
                 output.File = "test_out/PolyHighResolution.pdf",
                 output.dir = "plots/", 
                 summaryFile = "test_out/AxiomGT1.summary.txt", 
                 callFile = "test_out/AxiomGT1.calls.txt", 
                 posteriorFile = "test_out/AxiomGT1.snp-posteriors.txt", 
                 temp.dir = "test_out/temp", 
                 plot.ref = FALSE, max.num.SNP.draw = 6, col.AA = "green", col.AB = "blue", col.BB = "grey")



# changing the color of the heterozygous cluster
Ps_Visualization("PolyHighResolution.ps",
                 output.File=paste(output_dir, "PolyHighResolution_color.pdf",sep=""),
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""),
                 temp.dir=paste(output_dir,"Temp",sep=""), 
                 plot.ref=FALSE, max.num.SNP.draw=6, col.AB="green")

# plotting many SNPs at once
Ps_Visualization("PolyHighResolution.ps",
                 output.File=paste(output_dir, "PolyHighResolution_300.pdf",sep=""),
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), plot.ref=FALSE)




# OTV_Caller
OTV_Caller(summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
           posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
           callFile=paste(axiom_dir, "AxiomGT1.calls.txt",sep=""), 
           confidenceFile=paste(axiom_dir, "AxiomGT1.confidences.txt",sep=""), 
           pidFile=paste(output_dir,"OffTargetVariant.ps",sep=""), 
           output.dir=otv_dir, OTV.only=TRUE, keep.orig=FALSE)

# Ps_Visualization on OTV results
Ps_Visualization(pidFile=paste(otv_dir,"OTV.keep.ps",sep=""), 
                 output.File=paste(otv_dir,"Cluster_OTV.pdf",sep=""),
                 summaryFile=paste(otv_dir,"OTV.summary.txt",sep=""), 
                 callFile=paste(otv_dir,"OTV.calls.txt",sep=""), 
                 posteriorFile=paste(otv_dir,"OTV.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(otv_dir,"Temp",sep=""), plot.ref=FALSE)




# Ps_CallAdjust preparation
Ps_Visualization(pidFile="CallRateBelowThreshold.ps", 
                 output.File=paste(output_dir,"CRBT.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), 
                 keep.temp.dir=TRUE, plot.ref=FALSE, plot.intensity=FALSE)


confs <- read.table(paste(axiom_dir,"AxiomGT1.confidences.txt",sep=""),header=T)
probesets <- read.table(paste(axiom_dir,"call_adjust_probesets.ps",sep=""),header=T)
sum(confs[,1]%in%probesets[,1])
confs <- confs[confs[,1]%in%probesets[,1],]
quantile(unlist(confs[,2:dim(confs)[2]]))

Ps_Visualization(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
                 output.File=paste(output_dir,"CRBT_confs.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 confidenceFile=paste(axiom_dir,"AxiomGT1.confidences.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), 
                 keep.temp.dir=TRUE, plot.ref=FALSE, plot.intensity=FALSE, plot.calls=TRUE, 
                 plot.confs=TRUE, confs.thresh=c(0,0.00004,0.00159,0.01799,1))

Ps_Visualization(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
                 output.File=paste(output_dir,"CRBT_2.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 confidenceFile=paste(axiom_dir,"AxiomGT1.confidences.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), 
                 keep.temp.dir=TRUE, plot.ref=FALSE, plot.intensity=FALSE, 
                 confidences.cutoff=0.02)

Ps_Visualization(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
                 output.File=paste(output_dir,"CRBT_3.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 confidenceFile=paste(axiom_dir,"AxiomGT1.confidences.txt",sep=""),
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(output_dir,"Temp",sep=""), keep.temp.dir=TRUE, 
                 plot.ref=FALSE, plot.intensity=FALSE, confidences.cutoff=0.01)


# Ps_CallAdjust
Ps_CallAdjust(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
              callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
              confidenceFile=paste(axiom_dir,"AxiomGT1.confidences.txt",sep=""), 
              threshold=0.015, 
              outputFile=paste(adjust_dir,"call_adjust_calls.txt",sep=""))


# redo plots with new calls
Ps_Visualization(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
                 output.File=paste(adjust_dir,"call_adjust_new.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(adjust_dir,"call_adjust_calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(adjust_dir,"Temp",sep=""), 
                 keep.temp.dir=TRUE, plot.ref=FALSE, plot.intensity=FALSE)

# making plot of SNPs with old calls for comparison
Ps_Visualization(pidFile=paste(axiom_dir,"call_adjust_probesets.ps",sep=""), 
                 output.File=paste(adjust_dir,"call_adjust_old.pdf",sep=""), 
                 summaryFile=paste(axiom_dir,"AxiomGT1.summary.txt",sep=""), 
                 callFile=paste(axiom_dir,"AxiomGT1.calls.txt",sep=""), 
                 posteriorFile=paste(axiom_dir,"AxiomGT1.snp-posteriors.txt",sep=""), 
                 temp.dir=paste(adjust_dir,"Temp",sep=""), 
                 keep.temp.dir=TRUE, plot.ref=FALSE, plot.intensity=FALSE)






### EXAMPLE 3: B ALLELE FREQUENCY TEST

library("SNPolisher")

# running the B allele frequency test
BalleleFreq_Test(callFile.small="example_3_ballele/case_calls.txt", 
                 callFile.full="example_3_ballele/casecontrol_calls.txt", 
                 pidFile="example_3_ballele/pid.ps", 
                 sampleList="example_3_ballele/samples.txt", 
                 intercept = 0.02, 
                 output.file="example_3_ballele/BalleleFreq_test.txt")

# inspecting the results
tests <- read.table("example_3_ballele/BalleleFreq_test.txt",header=T)
tests[1:5,]
tests[tests[,4]=="fail",1]
tests$probeset_id[tests$passBalleleFreqTest=="fail"]

sum(tests[,4]=="fail")
sum(tests[,4]=="pass")

# saving the markers that failed the test
write.table(tests[tests[,4]=="fail",1], file="example_3_ballele/markers_fail.txt", 
            sep=" ",row.names=FALSE,col.names="probeset_id",quote=FALSE)

# saving the markers that passed the test
write.table(tests[tests[,4]=="pass",1],file="example_3_ballele/markers_pass.txt",
            sep=" ",row.names=FALSE,col.names="probeset_id",quote=FALSE)


# generating the plots for the markers that failed the test
Ps_Visualization(pidFile="example_3_ballele/markers_fail.txt", 
                 output.File="example_3_ballele/failed.pdf", 
                 sampleFile="example_3_ballele/samples.txt", 
                 summaryFile="example_3_ballele/casecontrol_summary.txt", 
                 callFile="example_3_ballele/casecontrol_calls.txt", 
                 posteriorFile="example_3_ballele/casecontrol_posts.txt", 
                 temp.dir="example_3_ballele/Temp", plot.ref=FALSE)

# generating the plots for the markers that passed the test
Ps_Visualization(pidFile="example_3_ballele/markers_pass.txt", 
                 output.File="example_3_ballele/passed.pdf", 
                 sampleFile="example_3_ballele/samples.txt", 
                 summaryFile="example_3_ballele/casecontrol_summary.txt", 
                 callFile="example_3_ballele/casecontrol_calls.txt", 
                 posteriorFile="example_3_ballele/casecontrol_posts.txt", 
                 temp.dir="example_3_ballele/Temp", plot.ref=FALSE)





### EXAMPLE 4: AUTOTETRAPLOID SPECIES

library("SNPolisher")

# reformat Axiom summary file for fitTetra use
fitTetra_Input(summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt",
               output.file="example_4_autotetraploid/AxiomGT1.summary.fitTetra.txt")

# specify the number of markers per file and use a pidFile
fitTetra_Input(summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt",
               output.file="example_4_autotetraploid/AxiomGT1.summary.fitTetra.txt", 
               output.count=10)

# use a pidFile
fitTetra_Input(summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt",
               output.file="example_4_autotetraploid/AxiomGT1.summary.fitTetra.small.txt", 
               pidFile="example_4_autotetraploid/keep_fitTetra.ps")

# using the gzipped file
fitTetra_Input(summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt.gz",
               output.file="example_4_autotetraploid/AxiomGT1.summary.fitTetra.small.txt", 
               pidFile="example_4_autotetraploid/keep_fitTetra.ps")


# if the fitTetra package is not installed, do it now.
install.packages("fitTetra")
library("fitTetra")


# run fitTetra function "saveMarkerModels" to produce score file
# example uses 4 cores
markerData<-read.delim("example_4_autotetraploid/AxiomGT1.summary.fitTetra.txt", 
                       sep="\t",stringsAsFactors=F)

saveMarkerModels(data=markerData, 
                 maxiter=500, ncores=1, 
                 try.HW=F, dip.filter=F, 
                 p.threshold=0.75, call.threshold=0.01, 
                 logfile="example_4_autotetraploid/polyploid.log.txt", 
                 modelfile="example_4_autotetraploid/polyploid.model.txt", 
                 allmodelsfile="example_4_autotetraploid/polyploid.allmodel.txt", 
                 scorefile="example_4_autotetraploid/polyploid.score.txt", 
                 diploscorefile="", plot="fitted", plot.type="png") 

# find the levels of the markers
names(markerData)
levels(as.factor(markerData$MarkerName))
levels(as.factor(markerData$MarkerName))[c(3,8,17,24)]

# reformat fitTetra files for SNPolisher use
fitTetra_Output(scoreFile="example_4_autotetraploid/polyploid.score.txt", 
                summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt", 
                output.callFile="AxiomGT1.calls.fitTetra.txt", 
                output.confFile="AxiomGT1.confidences.fitTetra.txt",
                output.postFile="AxiomGT1.snp-posteriors.fitTetra.txt",
                output.logFile="log.txt", 
                conf.threshold=1, 
                output.dir="example_4_autotetraploid/fitTetra_Output")

# using the gzipped file
fitTetra_Output(scoreFile="example_4_autotetraploid/polyploid.score.txt", 
                summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt.gz", 
                output.callFile="AxiomGT1.calls.fitTetra.txt",
                output.confFile="AxiomGT1.confidences.fitTetra.txt", 
                output.postFile="AxiomGT1.snp-posteriors.fitTetra.txt", 
                output.logFile="log.txt", conf.threshold=1, 
                output.dir="example_4_autotetraploid/fitTetra_Output")

# reformat fitTetra files and remove several markers
fitTetra_Output(scoreFile="example_4_autotetraploid/polyploid.score.txt", 
                summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt", 
                output.callFile="AxiomGT1.calls.fitTetra.txt", 
                output.confFile="AxiomGT1.confidences.fitTetra.txt", 
                output.postFile="AxiomGT1.snp-posteriors.fitTetra.txt", 
                output.logFile="log.txt", conf.threshold=1, 
                output.dir="example_4_autotetraploid/fitTetra_Output_2", 
                pidFile="example_4_autotetraploid/keep_fitTetra.ps")

# plot the new calls
Ps_Visualization(pidFile="example_4_autotetraploid/all_fitTetra.ps", 
                 output.File="example_4_autotetraploid/fitTetra_Output/fittetra_output.pdf", 
                 summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt", 
                 callFile="example_4_autotetraploid/fitTetra_Output/AxiomGT1.calls.fitTetra.txt", 
                 posteriorFile="example_4_autotetraploid/fitTetra_Output/AxiomGT1.snp-posteriors.fitTetra.txt", plot.ref=FALSE)

# cut down the summary file to the probesets that passed fitTetra
summary <- read.table("example_4_autotetraploid/AxiomGT1.summary.txt", header=T, stringsAsFactors=FALSE)
calls <- read.table("example_4_autotetraploid/fitTetra_Output/AxiomGT1.calls.fitTetra.txt", header=T, stringsAsFactors=FALSE)
check <- strsplit(summary[,1],"-A")
keep <- c(which(check%in%calls[,1]),which(check%in%calls[,1])+1)
keep <- keep[order(keep)]

write.table(summary[keep,],file="example_4_autotetraploid/fitTetra_Output/AxiomGT1.summary.small.txt",
            sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

write.table(calls[,1],file="example_4_autotetraploid/fitTetra_Output/fittetra.ps",
            sep="\t",row.names=FALSE,col.names="probeset_id",quote=FALSE)

# plot the new calls
Ps_Visualization(pidFile="example_4_autotetraploid/fitTetra_Output/fittetra.ps", 
                 output.File="example_4_autotetraploid/fitTetra_Output/fittetra_output.pdf",
                 summaryFile="example_4_autotetraploid/fitTetra_Output/AxiomGT1.summary.small.txt", 
                 callFile="example_4_autotetraploid/fitTetra_Output/AxiomGT1.calls.fitTetra.txt",
                 posteriorFile="example_4_autotetraploid/fitTetra_Output/AxiomGT1.snp-posteriors.fitTetra.txt", plot.ref=FALSE)

# rerun with higher p.threshold value
saveMarkerModels(data=markerData, maxiter=500, 
                 ncores=1, try.HW=F, dip.filter=F, 
                 p.threshold=0.99, call.threshold=0.01, 
                 logfile="example_4_autotetraploid/polyploid.log.2.txt", 
                 modelfile="example_4_autotetraploid/polyploid.model.2.txt", 
                 allmodelsfile="example_4_autotetraploid/polyploid.allmodel.2.txt", 
                 scorefile="example_4_autotetraploid/polyploid.score.2.txt", 
                 diploscorefile="", plot="fitted", plot.type="png") 

# reformat fitTetra files for SNPolisher use
fitTetra_Output(scoreFile="example_4_autotetraploid/polyploid.score.2.txt", 
                summaryFile="example_4_autotetraploid/AxiomGT1.summary.txt", 
                output.callFile="AxiomGT1.calls.fitTetra.txt", 
                output.confFile="AxiomGT1.confidences.fitTetra.txt", 
                output.postFile="AxiomGT1.snp-posteriors.fitTetra.txt", 
                output.logFile="log.txt", conf.threshold=1, 
                output.dir="example_4_autotetraploid/fitTetra_Output_3")

# plot the new calls
Ps_Visualization(pidFile="example_4_autotetraploid/fitTetra_Output/fittetra.ps", 
                 output.File="example_4_autotetraploid/fitTetra_Output_3/fittetra_output.pdf", 
                 summaryFile="example_4_autotetraploid/fitTetra_Output_3/AxiomGT1.summary.fitTetra.txt", 
                 callFile="example_4_autotetraploid/fitTetra_Output_3/AxiomGT1.calls.fitTetra.txt", 
                 posteriorFile="example_4_autotetraploid/fitTetra_Output_3/AxiomGT1.snp-posteriors.fitTetra.txt", plot.ref=FALSE)






### EXAMPLE 5: BATCH PLOTTING

library("SNPolisher")

# create batch directory names
call_files <- c("example_5_batch_plotting/batch1/calls.txt",
                "example_5_batch_plotting/batch2/calls.txt",
                "example_5_batch_plotting/batch3/calls.txt",
                "example_5_batch_plotting/batch4/calls.txt",
                "example_5_batch_plotting/batch5/calls.txt",
                "example_5_batch_plotting/batch6/calls.txt",
                "example_5_batch_plotting/batch7/calls.txt",
                "example_5_batch_plotting/batch8/calls.txt")

axiom_dirs <- c("example_5_batch_plotting/batch1/",
                "example_5_batch_plotting/batch2/",
                "example_5_batch_plotting/batch3/",
                "example_5_batch_plotting/batch4/",
                "example_5_batch_plotting/batch5/",
                "example_5_batch_plotting/batch6/",
                "example_5_batch_plotting/batch7/",
                "example_5_batch_plotting/batch8/")

call_files <- paste(axiom_dirs,"calls.txt",sep="")
confidence_files <- paste(axiom_dirs,"confidences.txt",sep="")
posterior_files <- paste(axiom_dirs,"posteriors.txt",sep="")
summary_files <- paste(axiom_dirs,"summary.txt",sep="")

# plot batch plots
Ps_Visualization(pidFile="example_5_batch_plotting/probesets.ps", 
                 output.File="example_5_batch_plotting/batches.pdf", 
                 callFile=call_files, summaryFile=summary_files, 
                 posteriorFile=posterior_files)

temp_dirs <- paste(axiom_dirs,"Temp",sep="")

Ps_Visualization(pidFile="example_5_batch_plotting/probesets.ps", 
                 output.File="example_5_batch_plotting/batches.pdf", 
                 callFile=call_files, summaryFile=summary_files, 
                 posteriorFile=posterior_files, temp.dir=temp_dirs)

# update the number of columns
Ps_Visualization(pidFile="example_5_batch_plotting/probesets.ps", 
                 output.File="example_5_batch_plotting/batches_8_cols.pdf", 
                 callFile=call_files, summaryFile=summary_files, 
                 posteriorFile=posterior_files, temp.dir=temp_dirs, num.cols=8)

# highlight samples and update plot labels for 3 batches
axiom_dirs <- c("example_5_batch_plotting/batch1/",
                "example_5_batch_plotting/batch2/",
                "example_5_batch_plotting/batch3/")

call_files <- paste(axiom_dirs,"calls.txt",sep="")
confidence_files <- paste(axiom_dirs,"confidences.txt",sep="")
posterior_files <- paste(axiom_dirs,"posteriors.txt",sep="")
summary_files <- paste(axiom_dirs,"summary.txt",sep="")
temp_dirs <- paste(axiom_dirs,"Temp",sep="")

labels_files <- paste(axiom_dirs,"labels_batch_",1:3,".txt",sep="")
samples_files <- c("example_5_batch_plotting/batch1/samples_batch_1.txt",NA,NA)
Ps_Visualization(pidFile="example_5_batch_plotting/probesets.ps", 
                 output.File="example_5_batch_plotting/batches_3.pdf", 
                 callFile=call_files, summaryFile=summary_files, 
                 posteriorFile=posterior_files, temp.dir=temp_dirs, 
                 sampleFile=samples_files, labelsFile=labels_files, 
                 plot.ref=FALSE, plot.intensity=FALSE)

# using Ps_Extract
for(i in 1:length(axiom_dirs)){
  Ps_Extract(pidFile="example_5_batch_plotting/multiallelic_SNPs.ps", 
             output.dir=axiom_dirs[i], 
             summaryFile=paste(axiom_dirs[i],"summary.txt",sep=""), 
             output.summary="summary_small.txt", 
             callFile=paste(axiom_dirs[i],"calls.txt",sep=""), 
             output.calls="calls_small.txt")
}

call_files <- paste(axiom_dirs,"calls_small.txt",sep="")
summary_files <- paste(axiom_dirs,"summary_small.txt",sep="")
posteriors_files <- rep(NA,3)

Ps_Visualization(pidFile="example_5_batch_plotting/multiallelic_SNPs.ps", 
                 output.File="example_5_batch_plotting/batches_multiallelic.pdf", 
                 callFile=call_files, summaryFile=summary_files, 
                 posteriorFile=posterior_files, temp.dir=temp_dirs, 
                 sampleFile=samples_files, labelsFile=labels_files, 
                 plot.ref=FALSE, plot.intensity=FALSE)




### EXAMPLE 6: 3D PLOTTING

library("SNPolisher")
install.packages("plotly")
library("plotly")

# read in data and organize it for use with plot.ly
calls <- read.table("example_6_3D_plotting/calls.txt", header=T, stringsAsFactors=F)
summary <- read.table("example_6_3D_plotting/summary.txt", header=T, stringsAsFactors=F)

# create the data object for SNP71
data <- as.data.frame(cbind(t(calls[which(calls$probeset_id=="SNP71"),-1]),
                            t(summary[which(summary$probeset_id=="SNP71-A"),-1]),
                            t(summary[which(summary$probeset_id=="SNP71-B"),-1]),
                            t(summary[which(summary$probeset_id=="SNP71-C"),-1])))

names(data) <- c("calls","summary_A","summary_B","summary_C")

# set up the genotypes to be a factor
data$genotypes <- rep(NA,length(data$calls))
data$genotypes[data$calls==2] <- "BB"
data$genotypes[data$calls==8] <- "BC"
data$genotypes[data$calls==9] <- "CC"
data$genotypes <- as.factor(data$genotypes)

# transform the data in log2 space
data$A <- log2(data$summary_A)
data$B <- log2(data$summary_B)
data$C <- log2(data$summary_C)

# set up the default colors
multi.col <- c("#00FF00","#FF00FF","#000080","#006400",
               "#BA55D3","#B8860B","#FF1493","#556B2F",
               "#FF8C00","#DC143C","#000000","#8B4513")

multi.pch <- c(24,21,25,23,22,24,21,25,23,22,24,21,25,23,22)

# match the colors up to the genotypes
# this is only necessary if you want to control which factor each color is assigned to
# otherwise just use colors=multi.col[1:3] in the command

cols <- multi.col[1:3]
cols <- setNames(cols, c("BB", "BC","CC"))

# plot the 3D scatterplot
plot_ly(data, x=~A, y=~B, z=~C) %>%
  add_markers(type = 'scatter', mode = 'markers', marker = list(size = 9), 
              symbol = ~genotypes, symbols = I(multi.pch[1:3]), 
              color=~genotypes, colors=cols)

