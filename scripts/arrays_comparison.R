# Loading libraries
library(lmtest)
library(stats)
library(cowplot)
library(plyr)
library(ggplot2)
# Loading custom scripts
source("scatterPlot.R")
source("diagnosticPlots.R")
source("barPlot.R")

### Initial info section

# Load studies description file
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
studies <- studies[-4,]

# Hgu versus hgu, no pipeline
pdata1 <-read.table(paste("../pdata/pdata_", studies[1,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
pdata2 <-read.table(paste("../pdata/pdata_", studies[2,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")

eset1<-read.table(paste("../preprocessed/", studies[1,]$ID, "_preprocessed_affymetrix.tsv", sep=""), header=TRUE)
eset2<-read.table(paste("../preprocessed/", studies[2,]$ID, "_preprocessed_affymetrix.tsv", sep=""), header=TRUE)

pdata1.tnbc<-pdata1[which(pdata1$CancerType=="TNBC"),]
pdata1.tnbc<-pdata1.tnbc[which(is.na(pdata1.tnbc$Outliers)),]
eset1.tnbc <- eset1[,colnames(eset1) %in% pdata1.tnbc$SampleAccessionNumber]

pdata2.tnbc<-pdata2[which(pdata2$CancerType=="TNBC"),]
pdata2.tnbc<-pdata2.tnbc[which(is.na(pdata2.tnbc$Outliers)),]
eset2.tnbc <- eset2[,colnames(eset2) %in% pdata2.tnbc$SampleAccessionNumber]

eset1.tnbc <- data.frame(val=rowMeans(eset1.tnbc), row.names=rownames(eset1.tnbc))
eset2.tnbc <- data.frame(val=rowMeans(eset2.tnbc), row.names=rownames(eset2.tnbc))

df1 <- data.frame(eset=eset1.tnbc$val, rnaseq=eset2.tnbc$val)

# This has nothing to do with rna-seq data, it' just the name of the variable
rnaseq <- df1$rnaseq
eset <- df1$eset
lin<-lm(rnaseq~eset)
cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  

pl1 <- scatterPlot(df1, paste(studies[1,]$ID, " versus ", studies[2,]$ID, sep=""),
                        lin, cub, c(2, 15, 2, 15))

# hgu versus hugene, brainarray pipeline
eset1<-read.table(paste("../allsamples_exprs/", studies[1,]$ID, "_allsamples_brainarray.tsv", sep=""), header=TRUE)
eset1.tnbc <- eset1[,colnames(eset1) %in% pdata1.tnbc$SampleAccessionNumber]

pdata3 <-read.table(paste("../pdata/pdata_", studies[3,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
pdata3.tnbc<-pdata3[which(pdata3$CancerType=="TNBC"),]
pdata3.tnbc<-pdata3.tnbc[which(is.na(pdata3.tnbc$Outliers)),]

eset3<-read.table(paste("../allsamples_exprs/", studies[3,]$ID, "_allsamples_brainarray.tsv", sep=""), header=TRUE)
eset3.tnbc <- eset3[,colnames(eset3) %in% pdata3.tnbc$SampleAccessionNumber]

eset1.tnbc <- data.frame(val=rowMeans(eset1.tnbc), row.names=rownames(eset1.tnbc))
eset3.tnbc <- data.frame(val=rowMeans(eset3.tnbc), row.names=rownames(eset3.tnbc))

df1 <- data.frame(eset=eset1.tnbc$val, rnaseq=eset3.tnbc$val)

rnaseq <- df1$rnaseq
eset <- df1$eset
lin<-lm(rnaseq~eset)
cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  

pl2 <- scatterPlot(df1, paste(studies[1,]$ID, " versus ", studies[3,]$ID, ". Brainarray", sep=""),
                   lin, cub, c(2, 15, 2, 15))


pl <- plot_grid(pl1, pl2, ncol=2, align="v")
# Save plot for manual quality control
save_plot(paste("../plots/arrays_comparison.pdf", sep=""),
          pl, base_height=6, base_aspect_ratio=1.9, nrow=1)

