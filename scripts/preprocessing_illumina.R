# Load packages
library(illuminaHumanv4.db)
library(lumi)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(cowplot)
library(affy)

# Load studies description
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
illu <- which(grepl("Illu", studies$platform))
# Change this to appropriate index
i = illu[1]

# Read data
pd <- read.AnnotatedDataFrame(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""))
x.lumi <- lumiR(paste("../raws/", studies[i,]$ID, "/", studies[i,]$ID, "_non-normalized.tsv", sep=""),
                sep = "\t", QC=FALSE, columnNameGrepPattern = list(exprs='P'))

x.lumi <- x.lumi[, -which(grepl("replicate", colnames(x.lumi)))]
# Normalize
lumi.T <- lumiT(x.lumi, method = 'log2')
lumi.N<-lumiN(lumi.T, method = "quantile")
colnames(lumi.N) <- pd@data$SampleAccessionNumber
# Save results into file
write.table(exprs(lumi.N), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_illumina.tsv", sep=""),
            sep="\t", quote=FALSE)

# Add phenoData to ExpressionSet object
lumi.N@phenoData = pd
# Perform PCA and plot
pca = prcomp(t(exprs(lumi.N)))
title <- ggdraw() + draw_label(studies[i,]$ID, fontface='bold')  
pl <- autoplot(pca, data = pData(lumi.N), colour="CancerType")
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

# Save plot for manual quality control
save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA.pdf", sep=""),
          pl)

probesetsID <- rownames(exprs(lumi.N))
probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")

