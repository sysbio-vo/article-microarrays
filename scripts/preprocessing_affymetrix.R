# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
# Other packages
library(affy)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(cowplot)
# Load studies description
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
affy <- which(grepl("Affy", studies$platform))
# Change this to appropriate index
i = affy[3]

## Affymetrix data preprocessing
# Load pdata and CEL files
pd <- read.AnnotatedDataFrame(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""))
affyData = ReadAffy(phenoData=pd, sampleNames=pd$SampleAccessionNumber, filenames=as.character(rownames(pd)),
                    celfile.path=paste("../raws/", studies[i,]$ID, sep=""))
# Affymetrix probesets. RMA
#affyData@cdfName <- "HG-U133_Plus_2"
#affyData@cdfName <- "HuGene-1_0-st-v1"
eset = rma(affyData)
# Brainarray probesets. RMA
#affyData@cdfName <- "hgu133plus2hsentrezgcdf"
affyData@cdfName <- "hugene10sthsentrezgcdf"
eset.br = rma(affyData)
# Save affy and brain expression sets
write.table(exprs(eset), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
write.table(exprs(eset.br), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), sep="\t", quote=FALSE)

## QC
# Add ScanDate to pdata
pData(eset)$ScanDate <- str_replace_all(eset@protocolData@data$ScanDate, "T", " ")
pData(eset)$ScanDate <- sapply(strsplit(pData(eset)$ScanDate, split=' ', fixed=TRUE), function(x) (x[1]))
pData(eset.br)$ScanDate <- pData(eset)$ScanDate
#pData(eset.br)$ScanDate <- substr(pData(eset.br)$ScanDate, 1, 7)
#pData(eset)$ScanDate <- substr(pData(eset)$ScanDate, 1, 7)

# Perform PCA
pca = prcomp(t(exprs(eset)))
pca.br = prcomp(t(exprs(eset.br)))

title <- ggdraw() + draw_label(paste("Affymetrix Probesets Definitions.", studies[i,]$ID), fontface='bold')  
pl1 <- autoplot(pca, data = pData(eset), colour="CancerType")
pl2 <- autoplot(pca, data = pData(eset), colour="ScanDate")
pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label(paste("Brainarray Probesets Definitions.", studies[i,]$ID), fontface='bold')  
pl1 <- autoplot(pca.br, data = pData(eset.br), colour="CancerType")
pl2 <- autoplot(pca.br, data = pData(eset.br), colour="ScanDate")
pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))
pl <- plot_grid(pl, pl3, nrow=2, align="hv")

# Save plot for manual quality control
save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA_nobatch.pdf", sep=""),
          pl, base_width=10, nrow=2)

## In case smth is wrong with PCA plots perform next steps.
## WARNING. You'll have to change some parts/reload data/save data again, cause it's always depends on the
## particular dataset.

# Batch-effect removal

# Eliminate samples, which are the single representation of particular batch. Reload data after that
# with new phenoData file, cause ExpressionSet subsetting works in a weird way
pdata = pData(eset.br)
n_occur <- data.frame(table(pdata$ScanDate))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
pdata <- pdata[-which(pdata$ScanDate %in% uniques),]
pd@data <- pd@data[pd@data$SampleAccessionNumber %in% pdata$SampleAccessionNumber,]
write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)

# Remove batch effect
batch = pData(eset.br)$ScanDate
mod = model.matrix(~as.factor(CancerType), data=pData(eset.br))
combat_edata = ComBat(dat=exprs(eset.br), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs(eset.br) <- combat_edata

combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs(eset) <- combat_edata

# In case not all the TNBC samples are clustered together
# perform manual step for outliers determination within TNBC subtype based on PCA plots
pd.tnbc <- pd[pd@data$CancerType=="TNBC"]
pca.tnbc <- pca.br$x[which(rownames(pca.br$x) %in% pd.tnbc@data$SampleAccessionNumber),]
outliers <- rownames(pca.tnbc[which(pca.tnbc[,1]<(37)),])
pd@data$Outliers <- rep("NA", length(pd@data$SampleAccessionNumber))
pd@data$Outliers[which(pd@data$SampleAccessionNumber %in% outliers)] <- "Yes"
autoplot(pca.br, data = pData(pd), colour="Outliers")
write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)
