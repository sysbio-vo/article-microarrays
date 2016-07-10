library(cowplot)
library(plyr)
library(ggplot2)
library(lumi)
library(limma)
library(illuminaHumanv4.db)
library(WGCNA)
library(BioBase)
library(lumidat)
library(beadarray)
source("scatterPlot.R")
# General
studies <- read.table("../general/studies_idat.tsv", header = TRUE, sep = "\t")
i = 1
dir = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), sep="\t", head=TRUE)
idatfiles = dir(dir, pattern="idat", full.name=TRUE)
bgxfile = dir("../general", pattern="bgx", full.name=TRUE)
txtbgxfile = dir("../general", pattern="B.txt", full.name=TRUE)

# We've tried several packages, but BeadArray was the most easy to use, so scroll to BeadArray section

# Limma
raw.data <- read.idat(idatfiles, bgxfile)
N.data = neqc(raw.data)
probesetsID_EntrezID<-select(illuminaHumanv4.db, N.data$genes$Probe_Id, "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Subset EList
indices <- which(N.data$genes$Probe_Id %in% probesetsID_EntrezID$PROBEID)
N.data <- N.data[indices, ]
limma.data <- N.data
boxplot(limma.data$E)
datET <- as.data.frame(N.data$E, row.names=N.data$genes$Probe_Id)
datET <- datET[,which(pdata$Subs=="Yes")]

# Lumidat

raw.data <- lumiR.idat(idatfiles, probeID="ProbeID", manifestfile=txtbgxfile,
                    detectionTh=0.0001, QC=TRUE)
pData(raw.data) <- pdata
T.data   <- lumiT(raw.data, method="vst")
# The lumi package provides lumiB function for background correction. We sup-
# pose the BeadStudio output data has been background corrected. Therefore, no
# background correction used by default. A method ’bgAdjust’ is designed to ap-
# proximate what BeadStudio does for background adjustment. In the case when
# ’log2’ transform is used in the lumiT step, the background correction method
# (’forcePositive’) will be automatically used, which basically adds an offset (minus
# minimum value plus one) if there is any negative values to force all expression
# values to be positive.
N.data<-lumiN(T.data, method = "rsn")

probesetsID_EntrezID <- select(illuminaHumanv4.db, rownames(N.data), "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Subset EList
indices <- which(rownames(N.data) %in% probesetsID_EntrezID$PROBEID)
N.data <- exprs(N.data)[indices, ]
lumi.data <- N.data
boxplot(lumi.data)
datET <- N.data
datET <- datET[,which(pdata$Subs=="Yes")]

# BeadArray

raw.data <- readIdatFiles(idatfiles)
probesetsID<-rownames(exprs(raw.data))

probesetsID_EntrezID <- select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]

# Subset EList
methods <- data.frame(method=c("vsn", "quantile", "neqc"),
                    transform=c("none", "log2", "none"))
scatterplot_list <- c()
for (j in 1:nrow(methods)) {
  # Transform and normalize with one of the methods
  #N.data <- normaliseIllumina(channel(raw.data, "Green"),
                               method=as.character(methods[j,]$method), transform=as.character(methods[j,]$transform))
  # Subset the probesets and choose TNBC samples only
  #indices <- which(rownames(N.data) %in% probesetsID_EntrezID$PROBEID)
  #N.data <- exprs(N.data)[indices, ]
  #N.data <- N.data[,which(pdata$Subs=="Yes")]
  # Collapse rows
  #collapsed = collapseRows(N.data, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")
  # Intersect eset and rnaseq
  #eset <-collapsed$datETcollapsed
  #rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  #eset <- eset[rownames(eset) %in% rnaseq$entrezID,]
  #eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
  #rownames(rnaseq) <- rnaseq$entrezID
  #rnaseq <- rnaseq[which(rownames(rnaseq) %in% rownames(eset)),]
  #rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(eset))),]
  # Assemble joint dataset
  #exprs <- data.frame(row.names=rownames(eset), eset=eset$val, rnaseq=rnaseq$val)
  #write.table(exprs, paste("../exprs/", studies[i]$ID, "_exprs_max_", methods[j,]$method, ".tsv", sep=""), sep="\t", quote=FALSE)

   # If you need to re-render plots, just load the datasets saved before and comment all above
  exprs <- read.table(paste("../exprs/", studies[i]$ID, "_exprs_max_", methods[j,]$method, ".tsv", sep=""), sep="\t", header=TRUE)
  
  # Prepare scatterplot
  # Fit linear and cubic functions: 
  if (length(which(exprs$eset<0))>0)  
    exprs <- exprs[-which(exprs$eset<0),]
  exprs<-exprs[complete.cases(exprs),]
  rnaseq <- exprs$rnaseq
  eset <- exprs$eset
  lin<-lm(rnaseq~eset)
  cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  # Individual scatter plot
  new_plot <- scatterPlot(exprs, paste("Method: ", methods[j,]$method, ". Transform: ", methods[j,]$transform, sep=""),
                          lin, cub, c(5, 16, -8, 11))

  scatterplot_list <- c(scatterplot_list, list(new_plot))
}

## Scatter plots
title <- ggdraw() + draw_label(paste(studies[i,]$ID, "IDAT files processing. Max pipeline"),
                               fontface='bold', size=26)  
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = nrow(methods))
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

save_plot("../plots/qc/illumina_normalization_comparison.pdf", pl,
         ncol = nrow(methods),
         base_height=7,
         base_aspect_ratio=0.7)
