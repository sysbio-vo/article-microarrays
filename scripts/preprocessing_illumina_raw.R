library(cowplot)
library(plyr)
library(ggplot2)
library(limma)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(WGCNA)
library(beadarray)
library(sva)
source("scatterPlot.R")
source("plots_utils.R")
source("probes_utils.R")

## General
studies <- read.table("../general/studies_idat.tsv", header = TRUE, sep = "\t")
i = 1
dir = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), sep="\t", head=TRUE)
idatfiles = dir(dir, pattern="idat", full.name=TRUE)
bgxfile = dir("../general", pattern="bgx", full.name=TRUE)
txtbgxfile = dir("../general", pattern="B.txt", full.name=TRUE)

## BeadArray
raw.data <- readIdatFiles(idatfiles)
raw.data@phenoData = AnnotatedDataFrame(pdata)

## QC
N.data <- normaliseIllumina(channel(raw.data, "Green"), method="neqc", transform="none")
N.exprs <- exprs(N.data)
pdata$Sentrix_ID <- as.character(pdata$Sentrix_ID)

pl <- qcPCAilluraw(N.exprs, pdata, studies[i,]$ID)
# Save plot for manual quality control
save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA.pdf", sep=""),
          pl, base_width=10)

# Correct batch
batch = pdata$Sentrix_ID
mod = model.matrix(~as.factor(CancerType), data=pdata)
N.exprs.nobatch = ComBat(dat=N.exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

pl <- qcPCAilluraw(N.exprs.nobatch, pdata, studies[i,]$ID)
# Save plot for manual quality control
save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA_nobatch.pdf", sep=""),
          pl, base_width=10)

exl = which(pdata$Outliers=="Yes")
pdata <- pdata[-exl,]
N.exprs.nobatch <- N.exprs.nobatch[,-exl]

eset = ExpressionSet(assayData=as.matrix(N.exprs.nobatch), phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/aqm/AQM_report_nooutliers_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("CancerType"))

# Get EntrezID
probesetsID <- rownames(exprs(raw.data))
probesetsID_EntrezID <- select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Remove probes mapped to different genes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]

# Normalization methods
methods <- data.frame(method=c("vsn", "quantile", "neqc"),
                    transform=c("none", "log2", "none"))

scatterplot_list <- c()
for (j in 1:nrow(methods)) {
  # Transform and normalize with one of the methods
  #N.data <- normaliseIllumina(channel(raw.data, "Green"),
  #                             method=as.character(methods[j,]$method), transform=as.character(methods[j,]$transform))
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
