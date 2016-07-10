# Load HuGene and hgu133 BrainArray packages
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv4.db)
library(reshape2)
library(ggplot2)
library(cowplot)
library(plyr)
library(stats)
library(Hmisc)

### Initial info
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
freq_list <- data.frame()
varplot_list <- c()
simple=FALSE
varofvar=TRUE
bootstrap=FALSE

# In general it looks into many-to-one probesets groups variance

# Generate aggregated microarray-rnaseq expression files for plots feeding
for (i in 1:length(studies$ID)) {
  # Reading expression data
  if (grepl("Affy", studies[i,]$platform)) {
    platform <- "affymetrix"
  } else {
    platform <- "illumina"
  }
  eset<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)

  # Reading phenodata. Select TNBC only
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  
  probesetsID<-rownames(eset)
  if (studies[i,]$platformAbbr=="hugene10st") {
    probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
  } else {
    probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
  }
  # Remove NA
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  # Duplicate probes
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  # Filter eset based on left probesets
  eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID$PROBEID),]

  # Analyze probesets mapped to one gene
  n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
  dupls <- n_occur[n_occur$Freq > 1,]
#   freq <- data.frame(table(dupls$Freq))
#   colnames(freq)[1] <- "ProbesetsNum"
#   freq$Study <- rep(studies[i,]$ID, nrow(freq))
#   freq_list <- structure(rbind(freq_list,freq), .Names = names(freq))

  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID %in% dupls$Var1),]
  eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID$PROBEID),]
  if (simple) {
    eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
    eset$ENTREZID <- probesetsID_EntrezID$ENTREZID
  }
  if (bootstrap) {
    eset$ENTREZID <- probesetsID_EntrezID$ENTREZID
    eset <- melt(eset)
    eset <- eset[,-2] 
  }
  if (varofvar) {
    eset$ENTREZID <- probesetsID_EntrezID$ENTREZID
  }

  rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  rnaseq <- rnaseq[rnaseq$entrezID %in% eset$ENTREZID,]
  eset <- eset[eset$ENTREZID %in% rnaseq$entrezID,]

  uniques <- unique(eset$ENTREZID)
  if (simple) {
    variance <- c()
    for (j in 1:length(uniques)) {
        variance <- c(variance, sd(eset[eset$ENTREZID %in% uniques[j],]$val))
    }
    eset <- data.frame(row.names=uniques, ENTREZID=uniques, variance=variance)
  }
  if (bootstrap) {
    variance <- data.frame(Mean=double(), Lower=double(), Upper=double(), stringsAsFactors=FALSE)
    for (j in 1:length(uniques)) {
      variance[nrow(variance)+1,] <- smean.cl.boot(eset[eset$ENTREZID %in% uniques[j],]$value)
    }
    eset <- data.frame(row.names=uniques, ENTREZID=uniques, variance=(variance$Upper - variance$Lower))
  }
  if (varofvar) {
    variance <- c()
    eset.temp <- eset[,-ncol(eset)]
    for (j in 1:length(uniques)) {
        variance <- c(variance, sd(apply(eset.temp[eset$ENTREZID %in% uniques[j],],2,sd)))
    }
    eset <- data.frame(row.names=uniques, ENTREZID=uniques, variance=variance)
  }

  exprs <- merge(eset, rnaseq, by.x = "ENTREZID", by.y="entrezID")
  colnames(exprs) <- c("ENREZID", "eset", "rnaseq")

  source("varianceScatterPlot.R")
  new_plot <- scatterPlot(exprs, paste(studies[i,]$ID, sep=""), c(-7, 11, 0, 2))
  varplot_list <- c(varplot_list, list(new_plot))  
}

pl <- plot_grid(plotlist=varplot_list, ncol=2, nrow=2, align="v")
save_plot("../plots/onetomany_varofvar_scatterplot.pdf", pl,
          base_height=4, base_aspect_ratio=1.5, nrow=2, ncol=2)

save_plot("../plots/onetomany_varofvar_scatterplot.pdf", new_plot,
          base_height=10)

# pl <- ggplot(data=freq_list, aes(x=ProbesetsNum, y=Freq, fill=ProbesetsNum)) +
#   geom_bar(stat="identity") +
#   facet_grid(~Study) + 
#   geom_text(aes(y=Freq, ymax=Freq, label=Freq),
#             position = position_dodge(width=0.9), vjust=-.5, color="black") + 
#   theme_gray() +
#   theme(legend.position = "none") +
#   labs(x = "Number of probesets mapped to one gene",
#        y = "Overall frequency")
# 
# save_plot("../plots/probesets_per_gene_frequency.pdf", pl,
#           base_height=8, base_aspect_ratio=2.1)


