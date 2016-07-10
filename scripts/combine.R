# Load HuGene and hgu133 BrainArray packages
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv4.db)
library(WGCNA)
# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)


### Initial info
pipe_types <- c("brainarray", "mean", "random", "maxoverall", "scores")
pipe_type <- pipe_types[1]
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)

# Generate aggregated microarray-rnaseq expression files for plots feeding
# Didn't wrap it up in 'for (pipe_type in pipe_types)' cycle for testing purposes, but feel free
# to do it. Don't forget that there is no brainarray pipeline for Illumina
for (i in 1:length(studies$ID)) {
  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
  # Reading expression data
  if (grepl("Affy", studies[i,]$platform)) {
    platform <- "affymetrix"
  } else {
    platform <- "illumina"
  }
  if (pipe_type=="brainarray") {
    platform <- pipe_type
  }
  eset<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)
  if (!pipe_type=="scores") {
    # List of probesets IDs
    probesetsID<-rownames(eset)
    # List of corresponding gene IDs
    if (pipe_type=="brainarray") {
      probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
    } else if (studies[i,]$platformAbbr=="hugene10st") {
      probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
    } else {
      probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
    }
    # Replace probesetsIDs with gene IDs in expression data matrix
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
    # Duplicate probes
    n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
    # Filter eset based on left probesets
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID$PROBEID),]
  }
  # Collapse rows based on one of the algorithms
  if (pipe_type=="brainarray") {
    rownames(eset) <- probesetsID_EntrezID$ENTREZID
  }
  if (pipe_type=="maxoverall") {
    collapsed = collapseRows(eset, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
    eset <-collapsed$datETcollapsed
  }
  if (pipe_type=="mean") {
    collapsed = collapseRows(eset, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="Average")  
    eset <-collapsed$datETcollapsed
  }
  if (pipe_type=="random") {
    # Duplicate genes
    n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
    dupls <- n_occur[n_occur$Freq > 1,]$Var1
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    ind <- c()
    for (dupl in dupls) {
      ind <- c(ind, sample(which(probesetsID_EntrezID$ENTREZID==dupl),1))
    }
    probesetsID_EntrezID_random <- probesetsID_EntrezID[ind,]
    probesetsID_EntrezID_random <- rbind(probesetsID_EntrezID_random,
                                         probesetsID_EntrezID[probesetsID_EntrezID$ENTREZID %in% uniques,])
    
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID_random$PROBEID),]
    eset <- eset[order(match(rownames(eset),probesetsID_EntrezID_random$PROBEID)),]
    rownames(eset) <- probesetsID_EntrezID_random$ENTREZID
  }
  if (pipe_type=="scores") {
    probesetsID_EntrezID<-read.table(paste("../scores/", studies[i,]$platformAbbr, "_scores.tsv", sep=""), header=TRUE)
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
    n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    dupls <- n_occur[n_occur$Freq > 1,]$Var1
    ind <- c()
    for (dupl in dupls) {
      v <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID==dupl),]
      ind <- c(ind, as.numeric(rownames(v[which.max(v$Cov_Spec),])))
    }
    probesetsID_EntrezID_scores <- probesetsID_EntrezID[which(rownames(probesetsID_EntrezID) %in% ind),]
    probesetsID_EntrezID_scores <- rbind(probesetsID_EntrezID_scores,
                                         probesetsID_EntrezID[probesetsID_EntrezID$ENTREZID %in% uniques,])
    
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID_scores$PROBEID),]
    rownames(eset) <- probesetsID_EntrezID_scores[match(rownames(eset),probesetsID_EntrezID_scores$PROBEID),]$ENTREZID
    #eset <- eset[order(match(rownames(eset),probesetsID_EntrezID_scores$PROBEID)),]
    #rownames(eset) <- probesetsID_EntrezID_scores$ENTREZID
  }
  # Select only common genes for all the platforms and pipelines
  eset <- eset[rownames(eset) %in% genes_common$ENTREZID,]
  write.table(eset, paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}

# You need to generate all the datasets in allsamples_exprs folder before executing this
for (pipe_type in pipe_types) {
  eset.all <- data.frame()
  pdata.all <- data.frame()
  if (pipe_type=="brainarray") {
    len = length(studies$ID)-1
  } else {
    len = length(studies$ID)
  }
  for (i in 1:len) {
    pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
    eset<-read.table(paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_", pipe_type, ".tsv", sep=""), header=TRUE)
    if (nrow(eset.all)==0) {
      eset.all <- eset
      pdata.all <- pdata
    } else {
      # This subsetting is caused by scores pipelie
      if (nrow(eset.all)>nrow(eset)) {
        eset.all <- eset.all[which(rownames(eset.all) %in% rownames(eset)),]
      } else {
        eset <- eset[which(rownames(eset) %in% rownames(eset.all)),]
      }
      eset.all <- eset.all[order(match(rownames(eset.all),rownames(eset))),]
      eset.all <- cbind(eset.all, eset)
      pdata.all <- rbind(pdata.all, pdata)
    }
  }
  batch = pdata.all$DataSetAccesionNumber
  mod = model.matrix(~as.factor(CancerType), data=pdata.all)
  eset.all = ComBat(dat=eset.all, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  write.table(eset.all, paste("../combined/combined_expression_", pipe_type,".tsv", sep=""), sep="\t", quote=FALSE)
  write.table(pdata.all, paste("../pdata/combined_pdata.tsv", sep=""), sep="\t", quote=FALSE)
}

pdata <-read.table("../pdata/combined_pdata.tsv", header=TRUE, sep="\t")
for (pipe_type in pipe_types) {
  eset<-read.table(paste("../combined/combined_expression_", pipe_type,".tsv", sep=""), header=TRUE)
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  # Take mean across samples
  eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
  # Select common genes with rna-seq
  rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  rnaseq <- rnaseq[rnaseq$entrezID %in% rownames(eset),]
  rownames(rnaseq) <- rnaseq$entrezID
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(eset))),]
  rnaseq <- rnaseq[-1]
  # Assemble joint dataset
  exprs <- data.frame(row.names=rownames(eset), eset=eset$val, rnaseq=rnaseq$val)
  write.table(exprs, paste("../combined_exprs/combined_exprs_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}

# One-to-many genes here are union of one-to-many genes in all the platforms
otmGenes <- c()
for (i in 1:length(studies$ID)) {
  otm <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
  if (length(otmGenes)==0) {
    otmGenes <- as.vector(otm$ENTREZID)
  } else {
    otmGenes <- union(otmGenes, as.vector(otm$ENTREZID))
  }
}

# Create data frames for storing correlation results
pearson <- data.frame(pipeline=character(),
                      correlation=double(),
                      study=character(),  stringsAsFactors=FALSE) 
spearman <- data.frame(pipeline=character(),
                       correlation=double(),
                       study=character(),  stringsAsFactors=FALSE) 

source("scatterPlot.R")
scatterplot_list <- c()

plot_types <- c("all","onetoone","onetomany")
plot_type <- plot_types[3]
# Define plot names based on current plot type
if (plot_type=="all") {
  scatterplot.name = "../plots/combined_all_scatterplot.pdf"
}
if (plot_type=="onetoone") {
  scatterplot.name = "../plots/combined_onetoone_scatterplot.pdf"
}
if (plot_type=="onetomany") {
  scatterplot.name = "../plots/combined_onetomany_scatterplot.pdf"
}

for (pipe_type in pipe_types) {
  exprs <- read.table(paste("../combined_exprs/combined_exprs_", pipe_type, ".tsv", sep=""), sep="\t", header=TRUE)
  if (plot_type=="onetomany") {
    exprs <- exprs[rownames(exprs) %in% otmGenes,]        
  }
  if (plot_type=="onetoone") {
    exprs <- exprs[!(rownames(exprs) %in% otmGenes),]        
  }
  
  # Prepare scatterplot
  # Fit linear and cubic functions: 
  rnaseq <- exprs$rnaseq
  eset <- exprs$eset
  lin<-lm(rnaseq~eset)
  cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  # Individual scatter plot
  new_plot <- scatterPlot(exprs, paste(pipe_type, sep=""),
                          lin, cub, c(2, 14, -7, 10))
  
  scatterplot_list <- c(scatterplot_list, list(new_plot))
  
  # Correlation data frame
  pearson[nrow(pearson)+1, ] <- c(pipe_type, round(cor(eset, rnaseq), 3), plot_type)
  spearman[nrow(spearman)+1, ] <- c(pipe_type, round(cor(eset, rnaseq, method="spearman"), 3), plot_type)
}

## Scatter plots
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = length(pipe_types))
save_plot(scatterplot.name, pl,
          ncol = length(pipe_types),
          base_height=6,
          base_aspect_ratio=0.7)


## Bar plots
source("barPlot.R")
# Convert factor to numeric
pearson$correlation <- as.numeric(as.character(pearson$correlation))
spearman$correlation <- as.numeric(as.character(spearman$correlation))
# Pearson correlation bar plot
pl <- barPlot(pearson)
# Extract legend
g <- ggplotGrob(pl + theme(legend.position="left"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
title <- ggdraw() + draw_label("Pearson Correlation", fontface='bold')  
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
# Spearman correlation bar plot
pl1 <- barPlot(spearman)
title <- ggdraw() + draw_label("Spearman Correlation", fontface='bold')  
pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
pl <- plot_grid(pl, pl1, nrow=2, align="h")
# Combine two plots with the legend
pl <- plot_grid(pl, legend, rel_widths = c(3, .3))
# Save plot
save_plot("../plots/combined_barplot.pdf", pl,
          base_height=7, base_width=15)
