library(org.Hs.eg.db)
library(WGCNA)
library(ggplot2)
library(cowplot)
library(plyr)
source("plots_utils.R")
source("probes_utils.R")

pipe_types <- c("maxoverall", "mean", "random", "scores")
pipe_type <- pipe_types[2]
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
common_genes <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)

pdata <- read.table(paste("../pdata/pdata_", studies[4,]$ID, ".tsv", sep=""), header=TRUE, sep="\t", stringsAsFactors = FALSE)
rnaseq <- read.table("../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_short.tsv", sep="\t")
rnaseq <- rnaseq[rownames(rnaseq) %in% common_genes$ENTREZID, , drop=FALSE]

rnaseq.back <- rnaseq
arseq <- c()
for (i in 1:length(pipe_types)) {
  exprs <- read.table(paste("../allsamples_exprs/", studies[4,]$ID, "_allsamples_exprs_", pipe_types[i], ".tsv", sep=""), header=TRUE)
  exprs <- exprs[which(rownames(exprs) %in% rownames(rnaseq)),]
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(exprs))), , drop=FALSE]

  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  exprs <- exprs[colnames(exprs) %in% pdata$SampleAccessionNumber]

  exprs <- data.frame(val=rowMeans(exprs), row.names=rownames(exprs))
  arseq <- c(arseq, list(data.frame(row.names=rownames(exprs), exprs=exprs$val, rnaseq=rnaseq$val)))
}

## Scatter and bar plots
# Create data frames for storing correlation results
pearson <- data.frame(pipeline=character(),
                      correlation=double(),
                      study=character(),  stringsAsFactors=FALSE) 
spearman <- data.frame(pipeline=character(),
                       correlation=double(),
                       study=character(),  stringsAsFactors=FALSE) 

# Create empty list for storing indibidual scatter plots
scatterplot_list = list()
corplot_list = list()

for (i in 1:length(arseq)) {
  # Just for simplicity of usage
  eset <- arseq[[i]]$exprs
  rnaseq <- arseq[[i]]$rnaseq
  # Fit linear and cubic functions:  
  lin<-lm(rnaseq~eset)
  cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  
  # Individual scatter plot
  df <- data.frame(eset=eset, rnaseq=rnaseq)
  # Ugly constants, cause Illumina has really different X range
  # TODO: get max and min automatically
  new_plot <- scatterPlot(df, pipe_types[i], lin, cub, c(5, 14, -7, 12))
  scatterplot_list <- c(scatterplot_list, list(new_plot))
  
  # Correlation data frame
  pearson[nrow(pearson)+1, ] <- c(pipe_types[i], round(cor(eset, rnaseq), 3), "GSE60788")
  spearman[nrow(spearman)+1, ] <- c(pipe_types[i], round(cor(eset, rnaseq, method="spearman"), 3), "GSE60788")
}

## Scatter plots
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = length(pipe_types))
save_plot("../plots/scatterPlots/GSE60788_arseq_all_scatterplot.pdf", pl,
          ncol = length(pipe_types),
          nrow = 1,
          base_height=5,
          base_aspect_ratio=0.6)

## Bar plots
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
pl <- plot_grid(pl, legend, ncol=2, rel_widths = c(2, 1))
# Save plot
save_plot("../plots/barPlots/GSE60788_arseq_all_barplot.pdf", pl,
          base_height=7, base_width=6)

## Self: rnaseq vs rnaseq

rnaseq.first <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_short.tsv", sep="\t", header=TRUE)
rnaseq.second <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_subread_short.tsv", sep="\t")

rnaseq.first <- rnaseq.first[rownames(rnaseq.first) %in% common_genes$ENTREZID, , drop=FALSE]
rnaseq.second <- rnaseq.second[rownames(rnaseq.second) %in% common_genes$ENTREZID, , drop=FALSE]

rnaseq.first <- rnaseq.first[order(match(rownames(rnaseq.first), rownames(rnaseq.second))), , drop=FALSE]

rnaseq <- rnaseq.first$val
eset <- rnaseq.second$val

lin<-lm(rnaseq~eset)
cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
df <- data.frame(eset=eset, rnaseq=rnaseq)

pl <- scatterPlot(df, "GSE58135 RNA-seq: tophat vs subread", lin, cub, c(-7, 11, -7, 11), labx="log2 CPM")
save_plot("../plots/scatterPlots/GSE58135_tophat_vs_subread_scatterplot_new.pdf", pl,
          base_height=6,
          base_aspect_ratio=1)


######

rnaseq.first <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_short.tsv", sep="\t", header=TRUE)
rnaseq.second <- read.table("../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_short.tsv", sep="\t")

rnaseq.second <- rnaseq.second[rownames(rnaseq.second) %in% common_genes$ENTREZID, , drop=FALSE]
rnaseq.first <- rnaseq.first[rownames(rnaseq.first) %in% rownames(rnaseq.second), , drop=FALSE]

rnaseq.first <- rnaseq.first[order(match(rownames(rnaseq.first), rownames(rnaseq.second))), ,drop=FALSE]

rnaseq <- rnaseq.first$val
eset <- rnaseq.second$val

lin<-lm(rnaseq~eset)
cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
df <- data.frame(eset=eset, rnaseq=rnaseq)

pl <- scatterPlot(df, "GSE58135 vs GSE60788 RNA-seq (tophat)", lin, cub, c(-7, 12, -7, 11))
save_plot("../plots/scatterPlots/GSE58135_vs_GSE60788_rnaseq_tophat_scatterplot.pdf", pl,
          base_height=6,
          base_aspect_ratio=1)
