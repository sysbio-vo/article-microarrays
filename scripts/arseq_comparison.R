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
rnaseq <- read.table("../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_long.tsv", sep="\t")
rnaseq <- rnaseq[rownames(rnaseq) %in% common_genes$ENTREZID, ]

rnaseq.back <- rnaseq
arseq <- c()
for (i in 1:length(pipe_types)) {
  exprs <- read.table(paste("../allsamples_exprs/", studies[4,]$ID, "_allsamples_exprs_", pipe_types[i], ".tsv", sep=""), header=TRUE)
  rnaseq <- rnaseq.back

  exprs <- exprs[which(rownames(exprs) %in% rownames(rnaseq)),]
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(exprs))),]

  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  exprs <- exprs[colnames(exprs) %in% pdata$SampleAccessionNumber]
  rnaseq <- rnaseq[colnames(rnaseq) %in% pdata$SampleAccessionNumber]
  
  exprs <- data.frame(val=rowMeans(exprs), row.names=rownames(exprs))
  rnaseq <- data.frame(val=rowMeans(rnaseq), row.names=rownames(rnaseq))
  rnaseq$val <- rnaseq$val+0.1
  rnaseq$val <- log(rnaseq$val, 2)
  
  arseq <- c(arseq, list(data.frame(row.names=rownames(exprs), exprs=exprs$val, rnaseq=rnaseq$val)))
  
  #write.table(arseq, paste("../arseq/", studies[i,]$ID, "_arseq_self_rnaseq.tsv", sep=""), sep="\t", quote=FALSE)
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
  new_plot <- scatterPlot(df, pipe_types[i], lin, cub, c(6, 14, -5, 11))
  scatterplot_list <- c(scatterplot_list, list(new_plot))
  
  # Correlation data frame
  pearson[nrow(pearson)+1, ] <- c(pipe_types[i], round(cor(eset, rnaseq), 3), "GSE60788")
  spearman[nrow(spearman)+1, ] <- c(pipe_types[i], round(cor(eset, rnaseq, method="spearman"), 3), "GSE60788")
}

## Scatter plots
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = length(pipe_types))
save_plot("../plots/scatterPlots/arseq_all_scatterplot.pdf", pl,
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
save_plot("../plots/barPlots/arseq_all_barplot.pdf", pl,
          base_height=7, base_width=6)

## Self: rnaseq vs rnaseq

rnaseq.first <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_long.tsv", sep="\t", header=TRUE)
rnaseq.second <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_subread.tsv", sep="\t")

rnaseq.first <- rnaseq.first[rownames(rnaseq.first) %in% rownames(rnaseq.second), ]
rnaseq.second <- rnaseq.second[rownames(rnaseq.second) %in% rownames(rnaseq.first), ]

rnaseq.first <- rnaseq.first[rownames(rnaseq.first) %in% common_genes$ENTREZID, ]
rnaseq.second <- rnaseq.second[rownames(rnaseq.second) %in% common_genes$ENTREZID, ]

rnaseq.first <- rnaseq.first[order(match(rownames(rnaseq.first), rownames(rnaseq.second))),]

rnaseq.first <- data.frame(val=rowMeans(rnaseq.first), row.names=rownames(rnaseq.first))
rnaseq.second <- data.frame(val=rowMeans(rnaseq.second), row.names=rownames(rnaseq.second))
rnaseq.first$val <- rnaseq.first$val+0.1
rnaseq.first$val <- log(rnaseq.first$val, 2)

rnaseq <- rnaseq.first$val
eset <- rnaseq.second$val

lin<-lm(rnaseq~eset)
cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
df <- data.frame(eset=eset, rnaseq=rnaseq)

pl <- scatterPlot(df, "RNA-seq vs RNA-seq", lin, cub, c(-8, 12, -4, 12), labx="log2 CPM")
save_plot("../plots/scatterPlots/self_rnaseq_vs_rnaseq_scatterplot.pdf", pl,
          base_height=6,
          base_aspect_ratio=1)
