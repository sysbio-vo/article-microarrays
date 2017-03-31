library(org.Hs.eg.db)
source("plots_utils.R")

pipe_types <- c("maxoverall", "mean", "random", "scores")
pipe_type <- pipe_types[2]
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
common_genes <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)

arseq <- c()
for (i in 1:length(pipe_types)) {
  pdata <- read.table(paste("../pdata/pdata_", studies[4,]$ID, ".tsv", sep=""), header=TRUE, sep="\t", stringsAsFactors = FALSE)
  exprs <- read.table(paste("../allsamples_exprs/", studies[4,]$ID, "_allsamples_exprs_", pipe_types[i], ".tsv", sep=""), header=TRUE)
  
  rnaseq <- read.table("../rnaseq/GSE60788/GSE60788_rnaseq_gex_normalized.tsv", header=TRUE, sep="\t")
  rnaseq <- rnaseq[, -which(grepl("replicate", colnames(rnaseq)))]
  colnames(rnaseq) <- c("SYMBOL", pdata$SampleAccessionNumber)
  
  SYMBOLs <- rownames(rnaseq)
  EntrezID_SYMBOLs <- select(org.Hs.eg.db, as.character(common_genes$ENTREZID), c("SYMBOL"))
  rnaseq <- rnaseq[which(rnaseq$SYMBOL %in% EntrezID_SYMBOLs$SYMBOL),]
  
  rownames(rnaseq) <- EntrezID_SYMBOLs[match(rnaseq$SYMBOL, EntrezID_SYMBOLs$SYMBOL),]$ENTREZID
  
  exprs <- exprs[which(rownames(exprs) %in% rownames(rnaseq)),]
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(exprs))),]
  rnaseq <- rnaseq[-1]
  
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  exprs <- exprs[colnames(exprs) %in% pdata$SampleAccessionNumber]
  rnaseq <- rnaseq[colnames(rnaseq) %in% pdata$SampleAccessionNumber]
  
  exprs <- data.frame(val=rowMeans(exprs), row.names=rownames(exprs))
  rnaseq <- data.frame(val=rowMeans(rnaseq), row.names=rownames(rnaseq))
  rnaseq$val <- rnaseq$val+5
  
  arseq <- c(arseq, list(data.frame(row.names=rownames(exprs), exprs=exprs$val, rnaseq=rnaseq$val)))
  
  #write.table(arseq, paste("../arseq/", studies[i,]$ID, "_arseq_self_rnaseq.tsv", sep=""), sep="\t", quote=FALSE)
}

## Scatter plots
scatterplot_list = list()

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
  new_plot <- scatterPlot(df, pipe_types[i], lin, cub, c(6, 14, -7, 10))
  scatterplot_list <- c(scatterplot_list, list(new_plot))
  
  # Correlation data frame
  #pearson[nrow(pearson)+1, ] <- c(pipelines[[i]]@name, round(cor(eset, rnaseq), 3), pipelines[[i]]@dataset@ID)
  #spearman[nrow(spearman)+1, ] <- c(pipelines[[i]]@name, round(cor(eset, rnaseq, method="spearman"), 3), pipelines[[i]]@dataset@ID)

}

## Scatter plots
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = length(pipe_types))
save_plot("../plots/scatterPlots/arseq_all_scatterplot.pdf", pl,
          ncol = length(pipe_types),
          nrow = 1,
          base_height=5,
          base_aspect_ratio=0.7)
