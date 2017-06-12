library(cowplot)
library(plyr)
library(ggplot2)
library(ggfortify)
# Loading custom scripts
source("plots_utils.R")

pipe_types <- c("brainarray", "maxoverall", "mean", "scores", "random")
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")

# Define rnaseq type
#seq_type = "subread"; ymin = -7; laby="log2 CPM"; dir="../arseq_subread_merged/arseq_subread_merged_"
seq_type = "tophat"; ymin = -7; laby="log2 FPKM"; dir="../arseq_tophat_merged/arseq_tophat_merged_"

otmGenes <- c()
for (i in 1:length(studies$ID)) {
  otm <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
  if (length(otmGenes)==0) {
    otmGenes <- as.vector(otm$ENTREZID)
  } else {
    otmGenes <- union(otmGenes, as.vector(otm$ENTREZID))
  }
}

plot_types <- c("all","onetoone","onetomany")
plot_type <- plot_types[3]

for (plot_type in plot_types) {
  
# Create data frames for storing correlation results
pearson <- data.frame(pipeline=character(),
                      correlation=double(),
                      study=character(),  stringsAsFactors=FALSE) 
spearman <- data.frame(pipeline=character(),
                       correlation=double(),
                       study=character(),  stringsAsFactors=FALSE) 

# Define plot names based on current plot type
if (plot_type=="all") {
  scatterplot.name = paste("../plots/scatterPlots/merged_all_scatterplot_", seq_type, ".pdf", sep="")
  barplot.name = paste("../plots/barPlots/merged_all_barplot_", seq_type, ".pdf", sep="")
}
if (plot_type=="onetoone") {
  scatterplot.name = paste("../plots/scatterPlots/merged_onetoone_scatterplot_", seq_type, ".pdf", sep="")
  barplot.name = paste("../plots/barPlots/merged_onetoone_barplot_", seq_type, ".pdf", sep="")
}
if (plot_type=="onetomany") {
  scatterplot.name = paste("../plots/scatterPlots/merged_onetomany_scatterplot_", seq_type, ".pdf", sep="")
  barplot.name = paste("../plots/barPlots/merged_onetomany_barplot_", seq_type, ".pdf", sep="")
}

scatterplot_list <- c()
for (pipe_type in pipe_types) {
  exprs <- read.table(paste(dir, pipe_type, ".tsv", sep=""), sep="\t", header=TRUE)
  if (plot_type=="onetomany") {
    exprs <- exprs[rownames(exprs) %in% otmGenes,]        
  }
  if (plot_type=="onetoone") {
    exprs <- exprs[!(rownames(exprs) %in% otmGenes),]        
  }
  
  # Prepare scatterplot
  # Fit linear and cubic functions: 
  rnaseq <- exprs$rnaseq
  eset <- exprs$exprs
  lin<-lm(rnaseq~eset)
  cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  df <- data.frame(eset=eset, rnaseq=rnaseq)
  # Individual scatter plot
  new_plot <- scatterPlot(df, paste(pipe_type, sep=""),
                          lin, cub, c(2, 14, -7, 11), laby=laby)
  
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
pl <- plot_grid(pl, legend, rel_widths = c(3, .7))
# Save plot
save_plot(barplot.name, pl, base_height=7, base_width=7)

}

#PCA
pipe_type = "maxoverall"

pdata <-read.table("../pdata/combined_pdata.tsv", header=TRUE, sep="\t")
eset<-read.table(paste("../allsamples_exprs_merged/allsamples_exprs_merged_", pipe_type,".tsv", sep=""), header=TRUE)
pdata$Grade <- as.character(pdata$Grade)
pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]

pca = prcomp(t(eset))

title <- ggdraw() + draw_label("PCA plot of combined microarray dataset. Max pipeline", fontface='bold')  
pl1 <- autoplot(pca, data = pdata, colour="CancerType")
pl2 <- autoplot(pca, data = pdata, colour="Grade")
pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

# pdata.grade<-pdata[which(pdata$Grade!="NA"),]
# eset.grade <- eset[,colnames(eset) %in% pdata.grade$SampleAccessionNumber]
# pca.grade = prcomp(t(eset.grade))

# title <- ggdraw() + draw_label("Samples with grade present", fontface='bold')  
# pl1 <- autoplot(pca.grade, data = pdata.grade, colour="CancerType")
# pl2 <- autoplot(pca.grade, data = pdata.grade, colour="Grade")
# pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
# pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))

# pl <- plot_grid(pl, pl3, nrow=2)

# Save plot for manual quality control
save_plot(paste("../plots/merged_PCA.pdf", sep=""),
          pl, base_height=5,
          base_aspect_ratio=2.3,
          nrow=1)

# In case one is interested in per grade representation
# pdata.grade3 <-pdata[which(pdata$Grade=="3"),]
# eset.grade3 <- eset[,colnames(eset) %in% pdata.grade3$SampleAccessionNumber]
# pca.grade3 = prcomp(t(eset.grade3))
# pdata.grade2 <-pdata[which(pdata$Grade=="2"),]
# eset.grade2 <- eset[,colnames(eset) %in% pdata.grade2$SampleAccessionNumber]
# pca.grade2 = prcomp(t(eset.grade2))
# pdata.grade1 <-pdata[which(pdata$Grade=="1"),]
# eset.grade1 <- eset[,colnames(eset) %in% pdata.grade1$SampleAccessionNumber]
# pca.grade1 = prcomp(t(eset.grade1))
# 
# 
# title <- ggdraw() + draw_label("Samples with grade=1", fontface='bold')  
# pl1 <- autoplot(pca.grade1, data = pdata.grade1, colour="CancerType")
# pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
# title <- ggdraw() + draw_label("Samples with grade=2", fontface='bold')  
# pl2 <- autoplot(pca.grade2, data = pdata.grade2, colour="CancerType")
# pl2 <- plot_grid(title, pl2, ncol=1, rel_heights=c(0.1, 1))
# title <- ggdraw() + draw_label("Samples with grade=3", fontface='bold')  
# pl3 <- autoplot(pca.grade3, data = pdata.grade3, colour="CancerType")
# pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))
# 
# pl4 <- plot_grid(pl1, pl2, pl3, ncol=3)
# 
# # Save plot for manual quality control
# save_plot(paste("../plots/combined_PCA_pergrade_br.pdf", sep=""),
#           pl4, base_height=5,
#           ncol=3)

