# Loading libraries
library(lmtest)
library(stats)
library(cowplot)
library(plyr)
library(ggplot2)
# Loading custom scripts
source("scatterPlot.R")
source("diagnosticPlots.R")
source("barPlot.R")
source("partialCorPlot.R")

### Initial info section

# Load studies description file
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
# Define pipelines you want to analyze
pipe_types <- c("brainarray", "max", "maxoverall", "mean", "scores", "random")
#pipe_types <- c("brainarray", "bioconductor", "scores")
# Define, which subset of genes you want to represent on plots
plot_types <- c("all","onetoone","onetomany")
plot_type <- plot_types[3]
# Define plot names based on current plot type
if (plot_type=="all") {
  barplot.name = "../plots/all_barplot.pdf"
  scatterplot.name = "../plots/all_scatterplot.pdf"
  diagnosticplot.name = "_all_diagnosticplot_"
}
if (plot_type=="onetoone") {
  barplot.name = "../plots/onetoone_barplot.pdf"
  scatterplot.name = "../plots/onetoone_scatterplot.pdf"
  diagnosticplot.name = "_onetoone_diagnosticplot_"
}
if (plot_type=="onetomany") {
  barplot.name = "../plots/onetomany_barplot.pdf"
  scatterplot.name = "../plots/onetomany_scatterplot.pdf"
  diagnosticplot.name = "_onetomany_diagnosticplot_"
}

# Create classes for storing pipelines results
setClass("dataset",
         representation(ID="character", exprs="data.frame"))
setClass("pipeline",
         representation(name="character", dataset="dataset"))
pipelines <- c()

# Read exprs datasets into list of "pipeline" class elements
for (i in 1:length(studies$ID)) {
  for (pipe_type in pipe_types) {
    if (grepl("Illu", studies[i,]$platform) && pipe_type=="brainarray") {
    } else {
      if (pipe_type=="bioconductor") {
        exprs <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs_", "max", ".tsv", sep=""), header = TRUE, sep = "\t")
      } else {
        exprs <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs_", pipe_type, ".tsv", sep=""), header = TRUE, sep = "\t")
      }        
      if (plot_type=="onetomany") {
        otmGenes <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
        exprs <- exprs[rownames(exprs) %in% otmGenes$ENTREZID,]        
      }
      if (plot_type=="onetoone") {
        otoGenes <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
        exprs <- exprs[!(rownames(exprs) %in% otoGenes$ENTREZID),]        
      }
      dataset <- new("dataset", ID=as.character(studies[i,]$ID), exprs=exprs)
      pipelines <- c(pipelines, new("pipeline", name=pipe_type, dataset=dataset))
    }
  }
}

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

### Generating plots section

# Cycle through all the pipelines
for (i in 1:length(pipelines)) {
  # Just for simplicity of usage
  eset <- pipelines[[i]]@dataset@exprs$eset
  rnaseq <- pipelines[[i]]@dataset@exprs$rnaseq
  # Fit linear and cubic functions:  
  lin<-lm(rnaseq~eset)
  cub<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  
  # Individual scatter plot
  df <- data.frame(eset=eset, rnaseq=rnaseq)
  # Ugly constants, cause Illumina has really different X range
  # TODO: get max and min automatically
  if (grepl("Illu", studies[studies$ID==pipelines[[i]]@dataset@ID,]$platform)) {
    new_plot <- scatterPlot(df, paste(pipelines[[i]]@dataset@ID,pipelines[[i]]@name), lin, cub, c(6, 14, -7, 10))
  } else {
    new_plot <- scatterPlot(df, paste(pipelines[[i]]@dataset@ID,pipelines[[i]]@name), lin, cub, c(2, 14, -7, 10))
  }
  scatterplot_list <- c(scatterplot_list, list(new_plot))
  # Partial correlation plot. This was not included in the paper,
  # too unreliable because of almost random constants, see the partialCorPlot function
#  corPlot <- partialCorPlot(df, paste(pipelines[[i]]@dataset@ID,pipelines[[i]]@name))
#  corplot_list <- c(corplot_list, list(corPlot))

  # Commented this, because it is really computationally heavy, uncomment if needed
  ## Diagnostic plots
#  diagPlts<-diagPlot(lin)
#  pl1 <- plot_grid(plotlist=diagPlts, ncol=2, align="hv")
#  title <- ggdraw() + draw_label(paste("Linear Function Fit.", pipelines[[i]]@dataset@ID), fontface='bold')  
#  pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
#  diagPlts<-diagPlot(cub)
#  pl2 <- plot_grid(plotlist=diagPlts, ncol=2, align="hv")
#  title <- ggdraw() + draw_label(paste("Cubic Function Fit.", pipelines[[i]]@dataset@ID), fontface='bold')  
#  pl2 <- plot_grid(title, pl2, ncol=1, rel_heights=c(0.1, 1))
#  pl <- plot_grid(pl1, pl2, nrow=2, align="h")
#  save_plot(paste("../plots/", pipelines[[i]]@dataset@ID, diagnosticplot.name, pipelines[[i]]@name,".pdf", sep=""),
#            pl, base_width=5, nrow = 2)
  
  # Correlation data frame
  pearson[nrow(pearson)+1, ] <- c(pipelines[[i]]@name, round(cor(eset, rnaseq), 3), pipelines[[i]]@dataset@ID)
  spearman[nrow(spearman)+1, ] <- c(pipelines[[i]]@name, round(cor(eset, rnaseq, method="spearman"), 3), pipelines[[i]]@dataset@ID)
}

## Scatter plots
pl <- plot_grid(plotlist=scatterplot_list, align="v", ncol = length(pipe_types))
save_plot(scatterplot.name, pl,
         ncol = length(pipe_types),
         nrow = length(studies$ID),
         base_height=5,
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
pl <- plot_grid(pl, legend, rel_widths = c(3, .3))
# Save plot
save_plot(barplot.name, pl,
          base_height=7, base_width=15)