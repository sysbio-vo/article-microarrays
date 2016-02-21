library(lmtest)
library(stats)
library(cowplot)
library(plyr)
library(ggplot2)
source("scatterPlot.R")
source("diagnosticPlots.R")
### Initial info
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
#pipe_types <- c("brainarray", "max", "maxoverall", "mean", "scores", "random")
pipe_types <- c("brainarray", "mean", "scores")
onetomany = FALSE
onetoone = TRUE

# Create classes for storing pipelines results
setClass("dataset",
         representation(ID="character", exprs="data.frame"))
setClass("pipeline",
         representation(name="character", dataset="dataset"))
pipelines <- c()

for (i in 1:length(studies$ID)) {
  for (pipe_type in pipe_types) {
    if (grepl("Illu", studies[i,]$platform) && pipe_type=="brainarray") {
    } else {
      exprs <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs_", pipe_type, ".tsv", sep=""), header = TRUE, sep = "\t")
      if (onetomany) {
        otmGenes <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
        exprs <- exprs[rownames(exprs) %in% otmGenes$ENTREZID,]        
      }
      if (onetoone) {
        otoGenes <- read.table(paste("../general/onetomany_genes_", studies[i,]$platformAbbr, ".txt", sep=""), header = TRUE, sep = "\t")
        exprs <- exprs[!(rownames(exprs) %in% otoGenes$ENTREZID),]        
      }
      dataset <- new("dataset", ID=as.character(studies[i,]$ID), exprs=exprs)
      pipelines <- c(pipelines, new("pipeline", name=pipe_type, dataset=dataset))
    }
  }
}

pearson <- data.frame(pipeline=character(),
                      correlation=double(),
                      study=character(),  stringsAsFactors=FALSE) 
spearman <- data.frame(pipeline=character(),
                       correlation=double(),
                       study=character(),  stringsAsFactors=FALSE) 
plot_list = list()
for (ind in 1:length(pipelines)) {
  eset <- pipelines[[ind]]@dataset@exprs$eset
  rnaseq <- pipelines[[ind]]@dataset@exprs$rnaseq
  # Fit linear and polynomial functions:  
#   lin<-lm(rnaseq~eset)
#   quad<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  # Scatter plot
#   df <- data.frame(x=eset, y=rnaseq)
#   if (grepl("Illu", studies[studies$ID==pipelines[[ind]]@dataset@ID,]$platform)) {
#     new_plot <- scatterPlot(df, paste(pipelines[[ind]]@dataset@ID,pipelines[[ind]]@name), c(6, 14, -7, 10))
#   } else {
#     new_plot <- scatterPlot(df, paste(pipelines[[ind]]@dataset@ID,pipelines[[ind]]@name), c(2, 14, -7, 10))
#   }
#   plot_list <- c(plot_list, list(new_plot))
  # Diagnostic plots
#   diagPlts<-diagPlot(lin)
#   pl1 <- plot_grid(plotlist=diagPlts, ncol=2, align="hv")
#   title <- ggdraw() + draw_label(paste("Linear Function Fit.", pipelines[[ind]]@dataset@ID), fontface='bold')  
#   pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
#   diagPlts<-diagPlot(quad)
#   pl2 <- plot_grid(plotlist=diagPlts, ncol=2, align="hv")
#   title <- ggdraw() + draw_label(paste("Cubic Function Fit.", pipelines[[ind]]@dataset@ID), fontface='bold')  
#   pl2 <- plot_grid(title, pl2, ncol=1, rel_heights=c(0.1, 1))
#   pl <- plot_grid(pl1, pl2, nrow=2, align="h")
#   save_plot(paste("../plots/", pipelines[[ind]]@dataset@ID, "_onetoone_diagnostic_", pipelines[[ind]]@name,".pdf", sep=""),
#             pl, base_width=5, nrow = 2)
  # Correlation bar plots
  pearson[nrow(pearson)+1, ] <- c(pipelines[[ind]]@name, round(cor(eset, rnaseq), 3), pipelines[[ind]]@dataset@ID)
  spearman[nrow(spearman)+1, ] <- c(pipelines[[ind]]@name, round(cor(eset, rnaseq, method="spearman"), 3), pipelines[[ind]]@dataset@ID)
}

pearson$correlation <- as.numeric(as.character(pearson$correlation))
spearman$correlation <- as.numeric(as.character(spearman$correlation))

source("barPlot.R")
pl <- barPlot(pearson)

g <- ggplotGrob(pl + theme(legend.position="left"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

title <- ggdraw() + draw_label("Pearson Correlation", fontface='bold')  
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
pl1 <- barPlot(spearman)
title <- ggdraw() + draw_label("Spearman Correlation", fontface='bold')  
pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
pl <- plot_grid(pl, pl1, nrow=2, align="h")
pl <- plot_grid(pl, legend, rel_widths = c(3, .3))

plotname = ""
if (onetomany) {
  plotname = "../plots/onetomany_barplot.pdf"
} else if (onetoone) {
  plotname = "../plots/onetoone_barplot.pdf"
} else {
  plotname = "../plots/all_barplot.pdf"
}
save_plot(plotname, pl,
          base_height=7, base_width=15)


pl <- plot_grid(plotlist=plot_list, align="v", ncol = length(pipe_types))

plotname = ""
if (onetomany) {
  plotname = "../plots/onetomany_short_scatterplot.pdf"
} else if (onetoone) {
  plotname = "../plots/onetoone_short_scatterplot.pdf"
} else {
  plotname = "../plots/all_short_scatterplot.pdf"
}

save_plot(plotname, pl,
          ncol = length(pipe_types),
          nrow = length(studies$ID),
          base_height=5,
          base_aspect_ratio=0.7)


#dev.copy(pdf, paste("../plots/", pipe_type, "/", studies[affy[ind],]$ID, "_scatterplot_", pipe_type,".pdf", sep=""))
#dev.off()

# Test for linearity
pvalue <- anova(lin, quad)$"Pr(>F)"
stats:::plot.lm(lin, ask = FALSE, which=3)
stats:::plot.lm(quad, ask = FALSE, which=3)

# Another test for linearity
resettest(rnaseqExprs$val~eset$val, power=3, type="regressor")

# Breusch-Pagan test against heteroskedasticity
bptest(model1)

## From the car package
# Evaluate homoscedasticity
# non-constant error variance test
library(car)
ncvTest(model1)
# plot studentized residuals vs. fitted values
spreadLevelPlot(model1)
# Evaluate Nonlinearity
# component + residual plot
crPlots(model1)
# Ceres plots
ceresPlots(model2)

# Global test of model assumptions
library(gvlma)
gvmodel <- gvlma(quad)
summary(gvmodel)
