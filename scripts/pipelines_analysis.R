library(lmtest)
library(stats)
library(cowplot)
library(plyr)
library(ggplot2)
source("scatterPlot.R")
source("diagnosticPlots.R")
### Initial info
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
affy <- which(grepl("Affy", studies$platform))
pipe_types <- c("brainarray", "max", "mean", "random", "maxoverall")

# Create classes for storing pipelines results
setClass("dataset",
         representation(ID="character", exprs="data.frame"))
setClass("pipeline",
         representation(name="character", dataset="dataset"))
pipelines <- c()

for (i in 1:length(studies$ID)) {
  for (pipe_type in pipe_types) {
    exprs <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs_", pipe_type, ".tsv", sep=""), header = TRUE, sep = "\t")
    dataset <- new("dataset", ID=as.character(studies[i,]$ID), exprs=exprs)
    pipelines <- c(pipelines, new("pipeline", name=pipe_type, dataset=dataset))
  }
}

plot_list = list()
for (ind in 1:length(pipelines)) {
  eset <- pipelines[[ind]]@dataset@exprs$eset
  rnaseq <- pipelines[[ind]]@dataset@exprs$rnaseq
  # Fit linear and polynomial functions:  
  lin<-lm(rnaseq~eset)
  quad<-lm(rnaseq~eset+I(eset^2)+I(eset^3))  
  # Scatter plot
  df <- data.frame(x=eset, y=rnaseq)
  new_plot <- scatterPlot(df, paste(pipelines[[ind]]@dataset@ID,pipelines[[ind]]@name))
  plot_list <- c(plot_list, list(new_plot))
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
#   save_plot(paste("../plots/", pipelines[[ind]]@dataset@ID, "_diagnostic_", pipelines[[ind]]@name,".png", sep=""),
#             pl, base_height=8, nrow = 2)
}

pl <- plot_grid(plotlist=plot_list, ncol=5, align="hv")
save_plot(paste("../plots/all_scatterplot.png", sep=""), pl,
          ncol = 5,
          nrow = 3,
          base_height=6,
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
