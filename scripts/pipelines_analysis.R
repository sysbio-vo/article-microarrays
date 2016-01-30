library("MASS")
library(RColorBrewer)
library(lmtest)
library(stats)


### Initial info
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")
affy <- which(grepl("Affy", studies$platform))
pipe_type <- "brainarray"

# Create classes for storing pipelines results
setClass("dataset",
         representation(ID="character", exprs="data.frame"))
setClass("pipeline",
         representation(name="character", dataset="dataset"))
pipelines <- c()

for (i in affy) {
  exprs <- read.table(paste("../expression_data/", studies[affy[i],]$ID, "_exprs_", pipe_type, ".tsv", sep=""), header = TRUE, sep = "\t")
  dataset <- new("dataset", ID=as.character(studies[affy[i],]$ID), exprs=exprs)
  pipelines <- c(pipelines, new("pipeline", name=pipe_type, dataset=dataset))
}

ind = 1
eset <- pipelines[[ind]]@dataset@exprs$eset
rnaseq <- pipelines[[ind]]@dataset@exprs$rnaseq

# Fit linear and polynomial functions:
lin<-lm(rnaseq~eset)
quad<-lm(rnaseq~eset+I(eset^2)+I(eset^3))

# Scatter plot
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
z <- kde2d(eset, rnaseq, n=50)
plot(eset, rnaseq, xlab="log2 probeset intensity", ylab="log2 FPKM", main=paste(studies$ID[affy[ind]],pipe_type),pch="*", cex=.4)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
# Draw medians
abline(h=median(rnaseq), v=median(eset), lwd=1)
# Draw regression lines
r_lin <- summary(lin)$r.squared
r_quad <- summary(quad)$r.squared
index <- order(eset)
lines(eset[index], predict(lin)[index], col="red", lwd=2)
lines(eset[index], predict(quad)[index], col="blue", lwd=2)
n=min(eset)+1
m=max(eset)-1
text(n, 8, paste("R=",round(r_lin, digits = 3),sep="",collapse = ""),col="red",font=2)
text(m, 8, paste("R=",round(r_quad, digits = 3),sep="",collapse = ""),col="blue",font=2)

# Fit non-linear model
#loess_fit <- loess(rnaseqExprs$val~eset$val)
#lines(eset$val[index], predict(loess_fit)[index], col="blue")

dev.copy(pdf, paste("../plots/", pipe_type, "/", studies[affy[ind],]$ID, "_scatterplot_", pipe_type,".pdf", sep=""))
dev.off()

# ggplot2
library(ggplot2)
commonTheme = list(labs(x="log2 probeset intensity",
                        y="log2 FPKM"),
                   theme_bw()+theme(legend.position="none"))
ggplot(data=df,aes(x,y)) + 
  geom_point(alpha=0.3) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
  scale_fill_continuous(low="#2FFF71", high="#FF4C48") +
  geom_smooth(method=lm, linetype=1, colour="black", se=F, size=2) + 
  geom_smooth(method=lm, formula=y~x+I(x^2)+I(x^3), linetype=1, colour="blue", se=F, size=2) + 
  commonTheme


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
gvmodel <- gvlma(model1)
summary(gvmodel)
