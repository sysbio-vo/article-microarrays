# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
# Load Bioconductor annotation packages
library(hgu133plus2.db)
library(lumiHumanAll.db)
# Other libs
library(affy)

### Initial info
studies <- read.table("studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("expression_data/common_genes.txt", header = FALSE, sep = "\t")

### RNA-seq data
rnaseqExprs <- read.table("expression_data/rnaseq_data_processed_sum.tsv", sep="\t", header=TRUE)

### Brainarray
pipe_type <- "brainarray"
ind <- which(grepl("Affy", studies$platform))
# Reading phenodata
pdata <-read.table(paste("pdata/pdata_", studies[1,]$ID, ".txt", sep=""), header=TRUE, sep="\t")
# Reading expression data
eset<-read.table(paste("expression_data/", studies[1,]$ID, "_", pipe_type, ".txt", sep=""), header=TRUE)
# List of probesets IDs
probesetsID<-rownames(eset)
# List of corresponding gene IDs
detach("package:MASS", unload=TRUE)
probesetsID_EntrezID<-select(get(paste(studies[1,]$platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
# Replace probesetsIDs with gene IDs in expression data matrix
eset <- eset[which(probesetsID_EntrezID$ENTREZID!="NA"),]
rownames(eset) <- probesetsID_EntrezID$ENTREZID[which(probesetsID_EntrezID$ENTREZID!="NA")]
# Select only common genes for all the platforms and pipelines
eset <- eset[rownames(eset) %in% genes_common[,1],]
# Select TNBC samples only
pdata<-pdata[which(pdata$CancerType=="TNBC"),]
eset <- eset[colnames(eset) %in% pdata$SampleAccessionNumber]
# Select common genes with rna-seq
eset <- eset[rownames(eset) %in% rnaseqExprs$entrezID,]
eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
# Organize rna-seq dataset
rownames(rnaseqExprs) <- rnaseqExprs$entrezID
rnaseqExprs <- rnaseqExprs[-1]
rnaseqExprs$val <- rnaseqExprs[order(match(rownames(rnaseqExprs),rownames(eset))),]
rownames(rnaseqExprs) <- rownames(eset)
# Scatter plot
library("MASS")
library(RColorBrewer)
library(lmtest)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
z <- kde2d(eset$val, rnaseqExprs$val, n=50)
plot(eset$val,rnaseqExprs$val, xlab="Probeset intensity", ylab="log2 FPKM", main=paste(studies$ID[1],pipe_type),pch="*", cex=.4)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
abline(h=median(rnaseqExprs$val), v=median(eset$val), lwd=1)
# Fit linear model:
lm_exprs<-lm(rnaseqExprs$val~eset$val)
sumry_lm<-summary(lm_exprs)
value_add<-sumry_lm$r.squared
abline(lm_exprs,col="red")
n=min(eset$val)+1
text(n,8,paste("R=",round(value_add,digits = 3),sep="",collapse = ""),col="red",font=2)
# Non-linear fits
loess_fit <- loess(rnaseqExprs$val~eset$val)
nls_fit <- nls(rnaseqExprs$val ~ a + b * eset$val^(-c), start = list(a = 80, b = 20, c = 0.2))
lines(predict(nls_fit), col="blue")

dev.copy(pdf,paste(studies[1,]$ID, '_scatter_brain_sum_loess.pdf', sep=""))
dev.off()

# Test for linearity
model1<-lm(rnaseqExprs$val~eset$val)
model2<-lm(rnaseqExprs$val~eset$val+I(rnaseqExprs$val)^2)
anova(model1,model2)

# Another test for linearity
resettest(rnaseqExprs$val~eset$val, power=2, type="regressor")

# White's test for heteroscedisity
bptest(lm_exprs)
