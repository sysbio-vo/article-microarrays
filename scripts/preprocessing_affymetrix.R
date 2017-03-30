# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
# Other packages
library(arrayQualityMetrics)
library(affy)
library(hgu133plus2cdf)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(cowplot)
source("plots_utils.R")

# Load studies description
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
affy <- which(grepl("Affy", studies$platform))
# Change this to index you need
i = affy[3]

READ.RAWS = FALSE

if (READ.RAWS) {
## Affymetrix data preprocessing
# Load pdata and CEL files
pd <- read.AnnotatedDataFrame(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""))
pd.plain <- read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header = TRUE, sep = "\t")
affyData = ReadAffy(phenoData=pd, sampleNames=pd$SampleAccessionNumber, filenames=as.character(rownames(pd)),
                    celfile.path=paste("../raws/", studies[i,]$ID, sep=""))

# Add ScanDate to pdata
pd.plain$ScanDate <- str_replace_all(affyData@protocolData@data$ScanDate, "T", " ")
pd.plain$ScanDate <- sapply(strsplit(pd.plain$ScanDate, split=' ', fixed=TRUE), function(x) (x[1]))
# Sometimes the scan date is in strange format, try this also
#pd.plain$ScanDate <- substr(pd.plain$ScanDate, 1, 7)

write.table(pd.plain, paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)

if (studies$platformAbbr[i]=="hgu133plus2") {
  # Affymetrix probesets
  affyData@cdfName <- "HG-U133_Plus_2"
  eset = rma(affyData)
  # Brainarray probesets
  affyData@cdfName <- "hgu133plus2hsentrezgcdf"
  eset.br = rma(affyData)
}

if (studies$platformAbbr[i]=="hugene10st") {
  # Affymetrix probesets
  affyData@cdfName <- "HuGene-1_0-st-v1"
  eset = rma(affyData)
  # Brainarray probesets. Change according to Affymetrix platform of the dataset
  affyData@cdfName <- "hugene10sthsentrezgcdf"
  eset.br = rma(affyData)
}  

# Save affy and brain expression sets
write.table(exprs(eset), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix_orig.tsv", sep=""),
            sep="\t", quote=FALSE)
write.table(exprs(eset.br), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray_orig.tsv", sep=""),
            sep="\t", quote=FALSE)

}

exprs <-    read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix_orig.tsv", sep=""),
                       sep="\t", header = TRUE)
exprs.br <- read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray_orig.tsv", sep=""),
                       sep="\t", header = TRUE)
pd.plain <- read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header = TRUE, sep = "\t")


## QC

# Perform PCA
pl <- qcPCA(exprs, exprs.br, pd.plain, studies[i,]$ID)

# Save plot for manual quality control
save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA.pdf", sep=""),
          pl, base_width=10, nrow=2)

# Batch-effect removal. Perform only if needed.
BATCH = FALSE
if (BATCH) {
  batch = pd.plain$ScanDate
  mod = model.matrix(~as.factor(CancerType), data=pd.plain)
  
  exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  exprs.br.nobatch = ComBat(dat=exprs.br, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

  pl <- qcPCA(exprs.nobatch, exprs.br.nobatch, pd.plain, studies[i,]$ID)
  
  # Save plot for manual quality control
  save_plot(paste("../plots/qc/", studies[i,]$ID, "_PCA_nobatch.pdf", sep=""),
            pl, base_width=10, nrow=2)
  
  # Save affy and brain expression sets after batch-effect removal
  write.table(exprs.nobatch, paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
  write.table(exprs.br.nobatch, paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), sep="\t", quote=FALSE)
}

# AQM to detect outliers. Edit pdata manually if you want to mark particular samples.
exprs.br <- read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""),
                       sep="\t", header = TRUE)
pd.plain <- read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header = TRUE, sep = "\t")
rownames(pd.plain) <- pd.plain$SampleAccessionNumber

eset = ExpressionSet(assayData=as.matrix(exprs.br), phenoData = AnnotatedDataFrame(pd.plain))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/aqm/AQM_report_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("CancerType"))
