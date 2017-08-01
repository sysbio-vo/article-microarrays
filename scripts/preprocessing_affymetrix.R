# Load HuGene and hgu133 BrainArray packages
source("https://bioconductor.org/biocLite.R")


library(affycoretools)

# library(hugene10sthsentrezgprobe)
# library(hugene10sthsentrezgcdf)
# library(hugene10sthsentrezg.db)
#
# library(hgu133plus2hsentrezgprobe)
# library(hgu133plus2hsentrezgcdf)
# library(hgu133plus2hsentrezg.db)
# Other packages
library(affy)
library(sva)
library(stringr)
library(ggplot2)
library(ggfortify)
library(cowplot)
library(ArrayExpress)

setwd('/home/rstudio/r/article-microarrays')
getwd()

rawspath = 'raws/affymetrix'
prepath = 'preprocessed/affymetrix'
pdatapath = 'pdata/'
plotsqcpath = paste(getwd(), 'plots/qc/', sep='/')


# Load studies description
studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")

# load IGEA phenodata
igea = read.table('igea_tsv/samples.tsv',header = TRUE, sep = '\t', fill = TRUE)

# install cdf annotation files for all listed microarray platforms
for (array in levels(studies$platformAbbr)){
  install.brainarray(array)
}

i = 4

current_path = paste(rawspath, '/', studies$accession[[i]], sep='')
if (! dir.exists(current_path)){
  dir.create(current_path)
}


aeData = getAE(
  studies$accession[[i]],
  path = current_path,
  sourcedir=current_path,
  # local = TRUE,
  type = 'raw')



z <- ArrayExpress:::readPhenoData(aeData$sdrf, aeData$path)


# merge ArrayExpress phenodata with IGEA phenodata
pd = merge(z@data, igea, all.x = TRUE, by.x = 'Source.Name', by.y = 'Sample.Name')

rownames(pd) = pd$Array.Data.File

nrow(pd)


pd$Experiment


affyData = ReadAffy(phenoData=pd,
                           sampleNames=pd$Sample.Name,
                           filenames=pd$Array.Data.File,
                           celfile.path=paste("raws/affymetrix/",
                                              pd["Experiment"][[1]],
                                              sep="")
)


affyData@cdfName <- paste(studies$platformAbbr[[i]], 'hsentrezgcdf', sep="")
nrow(affyData)

affyData.rma = rma(affyData)
nrow(exprs(affyData.rma))


# Save affy and brain expression sets
write.table(exprs(affyData.rma), paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
t = read.table(paste(prepath, '/', studies$accession[[i]], "_preprocessed_affymetrix.tsv", sep=""), sep="\t")
t


# QC
# Add ScanDate to pdata
pData(eset)$ScanDate <- str_replace_all(eset@protocolData@data$ScanDate, "T", " ")
pData(eset)$ScanDate <- sapply(strsplit(pData(eset)$ScanDate, split=' ', fixed=TRUE), function(x) (x[1]))
pData(eset.br)$ScanDate <- pData(eset)$ScanDate

# Sometimes the scan date is in strange format, try this also
#pData(eset.br)$ScanDate <- substr(pData(eset.br)$ScanDate, 1, 7)
#pData(eset)$ScanDate <- substr(pData(eset)$ScanDate, 1, 7)

# Perform PCA
pca = prcomp(t(exprs(eset)))
pca.br = prcomp(t(exprs(eset.br)))

title <- ggdraw() + draw_label(paste("Affymetrix Probesets Definitions.", studies[i,]$ID), fontface='bold')

pl1 <- autoplot(pca, data = pData(eset), colour="Characteristics.condition.")
pl2 <- autoplot(pca, data = pData(eset), colour="ScanDate")
pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label(paste("Brainarray Probesets Definitions.", studies[i,]$ID), fontface='bold')

pl1 <- autoplot(pca.br, data = pData(eset.br), colour="Characteristics.condition.")
pl2 <- autoplot(pca.br, data = pData(eset.br), colour="ScanDate")
pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))
pl <- plot_grid(pl, pl3, nrow=2, align="hv")

# Save plot for manual quality control
save_plot(paste(plotsqcpath, studies[i,]$accession, "_PCA_nobatch.pdf", sep=""),
          pl1, base_width=10, nrow=2)


# Batch-effect removal. Perform only if needed.

# Eliminate samples, which are the single representation of particular batch. Reload data after that
# with new phenoData file, cause ExpressionSet subsetting works in a weird way
pdata = pData(eset.br)
n_occur <- data.frame(table(pdata$ScanDate))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
pdata <- pdata[-which(pdata$ScanDate %in% uniques),]
pd@data <- pd@data[pd@data$SampleAccessionNumber %in% pdata$SampleAccessionNumber,]
write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)

# Remove batch effect
batch = pData(eset.br)$ScanDate
mod = model.matrix(~as.factor(CancerType), data=pData(eset.br))
combat_edata = ComBat(dat=exprs(eset.br), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs(eset.br) <- combat_edata

combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs(eset) <- combat_edata

# Save affy and brain expression sets after batch-effect removal
write.table(exprs(eset), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_affymetrix.tsv", sep=""), sep="\t", quote=FALSE)
write.table(exprs(eset.br), paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), sep="\t", quote=FALSE)

# In case not all the TNBC samples are clustered together
# perform manual step for outliers determination within TNBC subtype based on PCA plots
# However, pdata files already contain Outlier variable, therefor you can skip this and use
# ready pdata files instead
pd.tnbc <- pd[pd@data$CancerType=="TNBC"]
pca.tnbc <- pca.br$x[which(rownames(pca.br$x) %in% pd.tnbc@data$SampleAccessionNumber),]
outliers <- rownames(pca.tnbc[which(pca.tnbc[,1]<(37)),])
pd@data$Outliers <- rep("NA", length(pd@data$SampleAccessionNumber))
pd@data$Outliers[which(pd@data$SampleAccessionNumber %in% outliers)] <- "Yes"
autoplot(pca.br, data = pData(pd), colour="Outliers")
write.table(pd@data, paste("../pdata/pd_", studies[i,]$ID, ".tsv", sep=""), sep="\t", quote=FALSE)

ReadAffy2 <- function (..., filenames = character(0), widget = getOption("BioC")$affy$use.widgets,
                       compress = getOption("BioC")$affy$compress.cel, celfile.path = NULL,
                       sampleNames = NULL, phenoData = NULL, description = NULL,
                       notes = "", rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE,
                       verbose = FALSE, sd = FALSE, cdfname = NULL)
{
    l <- AllButCelsForReadAffy(..., filenames = filenames, widget = widget,
                               celfile.path = celfile.path, sampleNames = sampleNames,
                               phenoData = phenoData, description = description)

    ret <- read.affybatch(filenames = l$filenames, phenoData = l$phenoData,
                          description = l$description, notes = notes, compress = compress,
                          rm.mask = rm.mask, rm.outliers = rm.outliers, rm.extra = rm.extra,
                          verbose = verbose, sd = sd, cdfname = cdfname)
    sampleNames(ret) <- l$sampleNames
    return(ret)
}
