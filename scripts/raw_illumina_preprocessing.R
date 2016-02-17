library(lumi)
library(limma)
library(illuminaHumanv4.db)
library(WGCNA)
library(BioBase)
library(lumidat)
library(beadarray)

# General
dir = "../Raw_Data/GSE66666/"
pdata = read.table(paste(dir, "GSE66666_pdata.tsv", sep=""), sep="\t", head=TRUE)
idatfiles = dir(paste(dir, "idat", sep=""), pattern="idat", full.name=TRUE)
bgxfile = dir(dir, pattern="bgx", full.name=TRUE)
txtbgxfile = dir(dir, pattern="txt", full.name=TRUE)

# Limma
raw.data <- read.idat(idatfiles, bgxfile)
N.data = neqc(raw.data)
probesetsID_EntrezID<-select(illuminaHumanv4.db, N.data$genes$Probe_Id, "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Subset EList
indices <- which(N.data$genes$Probe_Id %in% probesetsID_EntrezID$PROBEID)
N.data <- N.data[indices, ]
limma.data <- N.data
boxplot(limma.data$E)
datET <- as.data.frame(N.data$E, row.names=N.data$genes$Probe_Id)
datET <- datET[,which(pdata$Subs=="Yes")]

# Lumidat

raw.data <- lumiR.idat(idatfiles, probeID="ProbeID", manifestfile=txtbgxfile,
                    detectionTh=0.0001, QC=TRUE)
pData(raw.data) <- pdata
T.data   <- lumiT(raw.data, method="vst")
# The lumi package provides lumiB function for background correction. We sup-
# pose the BeadStudio output data has been background corrected. Therefore, no
# background correction used by default. A method ’bgAdjust’ is designed to ap-
# proximate what BeadStudio does for background adjustment. In the case when
# ’log2’ transform is used in the lumiT step, the background correction method
# (’forcePositive’) will be automatically used, which basically adds an offset (minus
# minimum value plus one) if there is any negative values to force all expression
# values to be positive.
N.data<-lumiN(T.data, method = "rsn")

probesetsID_EntrezID <- select(illuminaHumanv4.db, rownames(N.data), "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Subset EList
indices <- which(rownames(N.data) %in% probesetsID_EntrezID$PROBEID)
N.data <- exprs(N.data)[indices, ]
lumi.data <- N.data
boxplot(lumi.data)
datET <- N.data
datET <- datET[,which(pdata$Subs=="Yes")]

# BeadArray
raw.data <- readIdatFiles(idatfiles)
N.data <- normaliseIllumina(channel(raw.data, "Green"),
                                              method="vsn", transform="none")

probesetsID_EntrezID <- select(illuminaHumanv4.db, rownames(N.data), "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Duplicate probes
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Subset EList
indices <- which(rownames(N.data) %in% probesetsID_EntrezID$PROBEID)
N.data <- exprs(N.data)[indices, ]
bead.data <- N.data
boxplot(bead.data)

datET <- N.data
datET <- datET[,which(pdata$Subs=="Yes")]


# Collapse rows
collapsed = collapseRows(datET, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")

# Intersect eset and rnaseq
eset <-collapsed$datETcollapsed
rnaseq <- read.table("../rnaseq_data_processed_sum.tsv", sep="\t", header=TRUE)
eset <- eset[rownames(eset) %in% rnaseq$entrezID,]
eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
rownames(rnaseq) <- rnaseq$entrezID
rnaseq <- rnaseq[which(rownames(rnaseq) %in% rownames(eset)),]
rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(eset))),]

# Assemble joint dataset
exprs <- data.frame(row.names=rownames(eset), eset=eset$val, rnaseq=rnaseq$val)
write.table(exprs, paste("../GSE66666_exprs_max.tsv", sep=""), sep="\t", quote=FALSE)
