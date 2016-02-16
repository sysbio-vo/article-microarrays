# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
detach("package:MASS", unload=TRUE)

pipe_type <- "brainarray"

### Initial info
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("../expression_data/common_genes.txt", header = FALSE, sep = "\t")

ind <- which(grepl("Affy", studies$platform))
for (i in ind) {
  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[ind[i],]$ID, ".txt", sep=""), header=TRUE, sep="\t")
  # Reading expression data
  rnaseq <- read.table("../expression_data/rnaseq_data_processed_sum.tsv", sep="\t", header=TRUE)
  eset<-read.table(paste("../expression_data/", studies[ind[i],]$ID, "_", pipe_type, ".txt", sep=""), header=TRUE)
  # List of probesets IDs
  probesetsID<-rownames(eset)
  # List of corresponding gene IDs
  probesetsID_EntrezID<-select(get(paste(studies[ind[i],]$platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
  # Replace probesetsIDs with gene IDs in expression data matrix
  eset <- eset[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  rownames(eset) <- probesetsID_EntrezID$ENTREZID[which(probesetsID_EntrezID$ENTREZID!="NA")]
  # Select only common genes for all the platforms and pipelines
  eset <- eset[rownames(eset) %in% genes_common[,1],]
  # Select TNBC samples only
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  eset <- eset[colnames(eset) %in% pdata$SampleAccessionNumber]
  # Select common genes with rna-seq
  eset <- eset[rownames(eset) %in% rnaseq$entrezID,]
  eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
  # Organize rna-seq dataset
  rownames(rnaseq) <- rnaseq$entrezID
  rnaseq <- rnaseq[-1]
  rnaseq$val <- rnaseq[order(match(rownames(rnaseq),rownames(eset))),]
  rownames(rnaseq) <- rownames(eset)
  # Assemble joint dataset
  exprs <- data.frame(row.names=rownames(eset), eset=eset$val, rnaseq=rnaseq$val)
  write.table(exprs, paste("../exprs/", studies[ind[i],]$ID, "_exprs_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}

