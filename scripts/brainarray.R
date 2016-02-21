# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
detach("package:MASS", unload=TRUE)

### Initial info
pipe_type <- "brainarray"
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)
ind <- which(grepl("Affy", studies$platform))

# Generate aggregated microarray-rnaseq expression files for plots feeding
for (i in ind) {
  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
  # Reading expression data
  rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  eset<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", pipe_type, ".tsv", sep=""), header=TRUE)
  # List of probesets IDs
  probesetsID<-rownames(eset)
  # List of corresponding gene IDs
  probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
  # Replace probesetsIDs with gene IDs in expression data matrix
  eset <- eset[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  rownames(eset) <- probesetsID_EntrezID$ENTREZID[which(probesetsID_EntrezID$ENTREZID!="NA")]
  # Select only common genes for all the platforms and pipelines
  eset <- eset[rownames(eset) %in% genes_common$ENTREZID,]
  # Select TNBC samples only
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  eset <- eset[colnames(eset) %in% pdata$SampleAccessionNumber]
  # Take mean across samples
  eset <- data.frame(val=rowMeans(eset), row.names=rownames(eset))
  # Select common genes with rna-seq
  rnaseq <- rnaseq[rnaseq$entrezID %in% rownames(eset),]
  rownames(rnaseq) <- rnaseq$entrezID
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(eset))),]
  rnaseq <- rnaseq[-1]
  # Assemble joint dataset
  exprs <- data.frame(row.names=rownames(eset), eset=eset$val, rnaseq=rnaseq$val)
  write.table(exprs, paste("../exprs/", studies[i,]$ID, "_exprs_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}
