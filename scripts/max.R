detach("package:MASS", unload=TRUE)

pipe_type <- "random"

### Initial info
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("../expression_data/common_genes.txt", header = FALSE, sep = "\t")
genes_common_long <- read.table("../expression_data/Common_genes_commonPipe.txt", header = FALSE, sep = "\t")

for (i in 1:length(studies$ID)) {
  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".txt", sep=""), header=TRUE, sep="\t")
  # Reading expression data
  rnaseq <- read.table("../expression_data/rnaseq_data_processed_sum.tsv", sep="\t", header=TRUE)
  eset<-read.table(paste("../expression_data/", studies[i,]$ID, "_", pipe_type, ".txt", sep=""), header=TRUE)
  rownames(eset) <- genes_common_long[,1]
  # Select TNBC samples only
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  eset <- eset[colnames(eset) %in% rownames(pdata)]
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
  write.table(exprs, paste("../expression_data/", studies[i,]$ID, "_exprs_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}
