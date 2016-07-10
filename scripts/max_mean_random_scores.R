# Load HuGene and hgu133 BrainArray packages
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv4.db)
library(WGCNA)

### Initial info
pipe_types <- c("max", "mean", "random", "maxoverall", "scores")
pipe_type <- pipe_types[5]
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
genes_common <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)

# Generate aggregated microarray-rnaseq expression files for plots feeding
# Didn't wrap it up in 'for (pipe_type in pipe_types)' cycle for testing purposes and
# because it takes too much time on my laptop. So execute the cycle for all the pipelines.
for (i in 1:length(studies$ID)) {
  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
  # Reading expression data
  rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  if (grepl("Affy", studies[i,]$platform)) {
    platform <- "affymetrix"
  } else {
    platform <- "illumina"
  }
  eset<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)
  if (!pipe_type=="scores") {
    # List of probesets IDs
    probesetsID<-rownames(eset)
    # List of corresponding gene IDs
    if (studies[i,]$platformAbbr=="hugene10st") {
      probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
    } else {
      probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
    }
    # Replace probesetsIDs with gene IDs in expression data matrix
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
    # Duplicate probes
    n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
    # Filter eset based on left probesets
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID$PROBEID),]
    # Select TNBC samples only
    eset.all <- eset
  }
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  eset <- eset[colnames(eset) %in% pdata$SampleAccessionNumber]
  # Collapse rows based on one of the algorithms
  if (pipe_type=="maxoverall") {
    collapsed = collapseRows(eset.all, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
    eset <-collapsed$datETcollapsed
    pdata<-pdata[which(pdata$CancerType=="TNBC"),]
    pdata<-pdata[which(is.na(pdata$Outliers)),]
    eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  }
  if (pipe_type=="max") {
    collapsed = collapseRows(eset, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
    eset <-collapsed$datETcollapsed
  }
  if (pipe_type=="mean") {
    collapsed = collapseRows(eset, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="Average")  
    eset <-collapsed$datETcollapsed
  }
  if (pipe_type=="random") {
    # Duplicate genes
    n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
    dupls <- n_occur[n_occur$Freq > 1,]$Var1
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    ind <- c()
    for (dupl in dupls) {
      ind <- c(ind, sample(which(probesetsID_EntrezID$ENTREZID==dupl),1))
    }
    probesetsID_EntrezID_random <- probesetsID_EntrezID[ind,]
    probesetsID_EntrezID_random <- rbind(probesetsID_EntrezID_random,
                                         probesetsID_EntrezID[probesetsID_EntrezID$ENTREZID %in% uniques,])
    
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID_random$PROBEID),]
    eset <- eset[order(match(rownames(eset),probesetsID_EntrezID_random$PROBEID)),]
    rownames(eset) <- probesetsID_EntrezID_random$ENTREZID
  }
  if (pipe_type=="scores") {
    probesetsID_EntrezID<-read.table(paste("../scores/", studies[i,]$platformAbbr, "_scores.tsv", sep=""), header=TRUE)
    probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
    n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
    uniques <- n_occur[n_occur$Freq == 1,]$Var1
    dupls <- n_occur[n_occur$Freq > 1,]$Var1
    ind <- c()
    for (dupl in dupls) {
      v <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID==dupl),]
      ind <- c(ind, as.numeric(rownames(v[which.max(v$Cov_Spec),])))
    }
    probesetsID_EntrezID_scores <- probesetsID_EntrezID[which(rownames(probesetsID_EntrezID) %in% ind),]
    probesetsID_EntrezID_scores <- rbind(probesetsID_EntrezID_scores,
                                         probesetsID_EntrezID[probesetsID_EntrezID$ENTREZID %in% uniques,])
    
    eset <- eset[which(rownames(eset) %in% probesetsID_EntrezID_scores$PROBEID),]
    rownames(eset) <- probesetsID_EntrezID_scores[match(rownames(eset),probesetsID_EntrezID_scores$PROBEID),]$ENTREZID
    #eset <- eset[order(match(rownames(eset),probesetsID_EntrezID_scores$PROBEID)),]
    #rownames(eset) <- probesetsID_EntrezID_scores$ENTREZID
  }
  # Select only common genes for all the platforms and pipelines
  eset <- eset[rownames(eset) %in% genes_common$ENTREZID,]
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
