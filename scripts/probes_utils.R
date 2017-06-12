# Remove probes mapped to different genes
removeMultimappedProbes <- function(probesetsID, platformAbbr, brainarray=FALSE) {
  # Get EntrezID
  if ((platformAbbr=="hugene10st")&&(brainarray==FALSE)) {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
  } else if (brainarray==FALSE) {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
  }
  if (brainarray) {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
  }
  
  # Remove NA entrezID
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  # Remove probes mapped to different genes
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  return(probesetsID_EntrezID)
}

subsetNrenameRows <- function(exprs, probesetsID_EntrezID) {
  exprs <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]
  rownames(exprs) <- probesetsID_EntrezID[match(rownames(exprs), probesetsID_EntrezID$PROBEID),]$ENTREZID
  return(exprs)
}

mergeWithRNAseq <- function(exprs, rnaseq, pdata, common_genes) {
  # Select only common genes for all the platforms and pipelines
  exprs <- exprs[rownames(exprs) %in% common_genes$ENTREZID,]

  # Select TNBC samples only
  pdata<-pdata[which(pdata$CancerType=="TNBC"),]
  pdata<-pdata[which(is.na(pdata$Outliers)),]
  exprs <- exprs[colnames(exprs) %in% pdata$SampleAccessionNumber]

  # Take mean across samples
  exprs <- data.frame(val=rowMeans(exprs), row.names=rownames(exprs))

  # Select common genes with rna-seq
  rnaseq <- rnaseq[rownames(rnaseq) %in% rownames(exprs), , drop=FALSE]
  rnaseq <- rnaseq[order(match(rownames(rnaseq),rownames(exprs))), , drop=FALSE]
  # Assemble joint dataset
  arseq <- data.frame(row.names=rownames(exprs), exprs=exprs$val, rnaseq=rnaseq$val)
  return(arseq)
}