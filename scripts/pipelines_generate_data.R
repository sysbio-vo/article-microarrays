library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv4.db)
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(WGCNA)
source("probes_utils.R")

studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
common_genes <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)
pipe_types <- c("max", "mean", "random", "maxoverall", "scores", "brainarray")
pipe_type <- pipe_types[5]

# Reading rnaseq data
rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)

for (i in 1:length(studies$ID)) {
  if (grepl("Affy", studies[i,]$platform)) {
    platform <- "affymetrix"
  } else {
    platform <- "illumina"
  }

  # Reading phenodata
  pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")

  if ((platform=="affymetrix")&&(pipe_type=="brainarray")) {
    # Reading expression data
    exprs <- read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), header=TRUE)

    probesetsID <- rownames(exprs)
    probesetsID_EntrezID <- removeMultimappedProbes(probesetsID, studies[i,]$platformAbbr, brainarray=TRUE)
    exprs <- subsetNrenameRows(exprs, probesetsID_EntrezID)
    exprs <- exprs[rownames(exprs) %in% common_genes$ENTREZID,]

    write.table(exprs, paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_brainarray.tsv", sep=""), sep="\t", quote=FALSE)
    
    arseq <- mergeWithRNAseq(exprs, rnaseq, pdata, common_genes) 
    write.table(arseq, paste("../arseq/", studies[i,]$ID, "_arseq_brainarray.tsv", sep=""), sep="\t", quote=FALSE)
  }
  if (pipe_type!="brainarray") {
  
    # Reading expression data
    if (pipe_type!="scores") {
      exprs<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)
    
      probesetsID <- rownames(exprs)
      probesetsID_EntrezID <- removeMultimappedProbes(probesetsID, studies[i,]$platformAbbr, brainarray=FALSE)
      exprs <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]
    }
    if (pipe_type=="maxoverall") {
      collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
      exprs <- collapsed$datETcollapsed
      exprs <- as.data.frame(exprs)
    }
    if (pipe_type=="max") {
      pdata<-pdata[which(pdata$CancerType=="TNBC"),]
      pdata<-pdata[which(is.na(pdata$Outliers)),]
      exprs <- exprs[,colnames(exprs) %in% pdata$SampleAccessionNumber]

      collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
      exprs <- collapsed$datETcollapsed
      exprs <- as.data.frame(exprs)
    }
    if (pipe_type=="mean") {
      collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="Average")  
      exprs <- collapsed$datETcollapsed
      exprs <- as.data.frame(exprs)
    }
    if (pipe_type=="random") {
      # Duplicate genes
      n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
      dupls <- n_occur[n_occur$Freq > 1,]$Var1
      uniques <- n_occur[n_occur$Freq == 1,]$Var1
      ind <- c()
      # Choose randomly between probesets mapped to the same gene
      for (dupl in dupls) {
        ind <- c(ind, sample(which(probesetsID_EntrezID$ENTREZID==dupl),1))
      }
      probesetsID_EntrezID_random <- probesetsID_EntrezID[ind,]
      probesetsID_EntrezID_random <- rbind(probesetsID_EntrezID_random,
                                           probesetsID_EntrezID[probesetsID_EntrezID$ENTREZID %in% uniques,])
      
      exprs <- subsetNrenameRows(exprs, probesetsID_EntrezID_random)
    }
    if (pipe_type=="scores") {
      exprs <- read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)
      probesetsID_EntrezID <- read.table(paste("../scores/", studies[i,]$platformAbbr, "_scores.tsv", sep=""), header=TRUE)
      probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
      # Remove probes mapped to different genes
      n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
      uniques <- n_occur[n_occur$Freq == 1,]$Var1
      probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
      # Select genes with best BLAST score
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
      
      exprs <- subsetNrenameRows(exprs, probesetsID_EntrezID_scores)
    }
    
    if (pipe_type!="max") {
      exprs <- exprs[rownames(exprs) %in% common_genes$ENTREZID,]
      write.table(exprs, paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
    }
    
    arseq <- mergeWithRNAseq(exprs, rnaseq, pdata, common_genes) 
    write.table(arseq, paste("../arseq/", studies[i,]$ID, "_arseq_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
  }
}  