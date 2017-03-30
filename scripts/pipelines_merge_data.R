library(sva)
source("probes_utils.R")
### Initial info
pipe_types <- c("brainarray", "mean", "random", "maxoverall", "scores")
pipe_type <- pipe_types[2]
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
common_genes <- read.table("../general/common_genes_list.txt", sep = "\t", header = TRUE)

for (pipe_type in pipe_types) {
  eset.all <- data.frame()
  pdata.all <- data.frame()
  if (pipe_type=="brainarray") {
    len = length(studies$ID)-1
  } else {
    len = length(studies$ID)
  }
  for (i in 1:len) {
    pdata <-read.table(paste("../pdata/pdata_", studies[i,]$ID, ".tsv", sep=""), header=TRUE, sep="\t")
    eset<-read.table(paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_", pipe_type, ".tsv", sep=""), header=TRUE)
    if (nrow(eset.all)==0) {
      eset.all <- eset
      pdata.all <- pdata
    } else {
      # This subsetting is caused by scores pipeline
      if (nrow(eset.all)>nrow(eset)) {
        eset.all <- eset.all[which(rownames(eset.all) %in% rownames(eset)),]
      } else {
        eset <- eset[which(rownames(eset) %in% rownames(eset.all)),]
      }
      eset.all <- eset.all[order(match(rownames(eset.all),rownames(eset))),]
      eset.all <- cbind(eset.all, eset)
      pdata.all <- rbind(pdata.all, pdata)
    }
  }
  batch = pdata.all$DataSetAccesionNumber
  mod = model.matrix(~as.factor(CancerType), data=pdata.all)
  eset.all = ComBat(dat=eset.all, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  write.table(eset.all, paste("../allsamples_exprs_merged/allsamples_exprs_merged_", pipe_type,".tsv", sep=""), sep="\t", quote=FALSE)
  write.table(pdata.all, paste("../pdata/combined_pdata.tsv", sep=""), sep="\t", quote=FALSE)
}

for (pipe_type in pipe_types) {
  exprs<-read.table(paste("../allsamples_exprs_merged/allsamples_exprs_merged_", pipe_type,".tsv", sep=""), header=TRUE)
  rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
  
  arseq <- mergeWithRNAseq(exprs, rnaseq, pdata.all, common_genes)
  write.table(arseq, paste("../arseq_merged/arseq_merged_", pipe_type, ".tsv", sep=""), sep="\t", quote=FALSE)
}
