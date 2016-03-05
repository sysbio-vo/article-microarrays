# RNA-seq data
samples<-read.table("../RNA-seq/SamplesRNASeq.txt",header=F)
samples<-as.character(samples$V1)
Transcripts<-as.character(read.table("../RNA-seq/SRR1313132_cufflinks/genes.fpkm_tracking",header=T)$tracking_id)
vec<-Transcripts
for(samp in samples){
  file1<-paste("../RNA-seq/",samp,"_cufflinks/","genes.fpkm_tracking",sep="",collapse = "")
  data_in<-read.table(file1,header=T)
  data_in<-data_in[match(Transcripts,data_in$tracking_id),]
  FPKM<-as.numeric(as.character(data_in$FPKM))
  vec<-cbind(vec,FPKM)
  print(samp)
}
vec1<-as.data.frame(vec)
vec1<-vec1[,-1]
vec1<-vec1[which(duplicated(Transcripts)==FALSE),]
rownames(vec1)<-Transcripts[which(duplicated(Transcripts)==FALSE)] # Removes duplicated (few)
for(i in 1:ncol(vec1)) vec1[,i]<-as.numeric(as.character(vec1[,i]))

write.table(vec1, "../rnaseq/rnaseq_data.tsv", sep="\t", quote=FALSE)
rnaseqExprs <- read.table("../rnaseq/rnaseq_data.tsv", header = TRUE, sep = "\t")

rnaseqExprs <- rowMeans(rnaseqExprs)
# Takes only transcripts with at least 0.01 fpkm expression
rnaseqExprs <- rnaseqExprs[which(rnaseqExprs>=0.01)]
# Convert RefSeqID to EntrezID
library(org.Hs.eg.db)
entrez <- mget(names(rnaseqExprs), org.Hs.egREFSEQ2EG, ifnotfound = NA)
rnaseqExprs <- data.frame(entrezID=unlist(entrez), val=rnaseqExprs, row.names=names(unlist(entrez)), stringsAsFactors=FALSE)
rnaseqExprs <- rnaseqExprs[which(rnaseqExprs$entrezID!="NA"),]
# Max expression between transcripts
n_occur <- data.frame(table(rnaseqExprs$entrezID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
dupls <- n_occur[n_occur$Freq > 1,]$Var1
max_ind <- c()
for(geneID in dupls){
  # Max
  i <- which(rnaseqExprs$entrezID==geneID)
  max <- which.max(rnaseqExprs[i,]$val)
  max_ind <- c(max_ind, i[max])
  # Sum
  #i <- which(rnaseqExprs$entrezID==geneID)
  #sum <- sum(rnaseqExprs[i,]$val)
  #newrow <- rnaseqExprs[i[1],]
  #newrow$val <- sum
  #rnaseqExprs <- rnaseqExprs[-i,]
  #rnaseqExprs = rbind(rnaseqExprs,newrow)
}
 
rnaseqExprsNoDupls <- rnaseqExprs[rnaseqExprs$entrezID %in% uniques,]
rnaseqExprsNoDupls <- rbind(rnaseqExprsNoDupls, rnaseqExprs[max_ind,])
rnaseqExprs <- rnaseqExprsNoDupls

rnaseqExprs$val <- unlist(lapply(rnaseqExprs$val, log2))
write.table(rnaseqExprs, "../rnaseq/rnaseq_data_processed_max_long.tsv", sep="\t", quote=FALSE)
