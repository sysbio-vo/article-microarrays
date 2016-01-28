# RNA-seq data
samples<-read.table("RNA-Seq/SamplesRNASeq.txt",header=F)
samples<-as.character(samples$V1)
Transcripts<-as.character(read.table("RNA-Seq/SRR1313132_cufflinks/genes.fpkm_tracking",header=T)$tracking_id)
vec<-Transcripts
for(samp in samples){
  file1<-paste("RNA-Seq/",samp,"_cufflinks/","genes.fpkm_tracking",sep="",collapse = "")
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

write.table(vec1, "expression_data/rnaseq_data.tsv", sep="\t", quote=FALSE)
rnaseqExprs <- read.table("expression_data/rnaseq_data.tsv", header = TRUE, sep = "\t")

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
geneIDs <- unique(rnaseqExprs$entrezID[which(duplicated(rnaseqExprs$entrezID)==TRUE)])
for(geneID in geneIDs){
  # Mean
  #i <- which(rnaseqExprs$entrezID==geneIDs)
  #min <- which.min(rnaseqExprs[i,]$val)
  #rnaseqExprs <- rnaseqExprs[-i[min],]
  # Sum
  i <- which(rnaseqExprs$entrezID==geneID)
  sum <- sum(rnaseqExprs[i,]$val)
  newrow <- rnaseqExprs[i[1],]
  newrow$val <- sum
  rnaseqExprs <- rnaseqExprs[-i,]
  rnaseqExprs = rbind(rnaseqExprs,newrow)
}
 
rnaseqExprs1 <- rnaseqExprs
rnaseqExprs1[which(rnaseqExprs1$val==max(rnaseqExprs1$val)),]

rnaseqExprs$val <- unlist(lapply(rnaseqExprs$val, log2))

genes_common <- read.table("expression_data/common_genes.txt", header = FALSE, sep = "\t")
rnaseqExprs <- rnaseqExprs[rnaseqExprs$entrezID %in% genes_common[,1],]
#write.table(rnaseqExprs, "expression_data/rnaseq_data_processed.tsv", sep="\t", quote=FALSE)
write.table(rnaseqExprs, "expression_data/rnaseq_data_processed_sum.tsv", sep="\t", quote=FALSE)
