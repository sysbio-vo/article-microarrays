library(org.Hs.eg.db)
library(WGCNA)

## GSE58135

samples <- read.table("../rnaseq/GSE58135/tophat/SamplesRNASeq.txt",header=F)
samples <- as.character(samples$V1)
transcripts <- as.character(read.table("../rnaseq/GSE58135/tophat/SRR1313132_cufflinks/genes.fpkm_tracking",header=T)$tracking_id)
rnaseq <- transcripts
for (samp in samples) {
  file <- paste("../rnaseq/GSE58135/tophat/",samp,"_cufflinks/","genes.fpkm_tracking",sep="",collapse = "")
  data <- read.table(file, header=T)
  data <- data[match(transcripts, data$tracking_id), ]
  FPKM <- as.numeric(as.character(data$FPKM))
  rnaseq <- cbind(rnaseq, FPKM)
  print(samp)
}

rnaseq <- as.data.frame(rnaseq)
colnames(rnaseq)[2:ncol(rnaseq)] <- samples
colnames(rnaseq)[1] <- "REFSEQID"

n_occur <- data.frame(table(rnaseq$REFSEQID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
rnaseq <- rnaseq[which(rnaseq$REFSEQID %in% uniques),]
rownames(rnaseq) <- rnaseq$REFSEQID
rnaseq <- rnaseq[-1]

write.table(rnaseq, "../rnaseq/GSE58135/GSE58135_rnaseq_tophat.tsv", sep="\t", quote=FALSE)
rnaseq <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_tophat.tsv", sep="\t")

#rnaseq <- rnaseq[!apply(rnaseq, 1, function(row) all(row == 0)), ]
RefseqID_EntrezID <- select(org.Hs.eg.db, rownames(rnaseq), c("ENTREZID"), keytype = "REFSEQ")
RefseqID_EntrezID <- RefseqID_EntrezID[which(RefseqID_EntrezID$ENTREZID!="NA"),]
rnaseq <- rnaseq[which(rownames(rnaseq) %in% RefseqID_EntrezID$REFSEQ),]

collapsed.max = collapseRows(rnaseq, RefseqID_EntrezID$ENTREZID, RefseqID_EntrezID$REFSEQ, method="MaxMean")  
collapsed.sum = collapseRows(rnaseq, RefseqID_EntrezID$ENTREZID, RefseqID_EntrezID$REFSEQ,
                             method="function", methodFunction=colSums)

rnaseq <- collapsed.max$datETcollapsed
rnaseq <- as.data.frame(rnaseq)
write.table(rnaseq, "../rnaseq/GSE58135/GSE58135_rnaseq_processed_max_long.tsv", sep="\t", quote=FALSE)

rnaseq <- collapsed.sum$datETcollapsed
rnaseq <- as.data.frame(rnaseq)
write.table(rnaseq, "../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_long.tsv", sep="\t", quote=FALSE)

# GSE60788

studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
pdata <- read.table(paste("../pdata/pdata_", studies[4,]$ID, ".tsv", sep=""), header=TRUE, sep="\t", stringsAsFactors = FALSE)
IDREF_SYMBOL <- read.table("../rnaseq/GSE60788/geo/GSE60788_map_transcriptID_geneSymbol.csv", header=TRUE, sep="\t")
rnaseq <- read.table("../rnaseq/GSE60788/geo/GSE60788_rnaseq_gex_raw.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

rnaseq <- rnaseq[, -which(grepl("replicate", colnames(rnaseq)))]
colnames(rnaseq) <- c("ID_REF", pdata$SampleAccessionNumber)

rownames(rnaseq) <- rnaseq$ID_REF
rnaseq <- rnaseq[-1]

collapsed.max = collapseRows(rnaseq, IDREF_SYMBOL$geneSymbol, IDREF_SYMBOL$transcriptID, method="MaxMean")  
collapsed.sum = collapseRows(rnaseq, IDREF_SYMBOL$geneSymbol, IDREF_SYMBOL$transcriptID,
                             method="function", methodFunction=colSums)

n_occur <- data.frame(table(IDREF_SYMBOL$geneSymbol))
SYMBOLs_EntrezID <- select(org.Hs.eg.db, as.character(n_occur$Var1), c("ENTREZID"), keytype = "SYMBOL")
SYMBOLs_EntrezID <- SYMBOLs_EntrezID[which(SYMBOLs_EntrezID$ENTREZID != "NA"), ]
n_occur <- data.frame(table(SYMBOLs_EntrezID$SYMBOL))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
SYMBOLs_EntrezID <- SYMBOLs_EntrezID[which(SYMBOLs_EntrezID$SYMBOL %in% uniques), ]

rnaseq <- collapsed.max$datETcollapsed
rnaseq <- as.data.frame(rnaseq)
rnaseq <- rnaseq[which(rownames(rnaseq) %in% SYMBOLs_EntrezID$SYMBOL), ]
rnaseq <- rnaseq[order(match(rownames(rnaseq), SYMBOLs_EntrezID$SYMBOL)),]
rownames(rnaseq) <- SYMBOLs_EntrezID$ENTREZID
write.table(rnaseq, "../rnaseq/GSE60788/GSE60788_rnaseq_processed_max_long.tsv", sep="\t", quote=FALSE)

rnaseq <- collapsed.sum$datETcollapsed
rnaseq <- as.data.frame(rnaseq)
rnaseq <- rnaseq[which(rownames(rnaseq) %in% SYMBOLs_EntrezID$SYMBOL), ]
rnaseq <- rnaseq[order(match(rownames(rnaseq), SYMBOLs_EntrezID$SYMBOL)),]
rownames(rnaseq) <- SYMBOLs_EntrezID$ENTREZID
write.table(rnaseq, "../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_long.tsv", sep="\t", quote=FALSE)



### Generate short files

GSE58135.t <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_long.tsv", sep="\t", header=TRUE)
GSE58135.s <- read.table("../rnaseq/GSE58135/GSE58135_rnaseq_subread.tsv", sep="\t")
GSE60788 <- read.table("../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_long.tsv", sep="\t")

rnaseq <- GSE58135.t
rnaseq <- data.frame(val=rowMeans(rnaseq), row.names=rownames(rnaseq))
rnaseq <- rnaseq[which(rnaseq>0.01), , drop = FALSE]
rnaseq$val <- log(rnaseq$val, 2)
write.table(rnaseq, "../rnaseq/GSE58135/GSE58135_rnaseq_processed_sum_short.tsv", sep="\t", quote=FALSE)

rnaseq <- GSE58135.s
rnaseq <- data.frame(val=rowMeans(rnaseq), row.names=rownames(rnaseq))
rnaseq <- rnaseq[which(rnaseq>-6.643856), , drop = FALSE]
write.table(rnaseq, "../rnaseq/GSE58135/GSE58135_rnaseq_subread_short.tsv", sep="\t", quote=FALSE)

rnaseq <- GSE60788
rnaseq <- data.frame(val=rowMeans(rnaseq), row.names=rownames(rnaseq))
rnaseq <- rnaseq[which(rnaseq>0.01), , drop = FALSE]
rnaseq$val <- log(rnaseq$val, 2)
write.table(rnaseq, "../rnaseq/GSE60788/GSE60788_rnaseq_processed_sum_short.tsv", sep="\t", quote=FALSE)
