# Load HuGene and hgu133 BrainArray packages
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
# Load Bioconductor annotation packages
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv4.db)
#library(illuminaHumanv3.db)

# Load studies description
studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
platforms <- levels(studies$platformAbbr)
ind = c()
for (pl in platforms) {
  ind = c(ind, which.max(studies$platformAbbr==pl))
}

IDs <- c()
platform = ""
for (i in ind) {
  # Affymetrix/Illumina probesets
  if (grepl("Illu", studies[i,]$platform)) {
    platform = "illumina" 
  } else {
    platform = "affymetrix"
  }
  
  eset<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_", platform, ".tsv", sep=""), header=TRUE)
  probesetsID <- rownames(eset)
  # Get probesets ID to EntrezID mapping
  if (studies[i,]$platformAbbr=="hugene10st") {
    probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
  } else {
    probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
  }
  
  # Delete NAs
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  # Delete probesets, which have multiple mappings to genes
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  #doubles <- n_occur[n_occur$Freq > 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  IDs <- c(IDs, list(probesetsID_EntrezID))  
  # Generate list of genes, represented by multiple probesets
  n_occur <- data.frame(table(probesetsID_EntrezID$ENTREZID))
  dupls <- n_occur[n_occur$Freq > 1,]$Var1
  write.table(dupls, paste("../general/onetomany_genes_", studies[i,]$platformAbbr,".txt", sep=""), sep="\t", quote=FALSE,
              col.names=c("ENTREZID"), row.names=FALSE)
  
  # Brainarray probesets
  eset.br<-read.table(paste("../preprocessed/", studies[i,]$ID, "_preprocessed_brainarray.tsv", sep=""), header=TRUE)
  probesetsID <- rownames(eset.br)
  probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  IDs <- c(IDs, list(probesetsID_EntrezID))  
}

# Intersect among all the platforms
common <- IDs[[1]]$ENTREZID
for (i in 2:length(IDs)) {
  common <- intersect(common, IDs[[i]]$ENTREZID)
}

# Intersect with scored probesets
for (i in ind) {
  probesetsID_EntrezID<-read.table(paste("../scores/", studies[i,]$platformAbbr, "_scores.tsv", sep=""), header=TRUE)
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  common <- intersect(common, probesetsID_EntrezID$ENTREZID)
}

# Intersect with RNA-seq
rnaseq <- read.table("../rnaseq/rnaseq_data_processed_sum_long.tsv", sep="\t", header=TRUE)
common <- intersect(common, rnaseq$entrezID)
write.table(common, "../general/common_genes_list.txt", sep="\t", quote=FALSE,
            col.names=c("ENTREZID"), row.names=FALSE)
