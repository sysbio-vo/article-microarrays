# Remove probes mapped to different genes
removeMultimappedProbes <- function(probesetsID, platformAbbr, brainarray=FALSE) {
  # Get EntrezID
  if ((platformAbbr=="hugene10st")&&(brainarray==FALSE)) {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, "transcriptcluster.db", sep="")), probesetsID, "ENTREZID")
  } else {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")
  }
  if (brainarray) {
    probesetsID_EntrezID<-select(get(paste(platformAbbr, "hsentrezg.db", sep="")), probesetsID, "ENTREZID")
  }

  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  # Remove probes mapped to different genes
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  return(probesetsID_EntrezID)
}

