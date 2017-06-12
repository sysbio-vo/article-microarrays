library(limma)
library(org.Hs.eg.db)
library(VennDiagram)
source("degs_utils.R")

studies <- read.table("../general/studies.tsv", header = TRUE, sep = "\t")
pdata <-read.table(paste("../pdata/combined_pdata.tsv"), header=TRUE, sep="\t")
#pdata <- pdata[is.na(pdata$Outliers),]
pipe_types <- c("brainarray", "maxoverall", "mean", "scores", "random")
pipe_type <- pipe_types[1]

degs.list <- list()
no_illu = FALSE
for (pipe_type in pipe_types) {
  eset <- read.table(paste("../allsamples_exprs_merged/allsamples_exprs_merged_", pipe_type, ".tsv", sep=""),
                     stringsAsFactors = FALSE, header=TRUE)
  pdata <-read.table(paste("../pdata/combined_pdata.tsv"), header=TRUE, sep="\t")
  pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]
  pdata <- pdata[which(is.na(pdata$Outliers)),]
  if (no_illu) {
    pdata <- pdata[which(pdata$DataSetAccesionNumber != "GSE60785"),]
  }
  eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  
  degs <- getDEGS(c("TNBC", "lumA"), pdata, eset)
  degs <- filterDEGS(degs, 0.05, 1, adj=FALSE)[, c(1:9)]
  degs.list <- c(degs.list, list(degs))
}

if (no_illu) {
  save(degs.list, file = "degs_noillu.RData")
} else {
  save(degs.list, file = "degs.RData")
}



filename = "../plots/venn_merged_br_max_sc.pdf"
if (no_illu) {
  filename = "../plots/venn_merged_br_max_sc_noillu.pdf"
}
venn.diagram(list(brainarray = degs.list[[1]]$ENTREZID, maxoverall = degs.list[[2]]$ENTREZID, scores= degs.list[[4]]$ENTREZID),
             fill = c("red", "green", "blue"),
             alpha = 0.5, cex = 0.6, cat.fontface = 1, lty = 1, fontfamily = 3, 
             height = 1600, width = 1600, margin = 0.15, cat.cex = 0.7,
             filename = filename)

pipe_type <- pipe_types[5]

if (pipe_type=="brainarray") {
  len = length(studies$ID)-1
} else {
  len = length(studies$ID)
}

degs.list <- list()
for (i in 1:len) {
  pdata <- read.table(paste("../pdata/combined_pdata.tsv"), header=TRUE, sep="\t")
  eset<-read.table(paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_", pipe_type, ".tsv", sep=""), header=TRUE)
  pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]
  pdata <- pdata[which(is.na(pdata$Outliers)),]
  eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  
  degs <- getDEGS(c("TNBC", "lumA"), pdata, eset)
  degs <- filterDEGS(degs, 0.05, 1, adj=FALSE)[, c(1:9)]
  degs.list <- c(degs.list, list(degs))
} 
save(degs.list, file = paste("degs_", pipe_type, ".RData", sep=""))

load(file = paste("degs_", pipe_types[1], ".RData", sep=""))
degs.list.br <- degs.list
load(file = paste("degs_", pipe_types[2], ".RData", sep=""))
degs.list.max <- degs.list
load(file = paste("degs_", pipe_types[3], ".RData", sep=""))
degs.list.mean <- degs.list
load(file = paste("degs_", pipe_types[4], ".RData", sep=""))
degs.list.scores <- degs.list
load(file = paste("degs_", pipe_types[5], ".RData", sep=""))
degs.list.random <- degs.list


filename = "../plots/venn_GSE60785_r_max_sc.pdf"
venn.diagram(list(random = degs.list.random[[4]]$ENTREZID, maxoverall = degs.list.max[[4]]$ENTREZID, scores = degs.list.scores[[4]]$ENTREZID),
             fill = c("red", "green", "blue"),
             alpha = 0.5, cex = 0.6, cat.fontface = 1, lty = 1, fontfamily = 3, 
             height = 1600, width = 1600, margin = 0.15, cat.cex = 0.7,
             filename = filename)

filename = "../plots/venn_GSE19615_r_max_sc.pdf"
venn.diagram(list(random = degs.list.random[[2]]$ENTREZID, maxoverall = degs.list.max[[2]]$ENTREZID, scores = degs.list.scores[[2]]$ENTREZID),
             fill = c("red", "green", "blue"),
             alpha = 0.5, cex = 0.6, cat.fontface = 1, lty = 1, fontfamily = 3, 
             height = 1600, width = 1600, margin = 0.15, cat.cex = 0.7,
             filename = filename)

filename = "../plots/venn_GSE65194_GSE19615_GSE60785_max.pdf"
venn.diagram(list(GSE65194 = degs.list.max[[1]]$ENTREZID, GSE19615 = degs.list.max[[2]]$ENTREZID, GSE60785 = degs.list.max[[4]]$ENTREZID),
             fill = c("red", "green", "blue"),
             alpha = 0.5, cex = 0.6, cat.fontface = 1, lty = 1, fontfamily = 3, 
             height = 1600, width = 1600, margin = 0.15, cat.cex = 0.7,
             filename = filename)

filename = "../plots/venn_GSE65194_GSE19615_GSE60785_rand.pdf"
venn.diagram(list(GSE65194 = degs.list.random[[1]]$ENTREZID, GSE19615 = degs.list.random[[2]]$ENTREZID, GSE60785 = degs.list.random[[4]]$ENTREZID),
             fill = c("red", "green", "blue"),
             alpha = 0.5, cex = 0.6, cat.fontface = 1, lty = 1, fontfamily = 3, 
             height = 1600, width = 1600, margin = 0.15, cat.cex = 0.7,
             filename = filename)


library(AIPS)
library(Biobase)
openVignette()
data(mcgill.gq)
mcgill.aips <- mclapply.AIPS(mcgill.gq$D,
                             mcgill.gq$EntrezID)

table(mcgill.aips$cl)/ncol(mcgill.aips$cl)
save(mcgill.aips, file = paste("mcgill.aips.RData", sep=""))
load("mcgill.aips.RData")

brain.aips <- mclapply.AIPS(eset, rownames(eset))
table(brain.aips$cl)/ncol(eset)

table(mcgill.aips$cl[,1],brain.aips$cl[,1])


data(aips.models)
## Get an example gene signature esr1 (estrogen)
esr1.gs = aips.models[[1706]]$gs.bresat
mcgill.ROI95.esr1.gs <- ROIq (mcgill.gq$D,
                              mcgill.gq$EntrezID,
                              list(up=esr1.gs$ENTREZ$up,down=esr1.gs$ENTREZ$down))
## Print a summary of all the assignments
table(mcgill.ROI95.esr1.gs$cl)


brain.ROI95.esr1.gs <- ROIq(eset, rownames(eset),
                            list(up=esr1.gs$ENTREZ$up,down=esr1.gs$ENTREZ$down))
table(brain.ROI95.esr1.gs$cl)

colnames(mcgill.gq$D)
