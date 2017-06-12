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
data(aips.models)
## Get an example gene signature esr1 (estrogen)
#esr1.gs = aips.models[[1706]]$gs.bresat
ind = 0
for (i in 1:1733) {
  if (aips.models[[i]]$gs.bresat$NAME == "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL") {
    ind = i
  }
}
esr1.gs = aips.models[[ind]]$gs.bresat
esr1.gs$DESCRIPTION

esr1.list <- list()
#pipe_types <- c("brainarray", "maxoverall", "mean", "scores", "random")
pipe_types <- c("maxoverall", "mean", "scores", "random")
i=4
for (pipe_type in pipe_types) {
  pdata <- read.table(paste("../pdata/combined_pdata.tsv"), header=TRUE, sep="\t", stringsAsFactors = FALSE)
  eset<-read.table(paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_", pipe_type, ".tsv", sep=""), header=TRUE)
  pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]
  pdata <- pdata[which(is.na(pdata$Outliers)),]
  pdata <- pdata[which(pdata$CancerType != "Healthy"),]
  #pdata <- pdata[which(pdata$CancerType != "Her2"),]
  eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]
  brain.ROI95.esr1.gs <- ROIq(eset, rownames(eset),
                              list(up=esr1.gs$ENTREZ$up,down=esr1.gs$ENTREZ$down))
  esr1.list <- c(esr1.list, list(brain.ROI95.esr1.gs$cl))
  print(table(brain.ROI95.esr1.gs$cl)/ncol(eset))
} 

for (i in 1:length(pipe_types)) {
  t <- table(esr1.list[[i]], pdata$CancerType)
  t[,2] <- t[,2] + t[,3]
  t <- t[,-3]
  colnames(t)[2] <- "Luminal"
  fph <- t[1, 1] + t[1, 3] + t[2, 2] + t[3, 2]
  tnh <- t[2, 3]
  rh <- fph / (fph + tnh); 
  fpl <- t[2, 1] + t[2, 2] + t[1, 3] + t[3, 3]
  tnl <- t[1, 2]
  rl <- fpl / (fpl + tnl); 
  #print(c(rh, rl))
  print(rh+rl)
  #print(t)
}


table(esr1.list[[1]], pdata$Grade)

no_illu = TRUE
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
  brain.ROI95.esr1.gs <- ROIq(eset, rownames(eset),
                              list(up=esr1.gs$ENTREZ$up,down=esr1.gs$ENTREZ$down))
  esr1.list <- c(esr1.list, list(brain.ROI95.esr1.gs$cl))
  print(table(brain.ROI95.esr1.gs$cl)/ncol(eset))
}

### GENEFU
library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(caret)

i = 3
#pipe_types <- c("brainarray", "maxoverall", "mean", "scores", "random")
pipe_types <- c("maxoverall", "random")
pipe_type = pipe_types[1]
for (pipe_type in pipe_types) {
pdata <- read.table(paste("../pdata/combined_pdata.tsv"), header=TRUE, sep="\t", stringsAsFactors = FALSE)
eset<-read.table(paste("../allsamples_exprs/", studies[i,]$ID, "_allsamples_exprs_", pipe_type, ".tsv", sep=""), header=TRUE)
pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]
pdata <- pdata[which(is.na(pdata$Outliers)),]
pdata <- pdata[which(pdata$CancerType != "Healthy"),]
#pdata <- pdata[which(pdata$CancerType != "Her2"),]
eset <- eset[,colnames(eset) %in% pdata$SampleAccessionNumber]

#EntrezID_Symbol<-select(org.Hs.eg.db, rownames(eset), c("SYMBOL", "GENENAME"))
#annot <- data.frame(EntrezGene.ID = rownames(eset), Gene.Symbol = EntrezID_Symbol$SYMBOL)

annot <- data.frame(EntrezGene.ID = rownames(eset))
rownames(annot) <- annot$probe <- annot$EntrezGene.ID

pam50 <- molecular.subtyping(sbt.model = "pam50", data = t(eset),
                             annot = annot, do.mapping = TRUE)
#table(pam50$subtype)

scmod2 <- molecular.subtyping(sbt.model = "AIMS", data = t(eset),
                              annot = annot,do.mapping = TRUE)
#table(scmod2$subtype)

Basals<-names(which(scmod2$subtype == "ER-/HER2-"))
#Select samples pertaining to HER2 Subtype
HER2s<-names(which(scmod2$subtype == "HER2+"))
#Select samples pertaining to Luminal Subtypes
LuminalB<-names(which(scmod2$subtype == "ER+/HER2- High Prolif"))
LuminalA<-names(which(scmod2$subtype == "ER+/HER2- Low Prolif"))
scmod2$subtype[LuminalB]<-"LumB"
scmod2$subtype[LuminalA]<-"LumA"
scmod2$subtype[Basals]<-"Basal"
scmod2$subtype[HER2s]<-"Her2"

#table(pam50$subtype)
#table(scmod2$subtype)

pdata[pdata$CancerType=="lumA",]$CancerType <- "LumA"
pdata[pdata$CancerType=="Luminal B",]$CancerType <- "LumB"
pdata[pdata$CancerType=="TNBC",]$CancerType <- "Basal"
if (length(which(pam50$subtype == "Normal")) > 0) {
  ind = which(pam50$subtype == "Normal")
  pam50$subtype <- pam50$subtype[-ind]
  scmod2$subtype <- scmod2$subtype[-ind]
  pdata <- pdata[-ind,]
  eset <- eset[,-ind]
}
pam50$subtype <- as.character(pam50$subtype)
confusionMatrix(pam50$subtype, scmod2$subtype)
confusionMatrix(pam50$subtype, pdata$CancerType)

### pamr
library(pamr)
#pamr.data <- list(x = as.matrix(eset), y = t(pdata$CancerType), geneid=t(rownames(eset)))
pamr.data <- list(x = as.matrix(eset), y = t(pam50$subtype), geneid=t(rownames(eset)))
train <- pamr.train(pamr.data)
results <- pamr.cv(train, pamr.data)
pamr.confusion(results, threshold=4.0)
#fdr.obj <- pamr.fdr(train, pamr.data)
#pamr.plotfdr(fdr.obj)

}