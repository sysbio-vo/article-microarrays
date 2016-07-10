# Loading libraries
library(cowplot)
library(plyr)
library(ggplot2)
library(ggfortify)

pipe_type = "maxoverall"

pdata <-read.table("../pdata/combined_pdata.tsv", header=TRUE, sep="\t")
eset<-read.table(paste("../combined/combined_expression_", pipe_type,".tsv", sep=""), header=TRUE)
pdata$Grade <- as.character(pdata$Grade)
pdata <- pdata[pdata$SampleAccessionNumber %in% colnames(eset),]

pca = prcomp(t(eset))

title <- ggdraw() + draw_label("PCA plot of combined microarray dataset. Max pipeline", fontface='bold')  
pl1 <- autoplot(pca, data = pdata, colour="CancerType")
pl2 <- autoplot(pca, data = pdata, colour="Grade")
pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

pdata.grade<-pdata[which(pdata$Grade!="NA"),]
eset.grade <- eset[,colnames(eset) %in% pdata.grade$SampleAccessionNumber]
pca.grade = prcomp(t(eset.grade))

title <- ggdraw() + draw_label("Samples with grade present", fontface='bold')  
pl1 <- autoplot(pca.grade, data = pdata.grade, colour="CancerType")
pl2 <- autoplot(pca.grade, data = pdata.grade, colour="Grade")
pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))

pl <- plot_grid(pl, pl3, nrow=2)

# Save plot for manual quality control
save_plot(paste("../plots/combined_PCA.pdf", sep=""),
          pl, base_height=5,
          base_aspect_ratio=2.3,
          nrow=1)

# In case one is interested in per grade representation
pdata.grade3 <-pdata[which(pdata$Grade=="3"),]
eset.grade3 <- eset[,colnames(eset) %in% pdata.grade3$SampleAccessionNumber]
pca.grade3 = prcomp(t(eset.grade3))
pdata.grade2 <-pdata[which(pdata$Grade=="2"),]
eset.grade2 <- eset[,colnames(eset) %in% pdata.grade2$SampleAccessionNumber]
pca.grade2 = prcomp(t(eset.grade2))
pdata.grade1 <-pdata[which(pdata$Grade=="1"),]
eset.grade1 <- eset[,colnames(eset) %in% pdata.grade1$SampleAccessionNumber]
pca.grade1 = prcomp(t(eset.grade1))


title <- ggdraw() + draw_label("Samples with grade=1", fontface='bold')  
pl1 <- autoplot(pca.grade1, data = pdata.grade1, colour="CancerType")
pl1 <- plot_grid(title, pl1, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label("Samples with grade=2", fontface='bold')  
pl2 <- autoplot(pca.grade2, data = pdata.grade2, colour="CancerType")
pl2 <- plot_grid(title, pl2, ncol=1, rel_heights=c(0.1, 1))
title <- ggdraw() + draw_label("Samples with grade=3", fontface='bold')  
pl3 <- autoplot(pca.grade3, data = pdata.grade3, colour="CancerType")
pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))

pl4 <- plot_grid(pl1, pl2, pl3, ncol=3)

# Save plot for manual quality control
save_plot(paste("../plots/combined_PCA_pergrade_br.pdf", sep=""),
          pl4, base_height=5,
          ncol=3)