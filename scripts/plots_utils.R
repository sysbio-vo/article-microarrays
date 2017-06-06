# PCA to compare datasets after RMA, check batches etc
qcPCA <- function(af, br, pd, study.id) {
  pca = prcomp(t(af))
  pca.br = prcomp(t(br))
  
  title <- ggdraw() + draw_label(paste("Affymetrix Probesets.", study.id), fontface='bold')  
  pl1 <- autoplot(pca, data = pd, colour="CancerType")
  pl2 <- autoplot(pca, data = pd, colour="ScanDate")
  pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
  pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  title <- ggdraw() + draw_label(paste("Brainarray Probesets.", study.id), fontface='bold')  
  pl1 <- autoplot(pca.br, data = pd, colour="CancerType")
  pl2 <- autoplot(pca.br, data = pd, colour="ScanDate")
  pl3 <- plot_grid(pl1, pl2, ncol=2, align="hv")
  pl3 <- plot_grid(title, pl3, ncol=1, rel_heights=c(0.1, 1))
  pl <- plot_grid(pl, pl3, nrow=2, align="hv")
  
  return(pl)
}

# PCA to compare datasets after RMA, check batches etc
qcPCAilluraw <- function(exprs, pd, study.id) {
  pca = prcomp(t(exprs))

  title <- ggdraw() + draw_label(study.id, fontface='bold')  
  pl1 <- autoplot(pca, data = pdata, colour="CancerType")
  pl2 <- autoplot(pca, data = pdata, colour="Sentrix_ID")
  pl <- plot_grid(pl1, pl2, ncol=2, align="hv")
  pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

  return(pl)
}

barPlot <- function(df) {
  palette <- c("#30110F", "#7F2C28", "#7F514E", "#CC4740", "#FF594F", "#FFA19C")
  pl <- ggplot(data=df, aes(x=pipeline, y=correlation, fill=pipeline)) +
    geom_bar(stat="identity") +
    facet_grid(~study) + 
    geom_text(aes(y=correlation, ymax=correlation, label=correlation),
              position = position_dodge(width=0.9), vjust=-.5, color="black") + 
    scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) + 
    theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
    scale_fill_manual(values=palette)
  return(pl)
}

scatterPlot <- function(df, title, lin, cub, ranges, labx="log2 probeset intensity", laby="log2 FPKM") {
  # TODO: replace the constants with passing arguments
  minx=ranges[1]
  maxx=ranges[2]
  miny=ranges[3]
  maxy=ranges[4]
  # Fit linear and polynomial functions:  
  r_lin <- summary(lin)$r.squared
  r_cub <- summary(cub)$r.squared
  
  commonTheme = list(labs(x=labx, y=laby, size=16),
                     theme_bw()+
                       theme(legend.position="none",
                             plot.margin=unit(c(5,5,5,5),"mm"),
                             axis.text = element_text(size=18),
                             axis.title = element_text(size=16),
                             plot.title = element_text(face="bold", size=15)))
  
  p <- ggplot(data=df,aes(x=eset,y=rnaseq)) + 
    geom_point(alpha=0.3) +
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#2FFF71", high="#FF4C48") +
    geom_smooth(method=lm, linetype=1, colour="black", se=F, size=2) + 
    geom_smooth(method=lm, formula=(y~x+I(x^2)+I(x^3)), linetype=1, colour="blue", se=F, size=2, fullrange=FALSE) + 
    geom_vline(xintercept = median(df$eset), linetype=2) + 
    geom_hline(yintercept = median(df$rnaseq), linetype=2) + 
    annotate("text", x=min(df$eset)+1.5, y=max(df$rnaseq)-1, label=paste("R=",round(r_lin, digits = 3)), color="black", size=6) +
    annotate("text", x=min(df$eset)+1.5, y=max(df$rnaseq)-2, label=paste("R=",round(r_cub, digits = 3)), color="blue", size=6) +
    scale_x_continuous(breaks = seq(minx, maxx, by = 2), limits=c(minx, maxx)) +
    scale_y_continuous(breaks = round(seq(round_any(min(df$rnaseq),2), max(df$rnaseq), by = 2), 1), limits=c(miny, maxy)) +  
    ggtitle(title) +
    coord_fixed() +
    commonTheme
  
  return(p)
}

scatterPlotVariance <- function(df, title, ranges) {
  # TODO: replace the constants with passing arguments
  minx=ranges[1]
  maxx=ranges[2]
  miny=ranges[3]
  maxy=ranges[4]
  
  commonTheme = list(labs(x="log2 FPKM",
                          y="SD of probesets groups SD", size=16),
                     theme_bw()+
                       theme(legend.position="none",
                             plot.margin=unit(c(5,5,5,5),"mm"),
                             axis.text = element_text(size=18),
                             axis.title = element_text(size=16),
                             plot.title = element_text(face="bold", size=18)))
  
  p <- ggplot(data=df,aes(x=rnaseq,y=eset)) + 
    geom_point(alpha=0.3) +
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#3E9EFF", high="#B23830") +
    geom_vline(xintercept = median(df$rnaseq), linetype=2) + 
    geom_hline(yintercept = median(df$eset), linetype=2) + 
    scale_x_continuous(breaks = seq(minx, maxx, by = 2), limits=c(minx, maxx)) +
    scale_y_continuous(breaks = round(seq(round_any(min(df$rnaseq),2), max(df$rnaseq), by = 2), 1), limits=c(miny, maxy)) +  
    ggtitle(title) +
    #    coord_fixed() +
    commonTheme
  
  return(p)
}

require(ggplot2)
diagPlot<-function(model){
  model <- fortify(model)
  p1<-ggplot(model, aes(.fitted, .resid))+geom_point(alpha=0.3)+
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#22B2AA", high="#FF9863")+
    stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")+
    xlab("Fitted values")+ylab("Residuals")+
    ggtitle("Residual vs Fitted Plot")+theme_bw()+theme(legend.position="none")
  
  p2<-ggplot(model, aes(.fitted, sqrt(abs(.stdresid))))+geom_point(na.rm=TRUE, alpha=0.3)+
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#22B2AA", high="#FF9863")+
    stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")+
    ylab(expression(sqrt("|Standardized residuals|")))+
    ggtitle("Scale-Location")+theme_bw()+theme(legend.position="none")
  
  #   p3<-ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE, alpha=0.3)+
  #       stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
  #       scale_fill_continuous(low="#22B2AA", high="#FF9863")+
  #       stat_smooth(method="loess", na.rm=TRUE)+
  #       xlab("Leverage")+ylab("Standardized Residuals")+
  #       ggtitle("Residual vs Leverage Plot")+
  #       scale_size_continuous("Cook's Distance", range=c(1,5))+
  #       theme_bw()+theme(legend.position="none")
  #   
  #   p4<-ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm=TRUE, alpha=0.3)+
  #       stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
  #       scale_fill_continuous(low="#22B2AA", high="#FF9863")+
  #       stat_smooth(method="loess", na.rm=TRUE)+
  #       xlab("Leverage hii")+ylab("Cook's Distance")+
  #       ggtitle("Cook's dist vs Leverage hii/(1-hii)")+
  #       geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")+
  #       theme_bw()+theme(legend.position="none")
  
  #return(list(rvfPlot=p1, sclLocPlot=p2, rvlevPlot=p3, cvlPlot=p4))
  return(list(rvfPlot=p1, sclLocPlot=p2))
}

partialCorPlot <- function(df, title) {
  df <- df[order(df$rnaseq),]
  
  step = 100
  window_size = trunc(nrow(df)/step)*5
  
  cor <- c()
  ind=1
  while (ind<(nrow(df)-window_size)) {
    sub <- df[ind:(ind+window_size),]
    cor <- c(cor, round(cor(sub$eset, sub$rnaseq, method="spearman"), 3))
    ind = ind + step
  }
  
  data <- data.frame(x=(1:length(cor)), y=cor)   
  
  commonTheme = list(labs(x="Rnaseq expression batch ->",
                          y="Spearman correlation", size=16),
                     theme_bw()+
                       theme(legend.position="none",
                             plot.margin=unit(c(5,5,5,5),"mm"),
                             axis.text = element_text(size=16),
                             plot.title = element_text(face="bold", size=18)))  
  
  pl <- ggplot(data, aes(x, y))+
    geom_smooth(span=0.5) +
    scale_y_continuous(limits=c(0, 0.3)) +
    commonTheme + 
    ggtitle(title)
}