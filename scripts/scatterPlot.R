scatterPlot <- function(df, title, lin, cub, ranges) {
  # TODO: replace the constants with passing arguments
  minx=ranges[1]
  maxx=ranges[2]
  miny=ranges[3]
  maxy=ranges[4]
  # Fit linear and polynomial functions:  
  r_lin <- summary(lin)$r.squared
  r_cub <- summary(cub)$r.squared
  
  commonTheme = list(labs(x="log2 probeset intensity",
                          y="log2 FPKM", size=16),
                     theme_bw()+
                     theme(legend.position="none",
                           plot.margin=unit(c(5,5,5,5),"mm"),
                           axis.text = element_text(size=18),
                           axis.title = element_text(size=16),
                           plot.title = element_text(face="bold", size=18)))
  
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