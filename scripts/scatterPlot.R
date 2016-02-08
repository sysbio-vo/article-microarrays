require(ggplot2)
scatterPlot <- function(df) {
  # Fit linear and polynomial functions:  
  lin<-lm(rnaseq~eset)
  quad<-lm(rnaseq~eset+I(eset^2)+I(eset^3))
  r_lin <- summary(lin)$r.squared
  r_quad <- summary(quad)$r.squared
  
  commonTheme = list(labs(x="log2 probeset intensity",
                          y="log2 FPKM", size=16),
                     theme_bw()+
                     theme(legend.position="none",
                           plot.margin=unit(c(5,5,5,5),"mm"),
                           axis.text = element_text(size=18),
                           axis.title = element_text(size=16),
                           plot.title = element_text(face="bold", size=22)))
  
  p <- ggplot(data=df,aes(x,y)) + 
    geom_point(alpha=0.3) +
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#2FFF71", high="#FF4C48") +
    geom_smooth(method=lm, linetype=1, colour="black", se=F, size=2) + 
    geom_smooth(method=lm, formula=y~x+I(x^2)+I(x^3), linetype=1, colour="blue", se=F, size=2) + 
    geom_vline(xintercept = median(df$x), linetype=2) + 
    geom_hline(yintercept = median(df$y), linetype=2) + 
    annotate("text", x=min(df$x)+1.5, y=max(df$y), label=paste("R=",round(r_lin, digits = 3)), color="black", size=6) +
    annotate("text", x=min(df$x)+1.5, y=max(df$y)-1, label=paste("R=",round(r_quad, digits = 3)), color="blue", size=6) +
    scale_x_continuous(breaks = round(seq(round_any(min(df$x),2), max(df$x), by = 2), 1)) +
    scale_y_continuous(breaks = round(seq(round_any(min(df$y),2), max(df$y), by = 2), 1)) +  
    ggtitle(paste(studies$ID[affy[ind]],pipe_type)) +
    coord_fixed() +
    commonTheme

  return(p)
}