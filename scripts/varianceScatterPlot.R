scatterPlot <- function(df, title, ranges) {
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