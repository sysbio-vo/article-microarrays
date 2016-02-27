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