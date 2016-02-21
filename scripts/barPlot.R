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
  