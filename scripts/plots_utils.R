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