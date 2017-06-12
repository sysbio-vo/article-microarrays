getDEGS <- function(meta.vars, pheno.data, exprs, noannotate=FALSE) {
  # Subset groups for comparison
  ind <- c()
  for (i in meta.vars) {
    ind <- c(ind, which(pheno.data$CancerType==i))
  }
  pdata.short <- pheno.data[ind, ]
  exprs.short <- exprs[,which(colnames(exprs) %in% pdata.short$SampleAccessionNumber)]
  exprs.short <- exprs.short[,order(match(colnames(exprs.short), rownames(pdata.short)))]
  
  # Create design matrix
  design = model.matrix(~factor(pdata.short$CancerType, levels=meta.vars), data=pdata.short)
  colnames(design) <- c(meta.vars[1], paste(meta.vars[1], "vs", meta.vars[2], sep=""))
  # Fit with linear models
  fit <- lmFit(exprs.short, design)
  fit <- eBayes(fit)
  # Get all the genes with logFC, p-values, no filtering
  degs <- topTable(fit, coef=paste(meta.vars[1], "vs", meta.vars[2], sep=""), adjust.method="fdr", number=nrow(fit))

  # Merge degs with expression matrix
  exprs.degs <- merge(degs, exprs.short, by="row.names")
  if (noannotate) {
    rownames(exprs.degs) <- exprs.degs[, 1]
    exprs.degs <- exprs.degs[, -1]
  } else {
    colnames(exprs.degs)[1] <- "ENTREZID"
    # Add information about gene names
    EntrezID_Symbol<-select(org.Hs.eg.db, exprs.degs$ENTREZID, c("SYMBOL", "GENENAME"))
    exprs.degs <- cbind(EntrezID_Symbol, exprs.degs)
    exprs.degs <- exprs.degs[,-4]
  }
  
  return(exprs.degs)
}

filterDEGS <- function(degs, pval, fc, adj) {
  # Filter by p-values
  if (missing(adj)) {
    degs <- degs[degs$adj.P.Val < pval,]
  } else if (adj==FALSE) {
    degs <- degs[degs$P.Value < pval,]
  } else if (adj==TRUE) {
    degs <- degs[degs$adj.P.Val < pval,]
  }
  
  # Sort by logFC
  degs <- degs[order(abs(degs$logFC), decreasing = TRUE),]
  # Filter by logFC
  degs <- degs[abs(degs$logFC) > fc,]  
  return (degs)
}
