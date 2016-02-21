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