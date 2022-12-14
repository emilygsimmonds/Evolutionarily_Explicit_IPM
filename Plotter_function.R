# Script with plotting function

#-----------------------------------------------------------------------------

plotter <- function(PS.store, y_label, x_label, main_label, x.axis=TRUE, y.axis=TRUE, adj,
                    ylim=c(0,450)){
  plot(1:91, PS.store[[1]]$PS, type="l", xlab="", ylab="", 
       main="",ylim=ylim, axes=F, col='grey', las=1, cex.lab=1.25)
  if(y.axis == TRUE){axis(2, las=1, cex.axis=1.25)}
  if(x.axis==TRUE){axis(1, at=c(1,91), labels=c("",2100), cex.axis=1.25)}
  title(main=main_label, adj=adj, cex.main=1.5)
  title(xlab=x_label, line=2.5, cex.lab=1.25)
  title(ylab=y_label, line=3.5, cex.lab=1.25)
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){
    test <- which(PS.store[[i]]$PS<1)
    if(length(test)>0){lines(1:test[1], PS.store[[i]]$PS[1:test[1]], type='l', col=alpha("grey",0.2))
      extinctions[i,] <- c(test[1], PS.store[[i]]$PS[test[1]])}
      else{lines(1:91, PS.store[[i]]$PS, type='l', col=alpha("grey",0.2))}}
  points(extinctions[,1], extinctions[,2], pch=4, cex=1, col=alpha("black",1))
  abline(h=131, col="black", lty=1)
  abline(h=219, col="black", lty=2)
  abline(h=40, col="black", lty=2)
  
}
