panel.hist <- function(x,...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x[1:(length(x))], plot = FALSE)
  #abline(v=x[length(x)],col='red')
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1],y, ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}
panel.scatter <- function(x,y, ...){
  id = seq(1,(length(x)),length.out=500)
  #points(x[id],y[id],xlim=c(0,1),ylim=c(0,1))
  data = data.frame(x=x[id],y=y[id])
  ggplot(data, aes(x=x, y=y) ) + geom_density_2d()
  #points(x[length(x)],y[length(y)],xlim=c(0,1),ylim=c(0,1),col='red',pch=16)
}
panel.dens = function(x,y,...){
  contour(MASS::kde2d(x, y),nlevels = 10,xlim = c(0,1),ylim = c(0,1),add=T)
}
plot_pairs = function(mcmc.out,mvc.data,start=1,nlevels=10){
  df.mcmc = data.frame(mcmc.out$t.samp[start:nrow(mcmc.out$t.samp),])
  #df.mcmc$type = 'sample'
  #df.mcmc[(nrow(df.mcmc)+1),1:mvc.data$XT.data$p.t] = t(mvc.data$XT.data$obs$T$trans[1,]); df.mcmc[nrow(df.mcmc),(mvc.data$XT.data$p.t+1)] = 'truth'
  #df.mcmc$type = as.factor(df.mcmc$type)
  pairs(df.mcmc[,1:mvc.data$XT.data$p.t],
        diag.panel = panel.hist,
        upper.panel = NULL,
        lower.panel = panel.dens,
        xlim=c(0,1))
}
trace_plots = function(mcmc.out,start=1,rows=1,cols=1){
  par(mfrow=c(rows,cols))
  for(i in 1:ncol(mcmc.out$t.samp)){
    plot(mcmc.out$t.samp[start:nrow(mcmc.out$t.samp),i],type='l')
  }
}