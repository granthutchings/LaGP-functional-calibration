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
plot_pairs = function(mcmc.out,mvc.data,start=1){
  df.mcmc = data.frame(mcmc.out$t.samp[start:nrow(mcmc.out$t.samp),])
  #df.mcmc$type = 'sample'
  #df.mcmc[(nrow(df.mcmc)+1),1:mvc.data$XT.data$p.t] = t(mvc.data$XT.data$obs$T$trans[1,]); df.mcmc[nrow(df.mcmc),(mvc.data$XT.data$p.t+1)] = 'truth'
  #df.mcmc$type = as.factor(df.mcmc$type)
  pairs(df.mcmc[,1:mvc.data$XT.data$p.t],
        diag.panel = panel.hist,
        lower.panel = NULL,
        xlim=c(0,1))
}
trace_plots = function(mcmc.out,start=1,rows=1,cols=1){
  par(mfrow=c(rows,cols))
  for(i in 1:ncol(mcmc.out$t.samp)){
    plot(mcmc.out$t.samp[start:nrow(mcmc.out$t.samp),i],type='l')
  }
}

design = as.matrix(read.table("data/design.txt",header = T)[-c(7,8,10,11)])
sim104 = as.matrix(10000*read.csv("data/ysim_104.csv")[,c(4,6,8,10)])
field104 <- as.matrix(read.csv("data/yobs_104.csv")[c(4,6,8,10)])
sim105 = as.matrix(10000*read.csv("data/ysim_105.csv")[,c(4,6,8,10)])
field105 <- as.matrix(read.csv("data/yobs_105.csv")[c(4,6,8,10)])
sim106 = as.matrix(10000*read.csv("data/ysim_106.csv")[,c(4,6,8,10)])
field106 <- as.matrix(read.csv("data/yobs_106.csv")[c(4,6,8,10)])

sim = cbind(sim104,sim105,sim106)
field = cbind(field104,field105,field106)

mvc.data = mvc_data(T.sim = design,
                    Y.sim = t(sim),
                    Y.obs = t(field),
                    y.ind.sim = seq(1,12), y.ind.obs = seq(1,12),
                    pct.var = .99, small=T, precomp = T, sc.subsample = F,
                    scaletype = 'scalar')
init = calib_init(as.matrix(expand.grid(rep(list(c(.25,.5,.75)),mvc.data$XT.data$p.t))),mvc.data,sample=F)
calib = mv_calib(mvc.data,init,nrestarts=1,optimize=T,sample=F)
prop.step = tune_step_sizes(mvc.data,n.samples=50,n.levels=20,target.accept.rate = .4,min.step = .01, max.step=1)
mcmc.out = do_mcmc(mvc.data,t.init=calib$theta.hat,n.samples=10000,n.burn=0,
                   prop.step=prop.step,verbose=T,adapt=T)
save(mvc.data,init,calib,prop.step,mcmc.out,file='results/data/calib_all.RData')

pdf('results/plots/calib_all_pairs.pdf')
plot_pairs(mcmc.out,mvc.data,start=1000)
dev.off()
pdf('results/plots/calib_all_trace.pdf')
trace_plots(mcmc.out,start=1,4,3)
dev.off()