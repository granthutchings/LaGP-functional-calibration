t_at_h = function(C,h,R,g){
  return(acosh(exp(C*h/R))/sqrt(C*g/R))
}
setwd('~/Desktop/laGP Project/examples/ballDrop/')
source("../../mv_calib_lagp.R")
load('simstudy/data/mvcDataSimStudy.RData')

m = 242

for(i in 1:1000){
  cat(i,' ')
  XTsim = create_oalhs(m,mvcDataSimStudy[[i]]$XTdata$pT+mvcDataSimStudy[[i]]$XTdata$pX,T,F)
  Xrange = c(.025,.3)
  Trange = c(.1-.05,.1+.05)
  Xsim=as.matrix(XTsim[,1]); Xsim = Xsim * (Xrange[2]-Xrange[1]) + Xrange[1]
  Tsim=as.matrix(XTsim[,2]); Tsim = Tsim * (Trange[2]-Trange[1]) + Trange[1]
  Ysim = matrix(nrow=length(mvcDataSimStudy[[i]]$YindSim),ncol=m)
  for(j in 1:m){
    Ysim[,j] = t_at_h(C=Tsim[j],h=mvcDataSimStudy[[i]]$YindSim,R=Xsim[j],g=9.8)
  }
  
  mvcDataSimStudy[[i]]$XTdata = transform_xt(Xsim,Tsim,mvcDataSimStudy[[i]]$XTdata$obs$X$orig,mvcDataSimStudy[[i]]$XTdata$obs$T$orig)
  mvcDataSimStudy[[i]]$Ydata = transform_y(Ysim,mvcDataSimStudy[[i]]$YindSim,mvcDataSimStudy[[i]]$Ydata$obs$orig,mvcDataSimStudy[[i]]$YindObs,
                                           center = T,scale = T, scaletype='scalar')
  mvcDataSimStudy[[i]]$simBasis = get_basis(mvcDataSimStudy[[i]]$Ydata$sim$trans,pctVar = .99)
  mvcDataSimStudy[[i]]$obsBasis = get_obs_basis(mvcDataSimStudy[[i]]$simBasis,mvcDataSimStudy[[i]]$Ydata$obs$trans,
                                                mvcDataSimStudy[[i]]$YindSim,mvcDataSimStudy[[i]]$YindObs)
  mvcDataSimStudy[[i]]$estLS = mv_lengthscales(mvcDataSimStudy[[i]]$XTdata,mvcDataSimStudy[[i]]$simBasis$Vt,g=1e-7)
  mvcDataSimStudy[[i]]$SCinputs = get_SC_inputs(mvcDataSimStudy[[i]]$estLS,mvcDataSimStudy[[i]]$XTdata,mvcDataSimStudy[[i]]$simBasis$nPC)
  mvcDataSimStudy[[i]]$bias = F
  mvcDataSimStudy[[i]]$design_score = sum(dist(XTsim,upper=T)^(-2))
}
save(mvcDataSimStudy,file=paste0('simstudy/mvcDataSimStudy_',m,'.RData'))
py_save_object(mvcDataSimStudy,file=paste0('simstudy/mvcDataSimStudy_',m,'.pkl'))
