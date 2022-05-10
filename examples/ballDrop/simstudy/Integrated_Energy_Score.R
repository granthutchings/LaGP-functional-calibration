source('../../mv_calib_lagp.R')
C_true = .1
g_true = 9.8
sd_true = 0.1 # sd of observation error - turns our this is very close to 5% of the variance of the sims at dist=5m, .15 is ~10%
# try making up a function where C,g have no interplay (additive)
t_at_h = function(C,h,R,g){
  return(acosh(exp(C*h/R))/sqrt(C*g/R))
}

load("~/Desktop/laGP Project/examples/ballDrop/simstudy/data/mvcDataSimStudy_50.RData")
load("~/Desktop/laGP Project/examples/ballDrop/simstudy/data/bd_ub_lagp_sim_results_m50.RData")
YindObs = mvcDataSimStudy[[1]]$YindObs
YindSim = mvcDataSimStudy[[1]]$YindSim
nsamples = 100
a = .025; b = .3
h = (b-a)/nsamples
k = seq(0,nsamples-1)
eval_points = c(a+k*h,b)

set.seed(11)
XobsNew = as.matrix(eval_points)
YobsNew = matrix(nrow=length(YindObs),ncol=nrow(XobsNew))
YobsTrue = matrix(nrow=length(YindObs),ncol=nrow(XobsNew))
for(j in 1:nrow(XobsNew)){
  YobsNew[,j] = t_at_h(C=C_true,h=YindObs,R=XobsNew[j,],g=g_true)# + rnorm(nrow(YobsNew),0,sd_true)
}

for(i in 1:1000){
  if(i %% 10 == 0){cat(i,' ')}
  XTdataNew = transform_xt(mvcDataSimStudy[[i]]$XTdata$sim$X$orig,
                           mvcDataSimStudy[[i]]$XTdata$sim$T$orig,
                           XobsNew)
  YdataNew = transform_y(mvcDataSimStudy[[i]]$Ydata$sim$orig,
                         YindSim,
                         YobsNew,
                         YindObs,
                         center = T,scale = T,scaletype = 'scalar')
  obsBasisNew = get_obs_basis(mvcDataSimStudy[[i]]$simBasis,YdataNew$obs$trans,YindSim,YindObs)
  SCinputsNew = get_SC_inputs(mvcDataSimStudy[[i]]$estLS,XTdataNew,nPC=mvcDataSimStudy[[i]]$simBasis$nPC)
  
  mvcDataNew = list(XTdata=XTdataNew,
                    Ydata=YdataNew,
                    simBasis=mvcDataSimStudy[[i]]$simBasis,
                    obsBasis=obsBasisNew,
                    estLS=mvcDataSimStudy[[i]]$estLS,
                    SCinputs=SCinputsNew,
                    bias=F)
  # get samples
  tsamp = matrix(sample(mcmc.out[[i]]$t.samp,100))
  # predict at eval_points
  yPred = ypred_mcmc(XpredOrig = matrix(eval_points),
                     mvcDataNew,tsamp,support='obs',returnSamples = T)
  # compute energy score at eval_points
  ES = energy_score(yPred$ysamp,mvcDataNew$Ydata$obs$orig)
  # compute integral over ES using trapezoidal rule
  #integrated_ES = h*(ES[1]/2 + sum(ES[2:(length(ES)-1)]) + ES[length(ES)]/2)
  integrated_ES = trapz(eval_points,ES)
}