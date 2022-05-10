library(parallel)
library(laGP)
setwd('~/Desktop/laGP Project/examples/ballDrop/simstudy/')
source("../../../mv_calib_lagp.R")

# MCMC sim study using the data saved during the point estimate study

N = 1000
MSE = numeric(N)
AbsErrMean = numeric(N)
MeanAbsErr = numeric(N)
acc.ratio = numeric(N)
prop.step = numeric(N)
mcmc.out = vector(mode='list',length=N)
cred.int = matrix(nrow=N,ncol=2)
cover = numeric(N)
load('data/tStar.RData')
load('data/mvcDataSimStudy_50.RData')
for(ii in 1:N){
  cat(ii)
  if(ii==1){
    ptm = proc.time()
  }
  prop.step[ii] = tune_step_sizes(mvcDataSimStudy[[ii]],50,20,end=50)
  mcmc.out[[ii]] = mcmc(mvcDataSimStudy[[ii]],tInit=.5,nsamples=2500,nburn=1000,prop.step=prop.step[ii],verbose = F,end=50)
  if(ii==1){
    mcmc.time = ptm - proc.time()
  }
  MSE[ii] = mean((mcmc.out[[ii]]$t.samp-tStar[ii])^2)
  AbsErrMean[ii] = abs(mean(mcmc.out[[ii]]$t.samp)-tStar[[ii]])
  MeanAbsErr[ii] = mean(abs(mcmc.out[[ii]]$t.samp-tStar[[ii]]))

  acc.ratio[ii] = mcmc.out[[ii]]$acpt.ratio
  cred.int[ii,] = quantile(mcmc.out[[ii]]$t.samp,c(.025,.975))
  cover[ii] = ifelse(tStar[ii]>=cred.int[ii,1] & tStar[ii]<=cred.int[ii,2],1,0)
}

save(tStar,N,MSE,AbsErrMean,MeanAbsErr,acc.ratio,mcmc.out,mcmc.time,cred.int,cover,prop.step,file='lagp_m50_results.RData')
