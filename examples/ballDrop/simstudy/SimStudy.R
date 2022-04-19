library(parallel)
library(laGP)
source("../../mv_calib_lagp.R")


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
load('mvcDataSimStudy.RData')
for(ii in 1:N){
  mvcDataSimStudy[[ii]]$bias=F
  prop.step[ii] = tune_step_sizes(mvcDataSimStudy[[ii]],50,20)
  mcmc.out[[ii]] = mcmc(mvcDataSimStudy[[ii]],tInit=.5,bias=F,nsamples=2500,nburn=1000,prop.step=prop.step[ii],verbose = T)

  MSE[ii] = mean((mcmc.out[[ii]]$t.samp-tStar[ii])^2)
  AbsErrMean[ii] = abs(mean(mcmc.out[[ii]]$t.samp)-tStar[[ii]])
  MeanAbsErr[ii] = mean(abs(mcmc.out[[ii]]$t.samp-tStar[[ii]]))

  acc.ratio[ii] = mcmc.out[[ii]]$acpt.ratio
  cred.int[ii,] = quantile(mcmc.out[[ii]]$t.samp,c(.025,.975))
  cover[ii] = ifelse(tStar[ii]>=cred.int[ii,1] & tStar[ii]<=cred.int[ii,2],1,0)
}

save(tStar,N,MSE,AbsErrMean,MeanAbsErr,acc.ratio,mcmc.out,cred.int,cover,prop.step,file='bd_ub_mcmc_simstudy.RData')
