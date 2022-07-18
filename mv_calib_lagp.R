library(Matrix)
library(tidyverse)
library(laGP)
library(tgp)
library(parallel)
library(pracma)
library(mvtnorm)
library(crs)

logit = function(x){
  log(x/(1-x))
}
inv.logit = function(x){
  exp(x)/(1+exp(x))
}

# ypred_samp n.samples x n.y x n
energy_score = function(y.samp,y,terms=F){
  n.samples = dim(y.samp)[1]
  n = dim(y.samp)[3]
  
  get_es = function(i,terms){
    diff = t(t(y.samp[,,i])-y[,i])
    error_term = mean(apply(diff,1,norm,type='2')) # mean over samples of norm of difference
    uq_term = (1/(2*n.samples^2))*sum(sapply(1:n.samples, function(k) 
      sum(sapply(1:n.samples, function(j)
        norm(y.samp[k,,i]-y.samp[j,,i],type='2')))))
    if(terms){
      return(list(error_term=error_term,uq_term=uq_term))
    } else{
      return(error_term - uq_term)
    }
    
  }
  if(terms){
    es = matrix(unlist(mclapply(1:n, function(i) get_es(i,terms))),ncol=2,byrow=T)
  } else{
    es = unlist(
      mclapply(1:n, function(i) get_es(i,terms))
    )
  }
  
  return(es)
}
interval_score = function(y.samp,y,terms=F){
  n.samples = dim(y.samp)[1]
  n = dim(y.samp)[3]
  ci = apply(y.samp,2:3,quantile,c(.025,.975))
  
  width = array(dim=dim(y.samp)[2:3])
  penalty = array(dim=dim(y.samp)[2:3])
  i_score = array(dim=dim(y.samp)[2:3])
  alpha = .05
  for(k in 1:n){ # loop over pred locations
    # i_score = (u-l)+2/alpha(l-x)I(x<l)+2/alpha(x-u)I(x>u)
    width[,k] = (ci[2,,k]-ci[1,,k])
    penalty[,k] = (2/alpha)*(ci[1,,k] - y[,k])*as.integer(y[,k]<ci[1,,k]) +
      (2/alpha)*(y[,k] - ci[2,,k])*as.integer(y[,k]>ci[2,,k])
    i_score[,k] = width[,k] + penalty[,k]
  }
  
  return(mean(i_score))
}

plot_wz_pairs = function(w,z=NULL,legend=T){
  n.pc = nrow(w)
  m = ncol(w)
  if(!is.null(z)){
    n = ncol(z)
    if(n.pc>1){
      pairs(rbind(t(w),t(z)),col=c(rep('blue',m),rep('orange',n)),pch=1,labels=paste0('PC',1:n.pc),
            oma=c(3,3,3,15))
      if(legend){
        par(xpd = TRUE)
        legend("bottomright", pch=1, col=c('blue','orange'), legend=c('w','z'))
      }
    } else{
      hist(t(w),col='blue');abline(v=t(z),col='orange')
      if(legend){
        legend("topright", pch=1, col=c('blue','orange'), legend=c('w','z'))
      }
    }
  } else{
    if(n.pc>1){
      pairs(t(w),col=rep('blue',m),pch=1,labels=paste0('PC',1:n.pc),
            oma=c(3,3,3,15))
      if(legend){
        par(xpd = TRUE)
        legend("bottomright", pch=1, col=c('blue'), legend=c('w'))
      }
    } else{
      hist(t(w),col='blue');
      if(legend){
        legend("topright", pch=1, col=c('blue'), legend=c('w'))
      }
    }
  }
}

w_to_y = function(w,B,mean=NULL,sd=NULL){
  Y = B%*%w
  if(!is.null(sd)){
    Y = Y * sd
  }
  if(!is.null(mean)){
    Y = Y + mean
  }
  return(Y)
}

# Returns predictions at X.pred.orig using the calibration parameters theta. To do this
# the emulator must be 'fit' at theta.
predict_w = function(X.pred.orig,mvc.data,theta=NULL,start=6,end=50,sample=F,n.samples=1){
  
  n.pc = nrow(mvc.data$sim.basis$V.t)
  y=NULL
  
  if(!is.null(theta)){
    X = transform_xt(X.sim = mvc.data$XT.data$sim$X$orig,
                     T.sim = mvc.data$XT.data$sim$T$orig,
                     X.obs = X.pred.orig,
                     T.obs = matrix(rep(theta,nrow(X.pred.orig)),nrow=nrow(X.pred.orig),byrow=T))
    X.pred.sc = get_SC_inputs(mvc.data$est.ls,X,n.pc)
    w = aGPsep_SC_mv(X=X.pred.sc$XT.sim,
                     Z=lapply(1:n.pc,function(i) mvc.data$sim.basis$V.t[i,]),
                     XX=lapply(1:n.pc,function(i) cbind(X.pred.sc$X.obs[[i]],X.pred.sc$T.obs[[i]])),
                     start=start,
                     end=end)
  } else{
    X = transform_xt(X.sim = mvc.data$XT.data$sim$X$orig,
                     X.obs = X.pred.orig)
    X.pred.sc = get_SC_inputs(mvc.data$est.ls,X,n.pc)
    w = aGPsep_SC_mv(X=X.pred.sc$XT.sim,
                     Z=lapply(1:n.pc,function(i) mvc.data$sim.basis$V.t[i,]),
                     XX=lapply(1:n.pc,function(i) X.pred.sc$X.obs[[i]]),
                     start=start,
                     end=end)
  }
  
  if(sample){
    # if(n.samples>1){
    #   w$sample = array(rmvnorm(n.samples,as.numeric(w$mean),diag(as.numeric(w$var),nrow=n.pc*nrow(X.pred.orig))),dim=c(n.samples,dim(w$mean)))
    # } else{
    #   w$sample = t(matrix(rmvnorm(n.samples,as.numeric(w$mean),diag(as.numeric(w$var),nrow=n.pc*nrow(X.pred.orig))),ncol=n.pc,byrow=T))
    # }
    w$sample = rmvnorm(n.samples,as.numeric(w$mean),diag(as.numeric(w$var)))
    dim(w$sample) = c(n.samples,dim(w$mean)); w$sample = drop(w$sample)
  } else{
    w$sample=w$mean
  }
  return(w)
}
predict_y = function(X.pred.orig,mvc.data,n.samples=1){

  n = nrow(X.pred.orig)
  n.y = nrow(mvc.data$sim.basis$B)
  
  w = predict_w(X.pred.orig,mvc.data,sample=T,n.samples=n.samples)

  y.samp = array(dim=c(n.y,n.samples,n))
  for(i in 1:n.samples){
      y.samp[,i,] = w_to_y(w$sample[i,,],mvc.data$sim.basis$B,mvc.data$Y.data$sim$mean,mvc.data$Y.data$sim$sd)
  }
  y.mean = apply(y.samp,c(1,3),mean)
  y.var = apply(y.samp,c(1,3),var)
  y.int = lapply(1:n,function(i) t(apply(y.samp[,,i],1,quantile,c(.025,.975))))
  return(list(mean=y.mean,var=y.var,conf.int=y.int,samp=y.samp))
    
  # if(n.samples>1){
  #   w$sample = rmvnorm(n.samples,as.numeric(w$mean),diag(as.numeric(w$var))); dim(w$sample) = c(n.samples,dim(w$mean))
  #   if(!is.null(v)){
  #     v$sample = rmvnorm(n.samples,as.numeric(v$mean),diag(as.numeric(v$var))); dim(v$sample) = c(n.samples,dim(v$mean))
  #   }
  #   y.samp = array(dim=c(n.y,n.samples,n))
  #   for(i in 1:n.samples){
  #       y.samp[,i,] = w_to_y(w$sample[i,,],B)
  #       if(!is.null(v)){
  #         y.samp[,i,] = y.samp[,i,] + w_to_y(v$sample[i,,],D)
  #       }
  #   }
  #   y.mean = apply(y.samp,c(1,3),mean)
  #   y.var = apply(y.samp,c(1,3),var)
  #   y.int = lapply(1:n,function(i) t(apply(y.samp[,,i],1,quantile,c(.025,.975))))
  #   return(list(mean=y.mean,var=y.var,conf.int=y.int,samp=y.samp))
  # } else{
  #   y = w_to_y(w$mean,B)
  #   if(!is.null(v)){
  #     y = y + w_to_y(v$mean,D)
  #   }
  #   return(y)
  # }
}

# Returns predictions from the bias model at X.pred.orig. n.pc GP objects must be given in a list
mv_delta_predict = function(X.pred.orig,GP.list,mvc.data,sample=F,n.samples=1){
  
  X.pred.std = unit_xform(X.pred.orig,X.min=mvc.data$XT.data$sim$X$min,X.range=mvc.data$XT.data$sim$X$range)$trans
  
  n.pc = ncol(mvc.data$obs.basis$D)
  tmp = lapply(1:n.pc, function(i) predGPsep(GP.list[[i]],X.pred.std,lite=T))
  v = NULL
  for(i in 1:n.pc){
    v$mean = rbind(v$mean,tmp[[i]]$mean)
    v$var = rbind(v$var,tmp[[i]]$s2)
  }
  if(sample){
    v$sample = rmvnorm(n.samples,as.numeric(v$mean),diag(as.numeric(v$var)))
    dim(v$sample) = c(n.samples,dim(v$mean)); v$sample = drop(v$sample)
  } else{
    v$sample=v$mean
  }
  return(v)
}

# FUNCTION: ypred_mle: get y predictions using MLE of t*
{## Parameters:
  # required
  # X.pred.orig : npred x px matrix containing the X locations for prediction on the original scale
  # mvc.data : Data object
  # theta: pt vector of MLE's of t*
  # ssq: sigma^2 associated with MLE theta
  # optional
  # n.samples: int, number of samples to be returned from the sampling distribution of y, if n.samples=1 the mean is returned
  # support: string, should samples of y be returned or just the mean and uncertainty?
  # start: int, initial laGP neighborhood size
  # end: int, final laGP neighborhood size
  ## Returns:
  # list containing y predictions on original scale
}
ypred_mle = function(X.pred.orig,mvc.data,theta,ssq,n.samples=1,support='obs',start=6,end=50,bias_GP_list=NULL){
  
  if(support=='obs'){
    B = mvc.data$obs.basis$B
    D = mvc.data$obs.basis$D
    ym = mvc.data$Y.data$obs$mean
    ysd = mvc.data$Y.data$obs$sd
    n.y = nrow(mvc.data$Y.data$obs$orig)
  } else{
    B = mvc.data$sim.basis$B
    D = mvc.data$sim.basis$D
    ym = mvc.data$Y.data$sim$mean
    ysd = mvc.data$Y.data$sim$sd
    n.y = nrow(mvc.data$Y.data$sim$orig)
  }
  n = nrow(X.pred.orig)
  
  # emulator predictions
  w = predict_w(X.pred.orig,mvc.data,theta,start=start,end=end,sample=T,n.samples)

  if(!mvc.data$bias){
    # unbiased prediction
    y.samp = array(dim=c(n.y,n.samples,n))
    for(i in 1:n){
      if(n.samples>1){
        y.samp[,,i] = t(mvtnorm::rmvnorm(n.samples,mean=B%*%w$mean[,i,drop=F],
                                        sigma = ssq*eye(n.y) + B%*%diag(w$var[,i])%*%t(B))) *
          ysd + ym
      } else{
        y.samp = B%*%w$mean[,i,drop=F] * ysd + ym
      }
      # for(j in 1:n.samples){
      #   y.samp[,j,i] = mvtnorm::rmvnorm(1,mean=B%*%w$sample[j,,i,drop=F],sigma = ssq*eye(n.y)) * ysd + ym
      # }
    }
    y.mean = apply(y.samp,c(1,3),mean)
    y.var = apply(y.samp,c(1,3),var)
    y.int = lapply(1:n,function(i) t(apply(y.samp[,,i],1,quantile,c(.025,.975))))
    return(list(y.samp=y.samp,y.mean=y.mean,y.var=y.var,y.int=y.int))
  } else{
    # biased prediction
    eta.samp = array(dim=c(n.y,n.samples,n))
    delta.samp = array(dim=c(n.y,n.samples,n))
    v = mv_delta_predict(X.pred.orig,bias_GP_list,mvc.data,F)
    for(j in 1:n){
      eta.samp[,,j] = t(mvtnorm::rmvnorm(n.samples,mean=B%*%w$mean[,j,drop=F],
                                       sigma = B%*%diag(w$var[,j])%*%t(B))) *
        ysd + ym
      delta.samp[,,j] = t(mvtnorm::rmvnorm(n.samples,mean=D%*%v$mean[,j,drop=F],
                                         sigma = D%*%diag(v$var[,j],nrow=length(v$var[,j]))%*%t(D) +
                                           ssq*eye(n.y) )) * ysd
    }

    y.samp = eta.samp + delta.samp
    y.mean = apply(y.samp,c(1,3),mean)
    y.var = apply(y.samp,c(1,3),var)
    eta.mean = apply(eta.samp,c(1,3),mean)
    delta.mean = apply(delta.samp,c(1,3),mean)
    y.int = aperm(array(sapply(1:n,function(i) t(apply(y.samp[,,i],1,quantile,c(.025,.975)))),dim=c(n.y,2,n)),c(2,1,3))
    return(list(y.samp=y.samp,mean=y.mean,var=y.var,conf.int=y.int,
                eta.mean=eta.mean,delta.mean=delta.mean))
  }
}

# FUNCTION: get_y_pred: get y predictions using posterior samples of t*
{## Parameters:
  # required
  # X.pred.orig : npred x px matrix containing the X locations for prediction on the original scale
  # mvc.data : Data object
  # t.samp: n.samples x pt matrix containing posterior samples of t* on the Original scale (mcmc function returns samples on the standardized scale)
  # ssq.samp: n.samples vector of posterior samples of sigma^2
  # optional
  # support: string, one of 'obs' or 'sim' indicating if predictions should be returned at the indices YindObs or YindSim
  # return.samples: boolean, should samples of y be returned or just the mean and uncertainty?
  # start: int, initial laGP neighborhood size
  # end: int, final laGP neighborhood size
  ## Returns:
  # list containing y predictions on original scale
}

# this function allows prediction at an X that was not passed to do_mcmc, but does require GP objects from do_mcmc
get_y_pred = function(X.pred.orig, mvc.data ,mcmc.out, samp.ids=0, support='obs',return.samples=F,start=6,end=50){
  if(support=='obs'){
    B = mvc.data$obs.basis$B
    D = mvc.data$obs.basis$D
    ym = mvc.data$Y.data$obs$mean
    ysd = mvc.data$Y.data$obs$sd
    n.y = nrow(mvc.data$Y.data$obs$orig)
  } else{
    B = mvc.data$sim.basis$B
    D = mvc.data$sim.basis$D
    ym = mvc.data$Y.data$sim$mean
    ysd = mvc.data$Y.data$sim$sd
    n.y = nrow(mvc.data$Y.data$sim$orig)
  }
  bias = mvc.data$bias
  if(samp.ids==0){
    samp.ids=mcmc.out$pred.samp.ids
  }
  if(bias & !all(samp.ids %in% mcmc.out$pred.samp.ids)){
    cat('For biased prediction samp.ids must be a subset of mcmc.out$pred.samp.ids. Setting samp.ids = mcmc.out$pred.samp.ids')
    break
  }
  n.samples = length(samp.ids)
  n.pred.X = nrow(X.pred.orig)
  t.samp = as.matrix(mcmc.out$t.samp)[samp.ids,,drop=F]
  ssq.samp = mcmc.out$ssq.samp[samp.ids]
  eta.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  delta.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  y.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  # put t.samp on original scale
  t.samp = t(t(t.samp) * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min)
  for(i in 1:n.samples){
    w = predict_w(X.pred.orig,mvc.data,t.samp[i,],sample=T,start=start,end=end)
    eta.samp[i,,] = B%*%w$sample
    if(bias){
      v = mv_delta_predict(X.pred.orig,mcmc.out$v.GPs[[i]],mvc.data,sample=T,n.samples=1)
      delta.samp[i,,] = D%*%drop(v$sample)
    }
    sigma = ssq.samp[i]*eye(n.y)
    for(j in 1:n.pred.X){
      y.samp[i,,j] = rmvnorm(1,mean=eta.samp[i,,j]+delta.samp[i,,j],sigma=sigma) * ysd + ym
    }
  }
  
  y.mean = apply(y.samp,2:3,mean)
  eta.mean = apply(eta.samp,2:3,mean)
  eta.conf.int = apply(eta.samp,2:3,quantile,c(.025,.975))
  delta.mean = apply(delta.samp,2:3,mean)
  delta.conf.int = apply(delta.samp,2:3,quantile,c(.025,.975))
  y.var = apply(y.samp,2:3,var)
  conf.int = apply(y.samp,2:3,quantile,c(.025,.975))
  returns = list(mean=y.mean,var=y.var,conf.int=conf.int,
                 eta.mean=eta.mean,delta.mean=delta.mean,
                 eta.conf.int=eta.conf.int,delta.conf.int=delta.conf.int,
                 ssq = ssq.samp)
  if(return.samples){
    returns$y.samp = y.samp
    returns$eta.samp = eta.samp
    returns$delta.samp = delta.samp
  }
  return(returns)
}

# this functions gets predictions of y using samples of w,v from MCMC. This function is 
# only applicable if X.pred was passed to do_mcmc.
get_y_pred_mcmc = function(mvc.data,mcmc.out,support='obs',return.samples=F){
  if(support=='obs'){
    B = mvc.data$obs.basis$B
    D = mvc.data$obs.basis$D
    ym = mvc.data$Y.data$obs$mean
    ysd = mvc.data$Y.data$obs$sd
    n.y = nrow(mvc.data$Y.data$obs$orig)
  } else{
    B = mvc.data$sim.basis$B
    D = mvc.data$sim.basis$D
    ym = mvc.data$Y.data$sim$mean
    ysd = mvc.data$Y.data$sim$sd
    n.y = nrow(mvc.data$Y.data$sim$orig)
  }
  bias = mvc.data$bias
  n.pred.X = dim(mcmc.out$w.pred)[3]
  n.samples = dim(mcmc.out$w.pred)[1]
  y.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  eta.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  delta.samp = array(0,dim=c(n.samples,n.y,n.pred.X))
  
  for(i in 1:n.samples){
    eta.samp[i,,] = B%*%matrix(mcmc.out$w.pred[i,,],nrow=dim(mcmc.out$w.pred)[2])
    if(bias){delta.samp[i,,] = D%*%mcmc.out$v.pred[i,,]}
    sigma = mcmc.out$ssq.pred[i]*eye(n.y)
    for(j in 1:n.pred.X){
      y.samp[i,,j] = rmvnorm(1,mean=eta.samp[i,,j]+delta.samp[i,,j],sigma=sigma) * ysd + ym
    }
  }
  y.mean = apply(y.samp,2:3,mean)
  y.var = apply(y.samp,2:3,var)
  conf.int = apply(y.samp,2:3,quantile,c(.025,.975))
  
  eta.mean = apply(eta.samp,2:3,mean)
  eta.conf.int = apply(eta.samp,2:3,quantile,c(.025,.975))
  
  delta.mean = apply(delta.samp,2:3,mean)
  delta.conf.int = apply(delta.samp,2:3,quantile,c(.025,.975))
  
  returns = list(mean=y.mean,var=y.var,conf.int=conf.int,
                 eta.mean=eta.mean,delta.mean=delta.mean,
                 eta.conf.int=eta.conf.int,delta.conf.int=delta.conf.int,
                 ssq = mcmc.out$ssq.pred)
  if(return.samples){
    returns$y.samp = y.samp
    returns$eta.samp = eta.samp
    returns$delta.samp = delta.samp
  }
  return(returns)
}

# FUNCTION: multivariate aGPsep for use with stretched and compressed inputs only
{## Parameters:
# required
# X : list of length n.pc which contains stretched and compressed training input matrices
# Z : list of length n.pc which contains training response vectors
# XX: list of length n.pc with contains stretched and compressed prediction input matrices
# optional
# start: integer, initial size of laPG neighborhood
# end: integer, final size of laGP neighborhood
# g: double, laGP nugget
## Returns:
# lagp_fit: list of n.pc outputs from aGPsep
# mean: matrix (n.pc x n) of prediction means
# var: matrix (n.pc x n) of prediction variances
## Calls:
# aGPsep from laGP package
# w_to_y to transform w -> y
## Called by:
# mv_calib.bias: biased calibration
# mv_calib.nobias: unbiased calibration
# mv.em.pred: laGP emulator prediction function at new inputs
}
aGPsep_SC_mv = function(X, Z, XX, start=6, end=50, g=1/10000, bias=F, sample=F){
  n.pc = length(Z)
  n.y = nrow(XX[[1]])

  if(bias){
    # if we use laGP for the bias, we cannot use stretched an compressed inputs anymore
    lagp_fit = lapply(1:n.pc,function(i) aGPsep(X = X[[i]],
                                               Z = Z[[i]],
                                               XX = XX[[i]],
                                               start = start,
                                               end = min(end,length(Z[[i]])-1),
                                               method = 'alc',
                                               verb=0))
  } else{
    lagp_fit = lapply(1:n.pc,function(i) aGPsep(X = X[[i]],
                                               Z = Z[[i]],
                                               XX = XX[[i]],
                                               d = list(mle = FALSE, start = 1),
                                               g = g,
                                               start = start,
                                               end = min(end,length(Z[[i]])-1),
                                               method = 'nn',
                                               verb=0))
  }

  mean = array(dim=c(n.pc,n.y))
  var = array(dim=c(n.pc,n.y))
  for (i in 1:n.pc) {
    mean[i,] = lagp_fit[[i]]$mean
    var[i,] = lagp_fit[[i]]$var
  }
  if(sample){
    sample = rmvnorm(1,as.numeric(mean),diag(as.numeric(var)))
    dim(sample) = dim(mean)
  } else{
    sample = mean
  }
  return(list(lagp_fit=lagp_fit,mean=mean,var=var,sample=sample))
}

# FUNCTION: get_basis
{### Accepts a matrix of response values Y, and returns the basis decomposition from SVD
## Parameters:
# Y (numeric)       : matrix of response values to decompose into basis representation. Y should be
#                     formatted such that each row represents the multivariate response for a single experiment
# n.pc (int)        : number of principal components desired from basis decomposition
# pct.var (numeric) : desired percentage of variability explained by PC's. Overrides n.pc if specified
## Returns: List containings the following
# B (numeric)       : matrix of n.pc basis vectors for Ytrans
# V.t (numeric)      : matrix of n.pc response vectors for Ytrans
# n.pc (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating sim basis
# mv_calib.bias for generating discrepancy basis
}
get_basis = function(Y, n.pc = 1, pct.var = NULL, full.basis=F, B=NULL){

  return.list = list()
  
  # Compute SVD
  nExp = ncol(Y)
  if(is.null(B)){
    svdY = svd(Y)
    return.list$B = svdY$u %*% diag(svdY$d) / sqrt(nExp)
    return.list$V.t = t(svdY$v) * sqrt(nExp)
  
    # Return only the required number of basis vectors
    percent.variance = zapsmall((svdY$d) / sum(svdY$d))
    if(!is.null(pct.var)){
      # percentage of variance explained specified
      n.pc = which.max(cumsum(percent.variance)>=pct.var)
    }
    
    return.list$n.pc = n.pc
    if(full.basis){
      return.list$pct.var = percent.variance
      return.list$B = as.matrix(return.list$B)
      return.list$V.t = as.matrix(return.list$V.t)
    } else{
      return.list$pct.var = percent.variance[1:return.list$n.pc]
      return.list$B = as.matrix(return.list$B[,1:return.list$n.pc,drop=F])
      return.list$V.t = as.matrix(return.list$V.t[1:return.list$n.pc,,drop=F])
    }
  } else{
    # manual basis matrix was passed
    return.list$B = B
    return.list$n.pc = ncol(B)
    return.list$V.t = solve(t(B)%*%B)%*%t(B)%*%Y
  }
  return(return.list)
}

# FUNCTION: get_obs_basis
{### Generated obs basis vectors by interpolating sim basis vectors
## Parameters:
# sim.basis       : output from get_basis
# Yobs        : matrix (n.y x n) field data
# YindSim : matrix of field data indices of observation
## Returns: List containings the following
# B (numeric)       : matrix of n.pc basis vectors
# V.t (numeric)      : matrix of n.pc response vectors
# n.pc (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating obs basis
}
get_obs_basis = function(sim.basis,Yobs,YindSim,YindObs,sigY=NULL){
  obs.basis = list()
  n.pc = sim.basis$n.pc
  obs.basis$n.pc = n.pc
  if(!isTRUE(all.equal(YindObs,YindSim))){
    B = matrix(nrow=nrow(YindObs),ncol=obs.basis$n.pc)
    if(ncol(YindSim)==1){
      # 1d interpolation
      for(j in 1:n.pc){
        B[,j] = interp1(as.numeric(YindSim),as.numeric(sim.basis$B[,j]),as.numeric(YindObs))
      }
    } else if(ncol(YindSim)==2){
      # 2d interpolation
      y = unique(YindSim[,1]); x = unique(YindSim[,2])
      yp = YindObs[,1]; xp = YindObs[,2]
      for(i in 1:n.pc){
        B[,i] = interp2(x,y,matrix(sim.basis$B[,i],length(y),length(x),byrow = F),xp,yp)
      }
    }
    obs.basis$B = B
  } else{
    cat('setting Bobs==Bsim')
    obs.basis$B = sim.basis$B
  }
  # compute V.t
  # Bridge = 1e-6 * diag(rep(1,n.pc)) # in case Bprod is ill-conditioned
  # if(!is.null(sigY)){
  #   lamY = (1/sigY)*eye(dim(obs.basis$B)[1])
  # } else{
  #   lamY = eye(dim(obs.basis$B)[1])
  # }
  # Bprod = t(obs.basis$B)%*%lamY%*%obs.basis$B
  # obs.basis$V.t = solve(Bprod + Bridge) %*% t(obs.basis$B)%*%lamY%*%Yobs
  obs.basis$V.t = solve(t(obs.basis$B)%*%obs.basis$B)%*%t(obs.basis$B)%*%Yobs
  return(obs.basis)
}

# FUNCTION: mv_lengthscales
# TO DO: this function needs to be updated for only using a subset of the data
### Computes estimated length scale parameters for use in stretching and compressing the inputs space
# Parameters:
# X (numeric) : untransformed matrix of inputs
# V.t (numeric) : matrix 
# ls_prior_ab (numeric) : 2 vector of parameters for gamma prior on length-scales
mv_lengthscales = function(XT.data,V.t,g=1e-7,ls_prior_ab=NULL,subsample=F,subsample.size=100,seed=NULL){

  XT = cbind(XT.data$sim$X$trans,XT.data$sim$T$trans)
  if(subsample){
    if(!is.null(seed)){
      set.seed(seed)
    }
    samp.id = sample(1:nrow(XT),subsample.size)
    XT = XT[samp.id,]
    V.t = V.t[,samp.id]
  }
  p.x = XT.data$p.x
  p.t = XT.data$p.t
  n.pc = nrow(V.t)
  
  dConfig = darg(list(mle = TRUE, max = 100), XT)
  if(is.null(ls_prior_ab)){
    ls_prior_ab = dConfig$ab
  }
  GPlist = vector(mode='list',length=n.pc)
  for(i in 1:n.pc){
    GPlist[[i]] = newGPsep(X = XT,
                           Z = V.t[i, ],
                           d = rep(dConfig$start, ncol(XT)),
                           g = g,
                           dK = TRUE)
  }
  
  estLenscales = lapply(1:n.pc, function(i) mleGPsep(GPlist[[i]],param='d',
                                                        tmin = dConfig$min, 
                                                        tmax = dConfig$max,
                                                        ab=ls_prior_ab,
                                                        maxit=200)$d)
  for(i in n.pc){
    deleteGPsep(GPlist[[i]])
  }
  
  return = list()
  return$XT = estLenscales
  if(p.x>0){
    return$X = lapply(1:n.pc, function(i) estLenscales[[i]][1:p.x])
  } else{
    return$X = list()
  }
  if(p.t>0){
    return$T = lapply(1:n.pc, function(i) estLenscales[[i]][(p.x+1):(p.x+p.t)])
  } else{
    return$T = list()
  }
  return(return)
}

# FUNCTION: unit_xform
### Transform columns to lie in unit interval
# Parameters
# X (numeric) : matrix to be transformed
# Returns
# Xtrans (numeric) : transformed matrix
unit_xform = function(X,X.min=NULL,X.range=NULL){
  X = as.matrix(X)
  if(is.null(X.min)){
    X.min = apply(X,2,min)
  }
  if(is.null(X.range)){
    Xmax = apply(X,2,max)
    X.range = Xmax - X.min
  }
  if(!isTRUE(all.equal(X.range,rep(0,length(X.range))))){
    Xtrans = t( (t(X) - X.min) / X.range )
  } else{
    # just return X if dummyx
    Xtrans = X
  }
  return(list(orig=X,
              trans=Xtrans,
              min=X.min,
              range=X.range))
}

transform_y = function(Ysim,YindSim=NULL,Yobs=NULL,YindObs=NULL,center=T,scale=T, scaletype='rowwise'){
  sim.list = list()
  obs.list = list()
  
  dim = ncol(YindSim)
  
  # scale Ysim
  sim.list$orig = Ysim
  if(center){
    sim.list$mean = rowMeans(Ysim)
  } else{
    sim.list$mean = rep(0,nrow(Ysim))
  }
  if(scale){
    if(scaletype=='rowwise'){
      sim.list$sd = apply(Ysim,1,sd)
      # dont scale places with var 0
      sim.list$sd[sim.list$sd==0]=1
    } else if(scaletype=='scalar') {
      sd = sd(Ysim)
      sim.list$sd = rep(sd,nrow(Ysim))
    }
  } else{
    sim.list$sd = rep(1,nrow(Ysim))
  }
  sim.list$trans = t(scale(t(Ysim),sim.list$mean,sim.list$sd))
  
  # scale Yobs
  if(!is.null(Yobs)){
    obs.list$orig = Yobs
    obs.list$n.y = dim(Yobs)[1]
    if(!isTRUE(all.equal(YindObs,YindSim))){
      # interpolation is needed
      if(dim==1){
        # 1d interpolation of sim mean onto yobs support
        obs.list$mean = interp1(as.numeric(YindSim),sim.list$mean,as.numeric(YindObs))
        obs.list$sd = interp1(as.numeric(YindSim),sim.list$sd,as.numeric(YindObs))
        obs.list$trans = t(scale(t(Yobs),obs.list$mean,obs.list$sd))
      } else if(dim==2){
        # 2d interpolation of sim mean onto yobs support
        x = unique(YindSim[,1]); y = unique(YindSim[,2])
        xp = YindObs[,1]; yp = YindObs[,2]
        obs.list$mean = interp2(x,y,matrix(sim.list$mean,length(y),length(x),byrow = T),xp,yp)
        if(scale){
          if(scaletype=='rowwise'){
            obs.list$sd = interp2(x,y,matrix(sim.list$sd,length(y),length(x),byrow = T),xp,yp)
          } else if(scaletype=='scalar'){
            obs.list$sd = rep(sd,nrow(Yobs))
          }
        } else{
          obs.list$sd = rep(1,nrow(Yobs))
        }
        obs.list$trans = t(scale(t(Yobs),obs.list$mean,obs.list$sd))
      }
    } else{
      obs.list$trans = t(scale(t(Yobs),sim.list$mean,sim.list$sd)) 
      obs.list$mean = sim.list$mean
      obs.list$sd = sim.list$sd
    }
  }
  return(list(sim=sim.list,obs=obs.list,
              n=ncol(Yobs),m=ncol(Ysim),n.y=nrow(Ysim)))
}

# FUNCTION: transform_xt
### Transforms inputs to lie on the unit hypercube
# T.obs are optional known calibration parameters for the observed data
# transforming these is useful for checking predictions at obs locations
transform_xt = function(X.sim=NULL,T.sim=NULL,X.obs=NULL,T.obs=NULL){
  sim.list = list()
  obs.list = list()
  if(is.null(X.sim) & is.null(T.sim)){
    cat('One of X.sim or T.sim must be specified')
    break
  }
  if(!is.null(X.sim)){
    p.x = ncol(X.sim)
    sim.list$X = unit_xform(X.sim)
  } else{
    p.x = 0
  }
  if(!is.null(T.sim)){
    p.t = ncol(T.sim)
    sim.list$T = unit_xform(T.sim)
  } else{
    p.t = 0
  }

  if(!is.null(X.obs)){
    obs.list$X = unit_xform(X.obs,X.min=sim.list$X$min,X.range=sim.list$X$range)
  }
  if(!is.null(T.obs)){
    obs.list$T = unit_xform(T.obs,X.min=sim.list$T$min,X.range=sim.list$T$range)
  }
  return(list(sim=sim.list,obs=obs.list,p.x=p.x,p.t=p.t))
}

# FUNCTOIN: mean_sd_xform
### Centers and/or scales a matrix
mean_sd_xform = function(X,center=T,scale=T){
  tmp = t(scale(t(X),center,scale))
  if(scale){
    sd = attr(tmp,'scaled:scale')
  }else{
    sd = rep(1,nrow(X))
  }
  return(list(orig=X,
              trans=as.matrix(tmp),
              mean=attr(tmp,'scaled:center'),
              sd=sd))
}

get_SC_inputs = function(est.ls,XT.data,n.pc){
  return = list()
  if(length(est.ls$X)>0){
    return$X.sim = lapply(1:n.pc, function(i) sc_inputs(XT.data$sim$X$trans,est.ls$X[[i]]))
  }
  if(length(est.ls$T)>0){
    return$T.sim =  lapply(1:n.pc, function(i) sc_inputs(XT.data$sim$T$trans,est.ls$T[[i]]))
  }
  return$XT.sim = lapply(1:n.pc, function(i) cbind(return$X.sim[[i]],return$T.sim[[i]]))
  
  if(!is.null(XT.data$obs$X)){
    return$X.obs = lapply(1:n.pc, function(i) sc_inputs(XT.data$obs$X$trans,est.ls$X[[i]]))
  }
  if(!is.null(XT.data$obs$T)){
    return$T.obs = lapply(1:n.pc, function(i) sc_inputs(XT.data$obs$T$trans,est.ls$T[[i]]))
  }
  return$p.x = XT.data$p.x
  return$p.t = XT.data$p.t
  return(return)
}
# FUNCTION: sc_inputs
### Function to apply stretching and compressing to inputs
# Parameters
# X (numeric) : matrix of inputs to be stretched and compressed by lengthscale parameters
# ls: vector of lengthscale parameters length(ls)==ncol(X)
# Returns
# Xtrans (numeric) : transformed matrix
# sc_inputs = function(X,ls){
#   return(t(t(X) / sqrt(ls))) # t(t(X)) is for the case where X is a vector (1d input space)
# }
sc_inputs = function(X,ls){
  return(sweep(X, 2, sqrt(ls), FUN = '/'))
}

neg_ll = function(theta,mvc.data,lite=T,sample=T,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  n.pc = mvc.data$sim.basis$n.pc
  n = mvc.data$Y.data$n
  n.y = mvc.data$Y.data$n.y
  p.t = mvc.data$XT.data$p.t
  p.x = mvc.data$XT.data$p.x
  #theta = matrix(theta,nrow=1)
  bias = mvc.data$bias
  
  w = NULL; v = NULL; delta = NULL
  ll = numeric(n)
  
  # transform theta by estimated length-scales
  #theta.sc = lapply(1:n.pc, function(jj) sc_inputs(theta,mvc.data$est.ls$T[[jj]]))
  # theta.sc = lapply(1:n.pc, function(jj) theta/mvc.data$est.ls$T[[jj]])
  # # make prediction matrix from X.obs and theta
  # if(p.t>1){
  #   X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(replicate(n,theta.sc[[jj]],simplify='matrix'))))
  # } else{
  #   X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(t(replicate(n,theta.sc[[jj]],simplify='matrix')))))
  # }
  # 
  theta.sc = lapply(1:n.pc, function(jj) matrix(theta/mvc.data$est.ls$T[[jj]],nrow=n,ncol=p.t,byrow = T))
  # make prediction matrix from X.obs and theta
  if(p.t>1){
    # XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(replicate(n,theta.sc[[jj]],simplify='matrix'))))
    if(p.x>0){
      XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],theta.sc[[jj]]))
    } else{
      XT.sc = theta.sc
    }
  } else{
    XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(t(replicate(n,theta.sc[[jj]],simplify='matrix')))))
  }
  
  # Predict from emulator at [X.obs,theta]
  w = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
                     Z=lapply(1:n.pc,function(jj) mvc.data$sim.basis$V.t[jj,]),
                     XX=XT.sc,
                     start=start,
                     end=end)
  
  if(sample){
    w$sample = rmvnorm(1,as.numeric(w$mean),diag(as.numeric(w$var)))
    dim(w$sample) = dim(w$mean)
  } else{
    w$sample = w$mean
  }

  # residuals
  y.resid = mvc.data$Y.data$obs$trans - w_to_y(w$sample,mvc.data$obs.basis$B)
  ssq.hat = mean(y.resid^2)
  if(!lite){
    y.resid.orig = mvc.data$Y.data$obs$orig - w_to_y(w$sample,mvc.data$obs.basis$B,mvc.data$Y.data$obs$mean,mvc.data$Y.data$obs$sd)
    ssq.hat.orig = mean(y.resid.orig^2)
  }
  if(bias){
    v$basis = get_basis(y.resid,pct.var=.95,B=mvc.data$obs.basis$D)
    if(laGP.bias){ # use laGP for a faster bias model
      v = append(v,aGPsep_SC_mv(X=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                       Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
                       XX=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                       start=start.bias,
                       end=end.bias,
                       bias=T))
    } else{       # use full GP for the bias model
      v$GPs = vector(mode='list',length=v$basis$n.pc)

      for(k in 1:v$basis$n.pc){
        d = darg(d=list(mle=T),X=mvc.data$XT.data$obs$X$trans)
        g = garg(g=list(mle=T),y=v$basis$V.t[k,])
        # Our responses have changed wrt to inputs, so we should not use d=1 an.ymore
        v$GPs[[k]] <- newGPsep(X=mvc.data$XT.data$obs$X$trans,
                               Z=v$basis$V.t[k,],
                               d=d$start, g=g$start, dK=TRUE)
        cmle <- jmleGPsep(v$GPs[[k]],
                          drange=c(d$min, d$max),
                          grange=c(g$min, g$max),
                          dab=d$ab,
                          gab=g$ab)
        pred = predGPsep(v$GPs[[k]], XX = mvc.data$XT.data$obs$X$trans, lite=T)
        v$mean = rbind(v$mean, pred$mean)
        v$var = rbind(v$var, pred$s2)
        if(lite){
          deleteGPsep(v$GPs[[k]])
        }
      }
    }
    # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
    if(sample){
      v$sample = rmvnorm(1,as.numeric(v$mean),diag(as.numeric(v$var)))
      dim(v$sample) = dim(v$mean)
    } else{
      v$sample = v$mean
    }
    y.resid = y.resid - w_to_y(v$sample,v$basis$B)
    ssq.hat = mean(y.resid^2)
    if(!lite){
      new.resid.orig = y.resid.orig - w_to_y(v$sample,v$basis$B,sd=mvc.data$Y.data$obs$sd)
      ssq.hat.orig = mean(new.resid.orig^2)
    }
    n.y = nrow(y.resid)
    BD = cbind(mvc.data$obs.basis$B,v$basis$B)
  } else{
    BD = mvc.data$obs.basis$B
  }
  
  if(sample){
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of y.resid
    ldetS = n.y*log(ssq.hat)
    Sinv = (1/ssq.hat)*eye(n.y)
    ll = apply(y.resid,2,function(d) -.5*(ldetS + d%*%Sinv%*%d))
  } else{
    r = 1e-8
    lambda = 1/ssq.hat
    tBD = t(BD)
    rankBD = rankMatrix(BD)
    BDtBD = tBD%*%BD
    if(rankBD != ncol(BD)){
      BDtBDinv = solve(BDtBD + r*eye(nrow(BDtBD)))
    } else{
      BDtBDinv = solve(BDtBD) # need to get rid of solves to make code faster
    }
    ldetBDtBD = determinant(BDtBD)$modulus
    I.n.y = eye(n.y)
    B.hat = BDtBDinv%*%tBD%*%y.resid
    
    SigB.hat = lapply(1:n, function(i) as.matrix(diag(c(w$var[,i],v$var[,i]),nrow=length(c(w$var[,i],v$var[,i]))) + ssq.hat*BDtBDinv))
    l_beta_hat = sapply(1:n, function(i) dmvnorm(x=as.numeric(B.hat[,i]),sigma = SigB.hat[[i]], log = T))
    ll = sapply(1:n, function(i) -.5*((n.y-rankBD)*log(2*pi*ssq.hat) + ldetBDtBD +
                                                  lambda*t(y.resid[,i])%*%(I.n.y-BD%*%BDtBDinv%*%tBD)%*%y.resid[,i]) + 
                                           l_beta_hat[i])
  }
  if(lite){
    return(-sum(ll) - log(ssq.hat) + sum(dbeta(theta,2,2,log=T)))
  } else{
    return = list()
    return$ll = -sum(ll) - log(ssq.hat) + sum(dbeta(theta,2,2,log=T))
    return$ssq.hat = ssq.hat
    return$ssq.hat.orig = ssq.hat.orig
    return$w = w
    return$v = v
    return(return)
  }
}

energy_score_snomadr = function(theta,mvc.data,sample=T,n.samples=100,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  n.pc = nrow(mvc.data$obs.basis$V.t)
  n = nrow(mvc.data$SC.inputs$X.obs[[1]])
  n.y = nrow(mvc.data$Y.data$obs$orig)
  theta = matrix(theta,nrow=1)
  bias = mvc.data$bias
  ym = mvc.data$Y.data$obs$mean
  ysd = mvc.data$Y.data$obs$sd
  w = NULL; v = NULL
  eta = NULL; delta = NULL
  
  # transform theta by estimated length-scales
  theta.sc = lapply(1:n.pc, function(jj) sc_inputs(theta,mvc.data$est.ls$T[[jj]]))
  # make prediction matrix from X.obs and theta
  if(mvc.data$SC.inputs$p.t>1){
    X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(replicate(n,theta.sc[[jj]],simplify='matrix'))))
  } else{
    X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(t(replicate(n,theta.sc[[jj]],simplify='matrix')))))
  }
  
  # Predict from emulator at [X.obs,theta]
  w = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
                   Z=lapply(1:n.pc,function(jj) mvc.data$sim.basis$V.t[jj,]),
                   XX=X,
                   start=start,
                   end=end)
  
  if(sample){
    w$sample = rmvnorm(1,as.numeric(w$mean),diag(as.numeric(w$var)))
    dim(w$sample) = dim(w$mean)
  } else{
    w$sample = w$mean
  }
  
  # residuals
  y.resid = mvc.data$Y.data$obs$trans - w_to_y(w$sample,mvc.data$obs.basis$B)
  ssq.hat = mean(y.resid^2)
  if(bias){
    v$basis = get_basis(y.resid,pct.var=.95,B=mvc.data$obs.basis$D)
    if(laGP.bias){ # use laGP for a faster bias model
      v = append(v,aGPsep_SC_mv(X=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                                Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
                                XX=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                                start=start.bias,
                                end=end.bias,
                                bias=T))
    } else{       # use full GP for the bias model
      v$GPs = vector(mode='list',length=v$basis$n.pc)
      
      for(k in 1:v$basis$n.pc){
        d = darg(d=list(mle=T),X=mvc.data$XT.data$obs$X$trans)
        g = garg(g=list(mle=T),y=v$basis$V.t[k,])
        # Our responses have changed wrt to inputs, so we should not use d=1 an.ymore
        v$GPs[[k]] <- newGPsep(X=mvc.data$XT.data$obs$X$trans,
                               Z=v$basis$V.t[k,],
                               d=d$start, g=g$start, dK=TRUE)
        cmle <- jmleGPsep(v$GPs[[k]],
                          drange=c(d$min, d$max),
                          grange=c(g$min, g$max),
                          dab=d$ab,
                          gab=g$ab)
        pred = predGPsep(v$GPs[[k]], XX = mvc.data$XT.data$obs$X$trans, lite=T)
        v$mean = rbind(v$mean, pred$mean)
        v$var = rbind(v$var, pred$s2)
        deleteGPsep(v$GPs[[k]])
      }
    }
    # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
    if(sample){
      v$sample = rmvnorm(1,as.numeric(v$mean),diag(as.numeric(v$var)))
      dim(v$sample) = dim(v$mean)
    } else{
      v$sample = v$mean
    }
    y.resid = y.resid - w_to_y(v$sample,v$basis$B)
    ssq.hat = mean(y.resid^2)
  }
  
  y.samp = array(0,dim=c(n.samples,n.y,n))
  eta = mvc.data$obs.basis$B%*%w$sample
  if(bias){
    delta = mvc.data$obs.basis$D%*%drop(v$sample)
  } else{
    delta = 0*eta
  }
  
  if(sample){
    B = mvc.data$obs.basis$B
    if(bias){
      D = mvc.data$obs.basis$D
      for(j in 1:n){
        y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j]+delta[,j],sigma=ssq.hat*eye(n.y) + 
                                B%*%diag(w$var[,j])%*%t(B) + D%*%diag(v$var[,j])%*%t(D)) * ysd + ym
      }
    } else{
      for(j in 1:n){
        y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j],sigma=ssq.hat*eye(n.y) + 
                                B%*%diag(w$var[,j])%*%t(B)) * ysd + ym
      }
    }
  } else{
    for(j in 1:n){
      y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j]+delta[,j],sigma=ssq.hat*eye(n.y)) * ysd + ym
    }
  }

  # compute energy score
  ES = mean(energy_score(y.samp,mvc.data$Y.data$obs$orig))
  return(ES)
}
interval_score_snomadr = function(theta,mvc.data,sample=T,n.samples=100,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50,
                                  lite=T){
  n.pc = nrow(mvc.data$obs.basis$V.t)
  n = nrow(mvc.data$SC.inputs$X.obs[[1]])
  n.y = nrow(mvc.data$Y.data$obs$orig)
  theta = matrix(theta,nrow=1)
  bias = mvc.data$bias
  ym = mvc.data$Y.data$obs$mean
  ysd = mvc.data$Y.data$obs$sd
  w = NULL; v = NULL
  eta = NULL; delta = NULL
  
  # transform theta by estimated length-scales
  theta.sc = lapply(1:n.pc, function(jj) sc_inputs(theta,mvc.data$est.ls$T[[jj]]))
  # make prediction matrix from X.obs and theta
  if(mvc.data$SC.inputs$p.t>1){
    X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(replicate(n,theta.sc[[jj]],simplify='matrix'))))
  } else{
    X = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(t(replicate(n,theta.sc[[jj]],simplify='matrix')))))
  }
  
  # Predict from emulator at [X.obs,theta]
  w = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
                   Z=lapply(1:n.pc,function(jj) mvc.data$sim.basis$V.t[jj,]),
                   XX=X,
                   start=start,
                   end=end)
  
  if(sample){
    w$sample = rmvnorm(1,as.numeric(w$mean),diag(as.numeric(w$var)))
    dim(w$sample) = dim(w$mean)
  } else{
    w$sample = w$mean
  }
  
  # residuals
  y.resid = mvc.data$Y.data$obs$trans - w_to_y(w$sample,mvc.data$obs.basis$B)
  ssq.hat = mean(y.resid^2)
  if(bias){
    v$basis = get_basis(y.resid,pct.var=.95,B=mvc.data$obs.basis$D)
    if(laGP.bias){ # use laGP for a faster bias model
      v = append(v,aGPsep_SC_mv(X=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                                Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
                                XX=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
                                start=start.bias,
                                end=end.bias,
                                bias=T))
    } else{       # use full GP for the bias model
      v$GPs = vector(mode='list',length=v$basis$n.pc)
      
      for(k in 1:v$basis$n.pc){
        d = darg(d=list(mle=T),X=mvc.data$XT.data$obs$X$trans)
        g = garg(g=list(mle=T),y=v$basis$V.t[k,])
        # Our responses have changed wrt to inputs, so we should not use d=1 an.ymore
        v$GPs[[k]] <- newGPsep(X=mvc.data$XT.data$obs$X$trans,
                               Z=v$basis$V.t[k,],
                               d=d$start, g=g$start, dK=TRUE)
        cmle <- jmleGPsep(v$GPs[[k]],
                          drange=c(d$min, d$max),
                          grange=c(g$min, g$max),
                          dab=d$ab,
                          gab=g$ab)
        pred = predGPsep(v$GPs[[k]], XX = mvc.data$XT.data$obs$X$trans, lite=T)
        v$mean = rbind(v$mean, pred$mean)
        v$var = rbind(v$var, pred$s2)
        if(lite){
          deleteGPsep(v$GPs[[k]])
        }
      }
    }
    # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
    if(sample){
      v$sample = rmvnorm(1,as.numeric(v$mean),diag(as.numeric(v$var)))
      dim(v$sample) = dim(v$mean)
    } else{
      v$sample = v$mean
    }
    y.resid = y.resid - w_to_y(v$sample,v$basis$B)
    ssq.hat = mean(y.resid^2)
  }
  
  # get samples of y ~ N(Bw+Dv,ssq*I)
  y.samp = array(0,dim=c(n.samples,n.y,n))
  eta = mvc.data$obs.basis$B%*%w$sample
  if(bias){
    delta = mvc.data$obs.basis$D%*%drop(v$sample)
  } else{
    delta = 0*eta
  }
  
  if(sample){
    B = mvc.data$obs.basis$B
    if(bias){
      D = mvc.data$obs.basis$D
      for(j in 1:n){
        y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j]+delta[,j],sigma=ssq.hat*eye(n.y) + 
                                B%*%diag(w$var[,j])%*%t(B) + D%*%diag(v$var[,j])%*%t(D)) * ysd + ym
      }
    } else{
      for(j in 1:n){
        y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j],sigma=ssq.hat*eye(n.y) + 
                                B%*%diag(w$var[,j])%*%t(B)) * ysd + ym
      }
    }
  } else{
    for(j in 1:n){
      y.samp[,,j] = rmvnorm(n.samples,mean=eta[,j]+delta[,j],sigma=ssq.hat*eye(n.y)) * ysd + ym
    }
  }
  
  # return interval score for samples
  IS = interval_score(y.samp,mvc.data$Y.data$obs$orig)
  if(lite){
    return(IS)
  } else{
    return(list(is=IS,ssq.hat=ssq.hat,w=w,v=v))
  }
}

calib_init = function(t.init,mvc.data,sample=F,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  p.t = mvc.data$XT.data$p.t
  
  ll = sapply(1:nrow(t.init), function(j) -neg_ll(t.init[j,,drop=F],mvc.data,lite=T,sample=sample,start=start,end=end,
                                                laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias))
  
  ll_df = as.data.frame(matrix(cbind(ll,as.matrix(t.init)),nrow=length(ll),ncol=p.t+1))
  names = c("ll")
  names = c(names,paste0("theta",1:p.t))
  colnames(ll_df) = names
  ll_df = ll_df[order(ll_df$ll,decreasing=T),]
  
  return(ll_df)
}
calib_init_es = function(t.init,mvc.data,sample=F,n.samples=100,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  p.t = mvc.data$XT.data$p.t

  es = sapply(1:nrow(t.init), function(j) energy_score_snomadr(t.init[j,,drop=F],mvc.data,sample=sample,
                                                               n.samples=n.samples,start=start,end=end,
                                                  laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias))
  
  es_df = as.data.frame(matrix(cbind(es,as.matrix(t.init)),nrow=length(es),ncol=p.t+1))
  names = c("es")
  names = c(names,paste0("theta",1:p.t))
  colnames(es_df) = names
  es_df = es_df[order(es_df$es,decreasing=T),]
  
  return(es_df)
}
calib_init_is = function(t.init,mvc.data,sample=F,n.samples=100,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  p.t = mvc.data$XT.data$p.t
  
  is = sapply(1:nrow(t.init), function(j) interval_score_snomadr(t.init[j,,drop=F],mvc.data,sample=sample,
                                                               n.samples=n.samples,start=start,end=end,
                                                               laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias))
  
  is_df = as.data.frame(matrix(cbind(is,as.matrix(t.init)),nrow=length(is),ncol=p.t+1))
  names = c("is")
  names = c(names,paste0("theta",1:p.t))
  colnames(is_df) = names
  is_df = is_df[order(is_df$is,decreasing=T),]
  
  return(is_df)
}

mv_calib = function(mvc.data,init,nrestarts=1,sample=T,optimize=T,
                    opts=list("MAX_BB_EVAL" = 1000, "INITIAL_MESH_SIZE" = .1,
                              "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0),
                    start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){

  if(optimize){
    # nrestarts cannot be greater than our set of initial values
    if(nrestarts>nrow(init)){
      cat('nrestarts>nrow(init), setting nrestarts=nrow(init)')
      nrestarts=nrow(init)
    }
    out = NULL
    for(i in 1:nrestarts){
      outi = snomadr(neg_ll, n = mvc.data$SC.inputs$p.t, bbin = rep(0, mvc.data$SC.inputs$p.t), bbout = 0, 
                     x0 = as.numeric(init[i,2:ncol(init)]),
                     lb = rep(.01, mvc.data$SC.inputs$p.t), 
                     ub = rep(.99, mvc.data$SC.inputs$p.t),
                     opts = opts,
                     print.output=F,
                     mvc.data = mvc.data,
                     lite=T,
                     sample=sample,
                     start=start,
                     end=end,
                     laGP.bias=laGP.bias,
                     start.bias=start.bias,
                     end.bias=end.bias)
      if(is.null(out) || outi$objective < out$objective) out <- outi
    }
    return = neg_ll(matrix(out$solution,nrow=1),mvc.data,lite=F,sample=sample,laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias)
    return$snomadr = out
    return$theta.hat = out$solution
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = out$solution * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  } else{
    return = neg_ll(init,mvc.data,lite=F,sample=sample,start=start,end=end,laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias)
    return$theta.hat = init
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = init * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  }
  return(return)
}

mv_calib_es = function(mvc.data,init,nrestarts=1,sample=T,optimize=T,n.samples=100,
                    opts=list("MAX_BB_EVAL" = 1000, "INITIAL_MESH_SIZE" = .1,
                              "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0),
                    start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  
  if(optimize){
    # nrestarts cannot be greater than our set of initial values
    if(nrestarts>nrow(init)){
      cat('nrestarts>nrow(init), setting nrestarts=nrow(init)')
      nrestarts=nrow(init)
    }
    out = NULL
    for(i in 1:nrestarts){
      outi = snomadr(energy_score_snomadr, n = mvc.data$SC.inputs$p.t, bbin = rep(0, mvc.data$SC.inputs$p.t), bbout = 0, 
                     x0 = as.numeric(init[i,2:ncol(init)]),
                     lb = rep(.01, mvc.data$SC.inputs$p.t), 
                     ub = rep(.99, mvc.data$SC.inputs$p.t),
                     opts = opts,
                     print.output=F,
                     mvc.data = mvc.data,
                     sample=sample,
                     n.samples=n.samples,
                     start=start,
                     end=end,
                     laGP.bias=laGP.bias,
                     start.bias=start.bias,
                     end.bias=end.bias)
      if(is.null(out) || outi$objective < out$objective) out <- outi
    }
    return = list()
    return$snomadr = out
    return$es = out$objective
    return$theta.hat = out$solution
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = out$solution * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  } else{
    return$es = energy_score_snomadr(init,mvc.data,sample=sample,n.samples=n.samples,
                                  start=start,end=end,laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias)
    return$theta.hat = init
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = init * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  }
  return(return)
}

mv_calib_is = function(mvc.data,init,nrestarts=1,sample=T,optimize=T,n.samples=100,
                       opts=list("MAX_BB_EVAL" = 1000, "INITIAL_MESH_SIZE" = .1,
                                 "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0),
                       start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
  
  if(optimize){
    # nrestarts cannot be greater than our set of initial values
    if(nrestarts>nrow(init)){
      cat('nrestarts>nrow(init), setting nrestarts=nrow(init)')
      nrestarts=nrow(init)
    }
    out = NULL
    for(i in 1:nrestarts){
      outi = snomadr(interval_score_snomadr, n = mvc.data$SC.inputs$p.t, bbin = rep(0, mvc.data$SC.inputs$p.t), bbout = 0, 
                     x0 = as.numeric(init[i,2:ncol(init)]),
                     lb = rep(.01, mvc.data$SC.inputs$p.t), 
                     ub = rep(.99, mvc.data$SC.inputs$p.t),
                     opts = opts,
                     print.output=F,
                     mvc.data = mvc.data,
                     sample=sample,
                     n.samples=n.samples,
                     start=start,
                     end=end,
                     laGP.bias=laGP.bias,
                     start.bias=start.bias,
                     end.bias=end.bias,
                     lite=T)
      if(is.null(out) || outi$objective < out$objective) out <- outi
    }
    return = interval_score_snomadr(out$solution,mvc.data,sample,n.samples,start,end,laGP.bias,start.bias,end.bias,lite=F)
    return$snomadr = out
    return$is = out$objective
    return$theta.hat = out$solution
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = out$solution * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  } else{
    return = interval_score_snomadr(init,mvc.data,sample=sample,n.samples=n.samples,
                                     start=start,end=end,laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias,lite=F)
    return$theta.hat = init
    # TO DO: will this work for multivariate theta?
    return$theta.hat.orig = init * mvc.data$XT.data$sim$T$range + mvc.data$XT.data$sim$T$min
  }
  return(return)
}

# mcmc = function(mvc.data,t.init,n.samples=100,n.burn=0,
#                 prop.step=rep(.25,mvc.data$XT.data$p.t),
#                 verbose=F,start=6,end=50,
#                 laGP.bias=F,start.bias=6,end.bias=50,
#                 X.pred=NULL,n.pred.samples=0){
#   ll_mcmc = function(theta,predict=F){
# 
#     theta = matrix(theta,nrow=1)
#     ll = numeric(n)
#     w = NULL; v = NULL;
#     if(!is.null(X.pred) & predict){
#       w.pred = NULL
#       v.pred = NULL
#     }
#     
#     # transform theta by estimated length-scales
#     theta.sc = lapply(1:n.pc, function(j) sc_inputs(theta,mvc.data$est.ls$T[[j]]))
#     # make prediction matrix from X.obs and theta
#     X = lapply(1:n.pc, function(j) cbind(mvc.data$SC.inputs$X.obs[[j]],matrix(rep(theta.sc[[j]],n),ncol=p.t,byrow=T)))
#     # Predict from emulator at [X.obs,theta]
#     w = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
#                        Z=lapply(1:n.pc,function(j) mvc.data$sim.basis$V.t[j,]),
#                        XX=X,
#                        start=start,
#                        end=end)
#     w$sample = rmvnorm(1,as.numeric(w$mean),diag(as.numeric(w$var)))
#     dim(w$sample) = dim(w$mean)
#     
#     # If a new X was passed in for prediction
#     if(!is.null(X.pred) & predict){
#       # cbind X.pred.sc to theta.sc
#       X = lapply(1:n.pc, function(j) cbind(X.pred.sc$X.obs[[j]],matrix(rep(theta.sc[[j]],nrow(X.pred.sc$X.obs[[j]])),ncol=p.t,byrow = T)))
#       w.pred = aGPsep_SC_mv(X=X.pred.sc$XT.sim,
#                          Z=lapply(1:n.pc,function(i) mvc.data$sim.basis$V.t[i,]),
#                          XX=X,
#                          start=start,
#                          end=end)
#       w.pred$sample = rmvnorm(1,as.numeric(w.pred$mean),diag(as.numeric(w.pred$var)))
#       dim(w.pred$sample) = dim(w.pred$mean)
#     }
#     
#     # residuals
#     y.resid = mvc.data$Y.data$obs$trans - w_to_y(w$sample,mvc.data$obs.basis$B)
#     
#     if(mvc.data$bias){
#       v$basis = get_basis(y.resid,pct.var=.95,B=mvc.data$obs.basis$D)
#       if(laGP.bias){ # use laGP for a faster bias model
#         v = append(v,aGPsep_SC_mv(X=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
#                                   Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
#                                   XX=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
#                                   start=start.bias,
#                                   end=end.bias,
#                                   bias=T))
#         if(!is.null(X.pred) & predict){
#           v.pred = aGPsep_SC_mv(X=rep(list(mvc.data$XT.data$obs$X$trans),v$basis$n.pc),
#                          Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
#                          XX=rep(list(X.pred.std$obs$X$trans),v$basis$n.pc),
#                          start=start.bias,
#                          end=end.bias,
#                          bias=T)
#         }
#       } else{
#         v$GPs = vector(mode='list',length=v$basis$n.pc)
#         
#         for(k in 1:v$basis$n.pc){
#           d = darg(d=list(mle=T),X=mvc.data$XT.data$obs$X$trans)
#           g = garg(g=list(mle=T),y=v$basis$V.t[k,])
#           # Our responses have changed wrt to inputs, so we should not use d=1 an.ymore
#           v$GPs[[k]] <- newGPsep(X=mvc.data$XT.data$obs$X$trans,
#                                      Z=v$basis$V.t[k,],
#                                      d=d$start, g=g$start, dK=TRUE)
#           cmle <- jmleGPsep(v$GPs[[k]],
#                             drange=c(d$min, d$max),
#                             grange=c(g$min, g$max),
#                             dab=d$ab,
#                             gab=g$ab)
#           pred = predGPsep(v$GPs[[k]], XX = mvc.data$XT.data$obs$X$trans, lite=T)
#           v$mean = rbind(v$mean, pred$mean)
#           v$var = rbind(v$var, pred$s2)
#           if(!is.null(X.pred) & predict){
#             pred = predGPsep(v$GPs[[k]], XX = X.pred.std$obs$X$trans, lite=T)
#             v.pred$mean = rbind(v.pred$mean, pred$mean)
#             v.pred$var = rbind(v.pred$var, pred$s2)
#           }
#           deleteGPsep(v$GPs[[k]])
#         }
#       }
#       # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
#       v$sample = rmvnorm(1,as.numeric(v$mean),diag(as.numeric(v$var)))
#       dim(v$sample) = dim(v$mean)
#       # update y.resid
#       y.resid = y.resid - w_to_y(v$sample,v$basis$B)
#       if(!is.null(X.pred) & predict){
#         v.pred$sample = rmvnorm(1,as.numeric(v.pred$mean),diag(as.numeric(v.pred$var)))
#         dim(v.pred$sample) = dim(v.pred$mean)
#       }
#     } 
#     ssq.hat = mean(y.resid^2)
#     
#     # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of y.resid
#     ldetS = n.y*log(ssq.hat)
#     Sinv = (1/ssq.hat)*eye(n.y)
#     ll = sum(apply(y.resid,2,function(d) -.5*(ldetS + d%*%Sinv%*%d))) - log(ssq.hat) + sum(dbeta(theta,2,2,log=T))
#     returns = list(ll=ll,ssq=ssq.hat,w.samp=w$sample,v.samp=v$sample)
#     if(!is.null(X.pred) & predict){
#       returns$w.pred = w.pred
#       returns$v.pred = v.pred
#     }
#     return( returns )
#   }
#   
#   n.pc = ncol(mvc.data$obs.basis$B)
#   n = nrow(mvc.data$SC.inputs$X.obs[[1]])
#   n.y = nrow(mvc.data$Y.data$obs$orig)
#   bias = mvc.data$bias
#   p.t = mvc.data$XT.data$p.t
#   get.pred = rep(F,n.samples)
#   
#   if(!is.null(X.pred)){
#     # standardize X.pred
#     X.pred.std = transform_xt(mvc.data$XT.data$sim$X$orig,
#                             mvc.data$XT.data$sim$T$orig,
#                             X.pred)
#     # stretch and compress X.pred
#     X.pred.sc = get_SC_inputs(mvc.data$est.ls,X.pred.std,n.pc)
#     w.pred.arr = array(dim=c(n.pred.samples,ncol(mvc.data$obs.basis$B),nrow(X.pred)))
#     if(bias){v.pred.arr = array(dim=c(n.pred.samples,ncol(mvc.data$obs.basis$D),nrow(X.pred)))}
#     samp.ids = as.integer(seq(n.burn+1,n.samples,length.out=n.pred.samples))
#     ssq.samp.pred = numeric(n.pred.samples)
#     samp.id = 1
#     get.pred[samp.ids] = T
#   }
#   
#   ptm = proc.time()
#   t = t.init
#   llt = ll_mcmc(t,T)
#   ssq = llt$ssq
#   w = llt$w.samp
#   v = llt$v.samp
#   if(!is.null(X.pred)){
#     w.pred = llt$w.pred$sample
#     v.pred = llt$v.pred$sample
#   }
#   ll.prop = list()
#   accept = matrix(0,nrow=1,ncol=p.t)
#   t.store = matrix(nrow=n.samples,ncol=mvc.data$XT.data$p.t);t.store[1,] = t
#   ssq.store = numeric(n.samples);ssq.store[1] = ssq
#   w.store = array(dim=c(n.samples,dim(llt$w.samp)));w.store[1,,] = w
#   v.store = array(dim=c(n.samples,dim(llt$v.samp)))
#   if(bias){v.store[1,,] = v}
#   if(!is.null(X.pred) & get.pred[1]){
#     w.pred.arr[samp.id,,] = w.pred
#     if(bias){v.pred.arr[samp.id,,] = v.pred}
#     ssq.samp.pred[samp.id] = ssq
#     samp.id = samp.id + 1
#   }
# 
#   for(i in 2:n.samples){
#     if(verbose){
#       if(mod(i,n.samples/10)==0){cat(100*(i/n.samples),'% ')}
#     }
#     for(j in 1:p.t){
#       t.prop = t # initialize t.prop to t so that the parameter we are not updating is fixed at its previous value
#       t.prop[j] = t[j] + (prop.step[j] * runif(1,-0.5, 0.5))
#       if(t.prop[j]<0 | t.prop[j]>1){
#         ll.prop$ll = -Inf
#       } else{
#         ll.prop = ll_mcmc(t.prop,get.pred[i])
#       }
#       llt = ll_mcmc(t,get.pred[i])
#       # 4. Accept or reject t
#       if(log(runif(1)) < (ll.prop$ll - llt$ll)){
#         t[j] = t.prop[j]
#         ssq = ll.prop$ssq
#         w = ll.prop$w.samp
#         v = ll.prop$v.samp
#         accept[j] = accept[j] + 1
#         if(!is.null(X.pred) & get.pred[i]){
#           w.pred = ll.prop$w.pred$sample
#           v.pred = ll.prop$v.pred$sample
#         }
#       }
#       t.store[i,j] = t[j]
#     }
#     # store w,v,ssq after both updates have been made
#     ssq.store[i] = ssq
#     w.store[i,,] = w
#     if(bias){v.store[i,,] = v}
#     
#     if(!is.null(X.pred) & get.pred[i]){
#       w.pred.arr[samp.id,,] = w.pred
#       if(bias){v.pred.arr[samp.id,,] = v.pred}
#       ssq.samp.pred[samp.id] = ssq
#       samp.id = samp.id + 1
#     }
#   }
#   mcmc.time = ptm-proc.time()
#   
#   returns = list(t.samp=t.store[(n.burn+1):n.samples,],
#                  ssq.samp=ssq.store[(n.burn+1):n.samples],
#                  w.samp=w.store[(n.burn+1):n.samples,,],
#                  v.samp=NA,
#                  acpt.ratio=accept/n.samples,
#                  time=mcmc.time,
#                  prop.step=prop.step,
#                  n.samples = n.samples,
#                  n.burn = n.burn)
#   if(bias){returns$v.samp=v.store[(n.burn+1):n.samples,,]}
#   if(!is.null(X.pred)){
#     returns$w.pred = w.pred.arr
#     if(bias){returns$v.pred = v.pred.arr}
#     returns$ssq.pred = ssq.samp.pred
#     returns$pred.samp.ids = samp.ids - n.burn
#   }
#   return(returns)
# }

# tune_step_sizes = function(mvc.data,n.burn,n.levels,target.accept.rate=NULL,start=6,end=50,laGP.bias=F,start.bias=6,end.bias=50){
#   p.t = mvc.data$XT.data$p.t
#   
#   # These are the step sizes SEPIA uses, but they seem to result is too small steps for ball drop
#   # step.default = .2
#   # step.sizes = matrix(nrow=n.levels,ncol=p.t)
#   # ex = seq(-(n.levels - 1)/2, (n.levels - 1)/2, length.out=n.levels)
#   # for(i in 1:n.levels){
#   #   for(j in 1:p.t){
#   #     if(ex[i] <= 0){
#   #       base = 2.0
#   #     } else{
#   #       base = 20^(2/(n.levels-1))
#   #     }
#   #     step.sizes[i,j] = step.default * base^ex[i]
#   #   }
#   # }
#   # Instead, I'll start the sequence at .01
#   step.sizes = matrix(logseq(.01,2,n=n.levels),nrow=n.levels,ncol=p.t,byrow = F)
#   
#   acc = matrix(nrow=n.levels,ncol=p.t)
#   for(i in 1:n.levels){
#     acc[i,] = mcmc(mvc.data,t.init = rep(.5,p.t),n.samples = n.burn,prop.step = step.sizes[i,], start=start, end=end,
#                    laGP.bias=laGP.bias,start.bias=start.bias,end.bias=end.bias)$acpt.ratio
#   }
# 
#   # Compute GLM for each parameter
#   step.tuned = numeric(p.t)
#   if(is.null(target.accept.rate)){
#     target.logit = rep(log(1 / (exp(1) - 1)),p.t)
#   } else{
#     target.logit = logit(target.accept.rate)
#   }
#   for(j in 1:p.t){
#       y = cbind(acc[,j]*n.burn,n.burn-(acc[,j]*n.burn))
#       x = cbind(log(step.sizes[,j]))
#       glm.model = glm(y~x,family=binomial(link='logit'))
#       coefs = as.numeric(coef(glm.model))
#       step.tuned[j] = exp((target.logit[j]-coefs[1])/coefs[2])
#   }
#   return(step.tuned)
# }

do_mcmc = function(mvc.data,t.init,n.samples=100,n.burn=0,
                   prop.step=rep(.25,mvc.data$XT.data$p.t),
                   verbose=F,start.eta=6,end.eta=50,
                   laGP.delta=F,start.delta=6,end.delta=50,
                   X.pred=NULL,n.pred.samples=0){
  
  bias = mvc.data$bias
  p.t = mvc.data$XT.data$p.t
  get.pred = rep(F,n.samples)
  X.pred.eta = NULL
  X.pred.delta = NULL
  
  if(!is.null(X.pred)){
    # standardize X.pred
    X.pred.std = transform_xt(mvc.data$XT.data$sim$X$orig,
                              mvc.data$XT.data$sim$T$orig,
                              X.pred)
    # delta needs standardized inputs
    X.pred.delta = X.pred.std$obs$X$trans
    # eta needs stretched and compressed inputs
    X.pred.eta = get_SC_inputs(mvc.data$est.ls,X.pred.std,n.pc)$X.obs
    w.pred.arr = array(dim=c(n.pred.samples,ncol(mvc.data$obs.basis$B),nrow(X.pred)))
    if(bias){
      v.pred.arr = array(dim=c(n.pred.samples,ncol(mvc.data$obs.basis$D),nrow(X.pred)))
      v.GPs = vector(mode='list',length=n.pred.samples)
    }
    samp.ids = as.integer(seq(n.burn+1,n.samples,length.out=n.pred.samples))
    ssq.samp.pred = numeric(n.pred.samples)
    samp.id = 1
    get.pred[samp.ids] = T
  }
  
  ptm = proc.time()
  t = t.init
  
  # some pre-computing to pass around during mcmc
  precomp = list()
  precomp$I.n.y = eye(mvc.data$Y.data$n.y)
  precomp$BD = cbind(mvc.data$obs.basis$B,mvc.data$obs.basis$D)
  precomp$tBD = t(precomp$BD)
  precomp$rankBD = rankMatrix(precomp$BD)
  precomp$BDtBD = precomp$tBD%*%precomp$BD
  if(precomp$rankBD != ncol(precomp$BD)){
    precomp$BDtBDinv = solve(precomp$BDtBD + 1e-8*eye(nrow(precomp$BDtBD)))
  } else{
    precomp$BDtBDinv = solve(precomp$BDtBD) # need to get rid of solves to make code faster
  }
  precomp$ldetBDtBD = determinant(precomp$BDtBD)$modulus
  precomp$Z = lapply(1:mvc.data$sim.basis$n.pc,function(jj) mvc.data$sim.basis$V.t[jj,])
  ####
  
  # remove the things from mvc.data that are expensive to pass around during MCMC and aren't needed
  mvc.data$Y.data$sim = NULL
  mvc.data$Y.data$obs$orig = NULL
  mvc.data$XT.data$sim = NULL; mvc.data$XT.data$obs = NULL
  mvc.data$sim.basis$V.t = NULL
  mvc.data$obs.basis$V.t = NULL
  mvc.data$est.ls$XT = NULL; mvc.data$est.ls$X = NULL
  mvc.data$SC.inputs$X.sim = NULL; mvc.data$SC.inputs$T.sim = NULL; mvc.data$SC.inputs$T.obs = NULL

  llt = fit_model(t,mvc.data,precomp,F,T,start.eta,end.eta,laGP.delta,start.delta,end.delta,X.pred.eta,X.pred.delta)
  
  ssq = llt$ssq.hat
  w = llt$eta$w$sample
  v = llt$delta$v$sample
  if(!is.null(X.pred)){
    w.pred = llt$eta$w.pred$sample
    if(bias){
      v.pred = llt$delta$v.pred$sample
      v.GP = llt$delta$v$GPs
    }
  }
  ll.prop = list()
  accept = matrix(0,nrow=1,ncol=p.t)
  t.store = matrix(nrow=n.samples,ncol=mvc.data$XT.data$p.t);t.store[1,] = t
  ssq.store = numeric(n.samples);ssq.store[1] = ssq
  w.store = array(dim=c(n.samples,dim(w)));w.store[1,,] = w
  v.store = array(dim=c(n.samples,dim(v)))
  if(bias){v.store[1,,] = v}
  if(!is.null(X.pred) & get.pred[1]){
    w.pred.arr[samp.id,,] = w.pred
    if(bias){
      v.pred.arr[samp.id,,] = v.pred
      v.GPs[[samp.id]] = v.GP
    }
    ssq.samp.pred[samp.id] = ssq
    samp.id = samp.id + 1
  }
  
  for(i in 2:n.samples){
    if(verbose){
      if(mod(i,n.samples/10)==0){cat(100*(i/n.samples),'% ')}
    }
    for(j in 1:p.t){
      t.prop = t # initialize t.prop to t so that the parameter we are not updating is fixed at its previous value
      t.prop[j] = t[j] + (prop.step[j] * runif(1,-0.5, 0.5))
      if(t.prop[j]<0 | t.prop[j]>1){
        ll.prop$ll = -Inf
      } else{
        if(get.pred[i]){
          ll.prop = fit_model(t.prop,mvc.data,precomp,F,T,start.eta,end.eta,laGP.delta,start.delta,end.delta,X.pred.eta,X.pred.delta)
        } else{
          ll.prop = fit_model(t.prop,mvc.data,precomp,F,T,start.eta,end.eta,laGP.delta,start.delta,end.delta,NULL,NULL)
        }
      }
      if(get.pred[i]){
        llt = fit_model(t,mvc.data,precomp,F,T,start.eta,end.eta,laGP.delta,start.delta,end.delta,X.pred.eta,X.pred.delta)
      } else{
        llt = fit_model(t,mvc.data,precomp,F,T,start.eta,end.eta,laGP.delta,start.delta,end.delta,NULL,NULL)
      }
      
      # 4. Accept or reject t
      if(log(runif(1)) < (ll.prop$ll - llt$ll)){
        t[j] = t.prop[j]
        ssq = ll.prop$ssq.hat
        w = ll.prop$eta$w$sample
        v = ll.prop$delta$v$sample
        accept[j] = accept[j] + 1
        if(!is.null(X.pred) & get.pred[i]){
          w.pred = ll.prop$eta$w.pred$sample
          if(bias){
            v.pred = ll.prop$delta$v.pred$sample
            v.GP = ll.prop$delta$v$GPs
          }
        }
      }
      t.store[i,j] = t[j]
    }
    # store w,v,ssq after both updates have been made
    ssq.store[i] = ssq
    w.store[i,,] = w
    if(bias){v.store[i,,] = v}
    
    if(!is.null(X.pred) & get.pred[i]){
      w.pred.arr[samp.id,,] = w.pred
      if(bias){
        v.pred.arr[samp.id,,] = v.pred
        v.GPs[[samp.id]] = v.GP
      }
      ssq.samp.pred[samp.id] = ssq
      samp.id = samp.id + 1
    }
  }
  mcmc.time = ptm-proc.time()
  
  returns = list(t.samp=t.store[(n.burn+1):n.samples,],
                 ssq.samp=ssq.store[(n.burn+1):n.samples],
                 w.samp=w.store[(n.burn+1):n.samples,,],
                 v.samp=NULL,
                 v.GPs=NULL,
                 acpt.ratio=accept/n.samples,
                 time=mcmc.time,
                 prop.step=prop.step,
                 n.samples = n.samples,
                 n.burn = n.burn)
  if(bias){returns$v.samp=v.store[(n.burn+1):n.samples,,]}
  if(!is.null(X.pred)){
    returns$w.pred = w.pred.arr
    if(bias){
      returns$v.pred = v.pred.arr
      returns$v.GPs = v.GPs
    }
    returns$ssq.pred = ssq.samp.pred
    returns$pred.samp.ids = samp.ids - n.burn
  }
  return(returns)
}

tune_step_sizes = function(mvc.data,n.samples,n.levels,target.accept.rate=NULL,start.eta=6,end.eta=50,laGP.delta=F,start.delta=6,end.delta=50){
  p.t = mvc.data$XT.data$p.t
  
  # These are the step sizes SEPIA uses, but they seem to result is too small steps for ball drop
  # step.default = .2
  # step.sizes = matrix(nrow=n.levels,ncol=p.t)
  # ex = seq(-(n.levels - 1)/2, (n.levels - 1)/2, length.out=n.levels)
  # for(i in 1:n.levels){
  #   for(j in 1:p.t){
  #     if(ex[i] <= 0){
  #       base = 2.0
  #     } else{
  #       base = 20^(2/(n.levels-1))
  #     }
  #     step.sizes[i,j] = step.default * base^ex[i]
  #   }
  # }
  # Instead, I'll start the sequence at .01
  step.sizes = matrix(logseq(.01,2,n=n.levels),nrow=n.levels,ncol=p.t,byrow = F)
  
  acc = matrix(nrow=n.levels,ncol=p.t)
  for(i in 1:n.levels){
    acc[i,] = do_mcmc(mvc.data,rep(.5,p.t),n.samples,0,step.sizes[i,],F,start.eta,end.eta,laGP.delta,start.delta,end.delta)$acpt.ratio
  }
  
  # Compute GLM for each parameter
  step.tuned = numeric(p.t)
  if(is.null(target.accept.rate)){
    target.logit = rep(log(1 / (exp(1) - 1)),p.t)
  } else{
    target.logit = logit(target.accept.rate)
  }
  for(j in 1:p.t){
    y = cbind(acc[,j]*n.samples,n.samples-(acc[,j]*n.samples))
    x = cbind(log(step.sizes[,j]))
    glm.model = glm(y~x,family=binomial(link='logit'))
    coefs = as.numeric(coef(glm.model))
    step.tuned[j] = exp((target.logit[j]-coefs[1])/coefs[2])
  }
  return(step.tuned)
}

# fit model at t=theta
# X.pred.eta; list of p_eta arrays of shape n.pred x p.x containing sc inputs for prediction
# X.pred.delta; list of p_delta arrays of shape n.pred x p.x containing standardized inputs for prediction
fit_model = function(theta,mvc.data,precomp,
                     lite=T,sample=T,
                     start.eta=6,end.eta=50,
                     laGP.delta=F,start.delta=6,end.delta=50,
                     X.pred.eta=NULL,
                     X.pred.delta=NULL){
  
  # n.pc = mvc.data$sim.basis$n.pc
  # n = mvc.data$Y.data$n
  # n.y = mvc.data$Y.data$n.y
  # bias = mvc.data$bias
  # theta = matrix(theta,nrow=1)
  
  # Predict from emulator at [X.obs,theta]
  eta = fit_eta(theta,mvc.data,precomp,sample,start.eta,end.eta,X.pred.eta)

  if(mvc.data$bias){
    delta = fit_delta(eta$y.resid,mvc.data$obs.basis$D,mvc.data$XT.data,sample=sample,lite=lite,laGP=laGP.delta,start=start.delta,end=end.delta,X.pred.delta)
  } else{
    delta = NULL
  }
  
  ll = compute_ll(theta,eta,delta,sample,precomp)
  
  if(lite){
    return(ll)
  } else{
    return(list(ll = ll,
                theta = theta,
                ssq.hat = ifelse(is.null(delta),eta$ssq.hat,delta$ssq.hat),
                eta = eta,
                delta = delta))
  }
}

# fit modular emulator with t=theta and optionally return predictions at centered and scaled inputs X.pred.sc
fit_eta = function(theta,mvc.data,precomp,
                   sample=T,start=6,end=50,X.pred.sc=NULL){
  
  n = mvc.data$Y.data$n
  n.pc = mvc.data$sim.basis$n.pc
  p.t = mvc.data$XT.data$p.t
  p.x = mvc.data$XT.data$p.x
  
  # transform theta by estimated length-scales
  # theta.sc = lapply(1:n.pc, function(jj) sc_inputs(theta,mvc.data$est.ls$T[[jj]]))
  theta.sc = lapply(1:n.pc, function(jj) matrix(theta/mvc.data$est.ls$T[[jj]],nrow=n,ncol=p.t,byrow = T))
  # make prediction matrix from X.obs and theta
  if(p.t>1){
    # XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(replicate(n,theta.sc[[jj]],simplify='matrix'))))
    if(p.x>0){
      XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],theta.sc[[jj]]))
    } else{
      XT.sc = theta.sc
    }
  } else{
    XT.sc = lapply(1:n.pc, function(jj) cbind(mvc.data$SC.inputs$X.obs[[jj]],t(t(replicate(n,theta.sc[[jj]],simplify='matrix')))))
  }

  w = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
                   Z=precomp$Z,
                   XX=XT.sc,
                   start=start,
                   end=end,
                   sample=sample)

  if(!is.null(X.pred.sc)){
    # need to fix this
    if(p.t>1){
      XT.pred.sc = lapply(1:n.pc, function(jj) cbind(X.pred.sc[[jj]],t(replicate(nrow(X.pred.sc[[jj]]),theta.sc[[jj]],simplify='matrix'))))
    } else{
      XT.pred.sc = lapply(1:n.pc, function(jj) cbind(X.pred.sc[[jj]],t(t(replicate(nrow(X.pred.sc[[jj]]),theta.sc[[jj]],simplify='matrix')))))
    }
    w.pred = aGPsep_SC_mv(X=mvc.data$SC.inputs$XT.sim,
                          Z=precomp$Z,
                          XX=XT.pred.sc,
                          start=start,
                          end=end)
  } else{
    w.pred = NULL
  }
  
  # residuals
  y.resid = mvc.data$Y.data$obs$trans - w_to_y(w$sample,mvc.data$obs.basis$B)
  ssq.hat = mean(y.resid^2)
  
  return(list(w=w,B=mvc.data$obs.basis$B,
              y.resid=y.resid,ssq.hat=ssq.hat,
              w.pred = w.pred))
}

# fit delta model to residuals with basis D
fit_delta = function(y.resid,D,XT.data,sample=T,lite=T,laGP=F,start=6,end=50,X.pred.std=NULL){
  v = NULL; v.pred = NULL
  v$basis = get_basis(y.resid,pct.var=.95,B=D)
  if(laGP){ # use laGP for a faster bias model
    v = append(v,aGPsep_SC_mv(X=rep(list(XT.data$obs$X$trans),v$basis$n.pc),
                              Z=lapply(1:v$basis$n.pc,function(jj) v$basis$V.t[jj,]),
                              XX=rep(list(XT.data$obs$X$trans),v$basis$n.pc),
                              start=start,
                              end=end,
                              bias=T))
  } else{       # use full GP for the bias model
    v$GPs = vector(mode='list',length=v$basis$n.pc)
    
    for(k in 1:v$basis$n.pc){
      d = darg(d=list(mle=T),X=XT.data$obs$X$trans)
      g = garg(g=list(mle=T),y=v$basis$V.t[k,])
      # Our responses have changed wrt to inputs, so we should not use d=1 an.ymore
      v$GPs[[k]] <- newGPsep(X=XT.data$obs$X$trans,
                             Z=v$basis$V.t[k,],
                             d=d$start, g=g$start, dK=TRUE)
      cmle <- jmleGPsep(v$GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      pred = predGPsep(v$GPs[[k]], XX = XT.data$obs$X$trans, lite=T)
      v$mean = rbind(v$mean, pred$mean)
      v$var = rbind(v$var, pred$s2)
      if(!is.null(X.pred.std)){
        pred = predGPsep(v$GPs[[k]], XX = X.pred.std, lite=T)
        v.pred$mean = rbind(v.pred$mean, pred$mean)
        v.pred$var = rbind(v.pred$var, pred$s2)
        if(sample){
          v.pred$sample = rmvnorm(1,as.numeric(v.pred$mean),diag(as.numeric(v.pred$var)))
          dim(v.pred$sample) = dim(v.pred$mean)
        } else{
          v.pred$sample = v.pred$mean
        }
      }
      if(lite){
        deleteGPsep(v$GPs[[k]])
      }
    }
  }
  # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
  if(sample){
    v$sample = rmvnorm(1,as.numeric(v$mean),diag(as.numeric(v$var)))
    dim(v$sample) = dim(v$mean)
  } else{
    v$sample = v$mean
  }
  y.resid = y.resid - w_to_y(v$sample,v$basis$B)
  ssq.hat = mean(y.resid^2)
  return(list(v=v,D=D,y.resid=y.resid,ssq.hat=ssq.hat,
              v.pred=v.pred))
}

# w,v are needed if sample=F
compute_ll = function(theta,eta,delta,sample,precomp){
  #BD = cbind(eta$B,delta$D)
  if(!is.null(delta)){
    ssq.hat = delta$ssq.hat
    y.resid = delta$y.resid
  } else{
    ssq.hat = eta$ssq.hat
    y.resid = eta$y.resid
  }
  if(sample){
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of y.resid
    ldetS = precomp$n.y*log(ssq.hat)
    Sinv = precomp$I.n.y/ssq.hat
    ll = apply(y.resid,2,function(d) -.5*(ldetS + d%*%Sinv%*%d))
  } else{
    # r = 1e-8
    lambda = 1/ssq.hat
    # tBD = t(BD)
    # rankBD = rankMatrix(BD)
    # BDtBD = tBD%*%BD
    # if(rankBD != ncol(BD)){
    #   BDtBDinv = solve(BDtBD + r*eye(nrow(BDtBD)))
    # } else{
    #   BDtBDinv = solve(BDtBD) # need to get rid of solves to make code faster
    # }
    # ldetBDtBD = determinant(BDtBD)$modulus
    # I.n.y = eye(n.y)
    B.hat = precomp$BDtBDinv%*%precomp$tBD%*%y.resid
    
    SigB.hat = lapply(1:n, function(i) as.matrix(diag(c(w$var[,i],v$var[,i]),nrow=length(c(w$var[,i],v$var[,i]))) + ssq.hat*precomp$BDtBDinv))
    l_beta_hat = sapply(1:n, function(i) dmvnorm(x=as.numeric(B.hat[,i]),sigma = SigB.hat[[i]], log = T))
    ll = sapply(1:n, function(i) -.5*((precomp$n.y-precomp$rankBD)*log(2*pi*ssq.hat) + precomp$ldetBDtBD +
                                        lambda*t(y.resid[,i])%*%(precomp$I.n.y-precomp$BD%*%precomp$BDtBDinv%*%precomp$tBD)%*%y.resid[,i]) + 
                  l_beta_hat[i])
  }
  return(sum(ll) + log(ssq.hat) - sum(dbeta(theta,2,2,log=T)))
}
