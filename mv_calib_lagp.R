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

plot_wz_pairs = function(w,z,legend=T){
  nPC = nrow(w)
  m = ncol(w)
  n = ncol(z)
  if(nPC>1){
    pairs(rbind(t(w),t(z)),col=c(rep('blue',m),rep('orange',n)),pch=1,labels=paste0('PC',1:nPC),
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
mv_eta_predict = function(theta,XpredOrig,mvcData){
  
  nPC = nrow(mvcData$simBasis$Vt)
  # transform and scale Xpred
  Xpred = transform_xt(mvcData$XTdata$sim$X$orig,mvcData$XTdata$sim$T$orig,XpredOrig,Tobs=matrix(rep(theta,nrow(XpredOrig)),nrow=nrow(XpredOrig),byrow=T))
  XpredSC = get_SC_inputs(mvcData$estLS,Xpred,nPC)
  eta = aGPsep_SC_mv(X=XpredSC$XTsim,
                       Z=lapply(1:nPC,function(i) mvcData$simBasis$Vt[i,]),
                       XX=lapply(1:nPC,function(i) cbind(XpredSC$Xobs[[i]],XpredSC$Tobs[[i]])),
                       returnYmean = F)
  return(eta)
}

mv_delta_predict = function(XpredOrig,calib,mvcData){
  
  Xpred = unit_xform(XpredOrig,Xmin=mvcData$XTdata$sim$X$min,Xrange=mvcData$XTdata$sim$X$range)
  
  delta = lapply(1:calib$delta$basis$nPC, function(i) predGPsep(calib$delta$GPs[[i]],Xpred$trans,lite=T))
  vMean = NULL; vVar = NULL
  for(i in 1:calib$delta$basis$nPC){
    vMean = rbind(vMean,delta[[i]]$mean)
    vVar = rbind(vVar,delta[[i]]$s2)
  }
  return(list(vMean=vMean,vVar=vVar))
}
ypred_mle = function(XpredOrig,mvcData,calib,nsamples=1,support='obs',bias=F,emOnly=F){
  
  # emulator predictions
  eta = mv_eta_predict(calib$theta.hat.orig,XpredOrig,mvcData)
  
  if(support=='obs'){
    B = mvcData$obsBasis$B
    D = mvcData$obsBasis$D
    ym = mvcData$Ydata$obs$mean
    ysd = mvcData$Ydata$obs$sd
    nY = nrow(mvcData$Ydata$obs$orig)
  } else{
    B = mvcData$simBasis$B
    D = mvcData$simBasis$D
    ym = mvcData$Ydata$sim$mean
    ysd = mvcData$Ydata$sim$sd
    nY = nrow(mvcData$Ydata$sim$orig)
  }
  n = nrow(XpredOrig)

  if(!bias){
    # unbiased prediction
    ySamp = array(dim=c(nY,nsamples,n))
    for(i in 1:n){
      ySamp[,,i] = t(mvtnorm::rmvnorm(nsamples,mean=B%*%eta$wMean[,i,drop=F],
                                      sigma = calib$ssq.hat*eye(nY) + B%*%diag(eta$wVar[,i])%*%t(B))) * 
        ysd + ym
    }
    yMean = apply(ySamp,c(1,3),mean)
    yVar = apply(ySamp,c(1,3),var)
    yInt = lapply(1:n,function(i) t(apply(ySamp[,,i],1,quantile,c(.025,.975))))
    return(list(ySamp=ySamp,yMean=yMean,yVar=yVar,yInt=yInt))
  } else{
    # biased prediction
    etaSamp = array(dim=c(nY,nsamples,n))
    deltaSamp = array(dim=c(nY,nsamples,n))
    delta = mv_delta_predict(XpredOrig,calib,mvcData)
    for(j in 1:n){
      etaSamp[,,j] = t(mvtnorm::rmvnorm(nsamples,mean=B%*%eta$wMean[,j,drop=F],
                                       sigma = B%*%diag(eta$wVar[,j])%*%t(B))) *
        ysd + ym
      deltaSamp[,,j] = t(mvtnorm::rmvnorm(nsamples,mean=D%*%delta$vMean[,j,drop=F],
                                         sigma = D%*%diag(delta$vVar[,j],nrow=length(delta$vVar[,j]))%*%t(D) *
                                           calib$ssq.hat*eye(nY) )) * ysd
    }

    ySamp = etaSamp + deltaSamp
    yMean = apply(ySamp,c(1,3),mean)
    yVar = apply(ySamp,c(1,3),var)
    etaMean = apply(etaSamp,c(1,3),mean)
    deltaMean = apply(deltaSamp,c(1,3),mean)
    yInt = lapply(1:n,function(i) t(apply(ySamp[,,i],1,quantile,c(.025,.975))))
    return(list(ySamp=ySamp,yMean=yMean,yVar=yVar,yInt=yInt,
                etaMean=etaMean,deltaMean=deltaMean))
  }
}
ypred_mcmc = function(XpredOrig,mvcData,tsamp,support='obs',bias=F){
  if(support=='obs'){
    B = mvcData$obsBasis$B
    D = mvcData$obsBasis$D
    ym = mvcData$Ydata$obs$mean
    ysd = mvcData$Ydata$obs$sd
    nY = nrow(mvcData$Ydata$obs$orig)
  } else{
    B = mvcData$simBasis$B
    D = mvcData$simBasis$D
    ym = mvcData$Ydata$sim$mean
    ysd = mvcData$Ydata$sim$sd
    nY = nrow(mvcData$Ydata$sim$orig)
  }
 
  ysamp = array(0,dim=c(nrow(tsamp),nY,nrow(XpredOrig)))
  for(i in 1:nrow(tsamp)){
    tmp = mv.calib(mvcData,tsamp[i,],bias=bias,sample=T,optimize=F)
    if(bias){
      for(j in 1:length(tmp$delta$GPs)){
        deleteGPsep(tmp$delta$GPs[[j]])
      }
      mu = B%*%tmp$eta$wSamp + D%*%tmp$delta$wSamp
    } else{
      mu = B%*%tmp$eta$wSamp
    }
    sigma = tmp$ssq.hat*eye(nY)
    for(j in 1:nrow(XpredOrig)){
        ysamp[i,,j] = rmvnorm(1,mean=mu[,j],sigma=sigma) * ysd + ym
    }
  }
  mean = apply(ysamp,2:3,mean)
  var = apply(ysamp,2:3,var)
  conf.int = apply(ysamp,2:3,quantile,c(.025,.975))
  return(list(mean=mean,var=var,conf.int=conf.int))
}
# FUNCTION: multivariate aGPsep for use with stretched and compressed inputs only
{## Parameters:
# required
# X : list of length nPC which contains stretched and compressed training input matrices
# Z : list of length nPC which contains training response vectors
# XX: list of length nPC with contains stretched and compressed prediction input matrices
# optional
# start: integer, initial size of laPG neighborhood
# end: integer, final size of laGP neighborhood
# g: double, laGP nugget
# returnYmean: bool, return GP mean on Y scale
# yMean: vector, mean of original y's, if nothing is passed yMean is returned unshifted
# ySD: vector, sd of original y's, if nothing is passed yMean is returned unscaled
# Bobs: matrix of basis vectors for transforming predicted basis weights to y scale
## Returns:
# lagp_fit: list of nPC outputs from aGPsep
# wMean: matrix (nPC x n) of prediction means
# wVar: matrix (nPC x n) of prediction variances
# yMean: matrix (nY x n) of prediction means on y scale (Bobs %*% wMean)
## Calls:
# aGPsep from laGP package
# w_to_y to transform w -> y
## Called by:
# mv.calib.bias: biased calibration
# mv.calib.nobias: unbiased calibration
# mv.em.pred: laGP emulator prediction function at new inputs
}
aGPsep_SC_mv = function(X, Z, XX,
                        start=6, end=50, 
                        g=1/10000,
                        returnYmean=F, yMean = NULL, ySD = NULL, Bobs=NULL){
  nPC = length(Z)
  nY = nrow(XX[[1]])

  lagp_fit = lapply(1:nPC,function(i) aGPsep(X = X[[i]],
                                             Z = Z[[i]],
                                             XX = XX[[i]],
                                             d = list(mle = FALSE, start = 1),
                                             g = g,
                                             start = start,
                                             end = min(end,length(Z[[i]])-1),
                                             method = 'nn',
                                             verb=0))
  wMean = array(dim=c(nPC,nY))
  wVar = array(dim=c(nPC,nY))
  for (i in 1:nPC) {
    wMean[i,] = lagp_fit[[i]]$mean
    wVar[i,] = lagp_fit[[i]]$var
  }
  return = list(lagp_fit=lagp_fit,wMean=wMean,wVar=wVar)

  # return results on Y scale?
  if(returnYmean){
    return$yMean = w_to_y(wMean,Bobs,yMean,ySD)
  }
  
  return(return)
}

# FUNCTION: get_basis
{### Accepts a matrix of response values Y, and returns the basis decomposition from SVD
## Parameters:
# Y (numeric)       : matrix of response values to decompose into basis representation. Y should be
#                     formatted such that each row represents the multivariate response for a single experiment
# nPC (int)        : number of principal components desired from basis decomposition
# pctVar (numeric) : desired percentage of variability explained by PC's. Overrides nPC if specified
## Returns: List containings the following
# B (numeric)       : matrix of nPC basis vectors for Ytrans
# Vt (numeric)      : matrix of nPC response vectors for Ytrans
# nPC (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating sim basis
# mv.calib.bias for generating discrepancy basis
}
get_basis = function(Y, nPC = 1, pctVar = NULL, fullBasis=F, B=NULL){

  retList = list()
  
  # Compute SVD
  nExp = ncol(Y)
  if(is.null(B)){
    svdY = svd(Y)
    retList$B = svdY$u %*% diag(svdY$d) / sqrt(nExp)
    retList$Vt = t(svdY$v) * sqrt(nExp)
  
    # Return only the required number of basis vectors
    percentVar = zapsmall((svdY$d) / sum(svdY$d))
    if(!is.null(pctVar)){
      # percentage of variance explained specified
      nPC = which.max(cumsum(percentVar)>=pctVar)
    }
    
    retList$nPC = nPC
    if(fullBasis){
      retList$pctVar = percentVar
      retList$B = as.matrix(retList$B)
      retList$Vt = as.matrix(retList$Vt)
    } else{
      retList$pctVar = percentVar[1:retList$nPC]
      retList$B = as.matrix(retList$B[,1:retList$nPC,drop=F])
      retList$Vt = as.matrix(retList$Vt[1:retList$nPC,,drop=F])
    }
  } else{
    # manual basis matrix was passed
    retList$B = B
    retList$nPC = ncol(B)
    retList$Vt = solve(t(B)%*%B)%*%t(B)%*%Y
  }
  return(retList)
}

# FUNCTION: get_obs_basis
{### Generated obs basis vectors by interpolating sim basis vectors
## Parameters:
# simBasis       : output from get_basis
# Yobs        : matrix (nY x n) field data
# YindSim : matrix of field data indices of observation
## Returns: List containings the following
# B (numeric)       : matrix of nPC basis vectors
# Vt (numeric)      : matrix of nPC response vectors
# nPC (integer)    : number of PC's used for the basis decomposition
## Calls:
## Called by:
# user when generating obs basis
}
get_obs_basis = function(simBasis,Yobs,YindSim,YindObs,sigY=NULL){
  obsBasis = list()
  nPC = simBasis$nPC
  obsBasis$nPC = nPC
  if(!isTRUE(all.equal(YindObs,YindSim))){
    B = matrix(nrow=nrow(YindObs),ncol=obsBasis$nPC)
    if(ncol(YindSim)==1){
      # 1d interpolation
      for(j in 1:nPC){
        B[,j] = interp1(as.numeric(YindSim),as.numeric(simBasis$B[,j]),as.numeric(YindObs))
      }
    } else if(ncol(YindSim)==2){
      # 2d interpolation
      y = unique(YindSim[,1]); x = unique(YindSim[,2])
      yp = YindObs[,1]; xp = YindObs[,2]
      for(i in 1:nPC){
        B[,i] = interp2(x,y,matrix(simBasis$B[,i],length(y),length(x),byrow = F),xp,yp)
      }
    }
    obsBasis$B = B
  } else{
    cat('setting Bobs==Bsim')
    obsBasis$B = simBasis$B
  }
  # compute Vt
  # Bridge = 1e-6 * diag(rep(1,nPC)) # in case Bprod is ill-conditioned
  # if(!is.null(sigY)){
  #   lamY = (1/sigY)*eye(dim(obsBasis$B)[1])
  # } else{
  #   lamY = eye(dim(obsBasis$B)[1])
  # }
  # Bprod = t(obsBasis$B)%*%lamY%*%obsBasis$B
  # obsBasis$Vt = solve(Bprod + Bridge) %*% t(obsBasis$B)%*%lamY%*%Yobs
  obsBasis$Vt = solve(t(obsBasis$B)%*%obsBasis$B)%*%t(obsBasis$B)%*%Yobs
  return(obsBasis)
}

# FUNCTION: mv_lengthscales
# TO DO: this function needs to be updated for only using a subset of the data
### Computes estimated length scale parameters for use in stretching and compressing the inputs space
# Parameters:
# X (numeric) : untransformed matrix of inputs
# Vt (numeric) : matrix 
# ls_prior_ab (numeric) : 2 vector of parameters for gamma prior on length-scales
mv_lengthscales = function(XTdata,Vt,g=1e-7,ls_prior_ab=NULL){

  XT = cbind(XTdata$sim$X$trans,XTdata$sim$T$trans)
  pX = XTdata$pX
  pT = XTdata$pT
  nPC = nrow(Vt)
  
  dConfig = darg(list(mle = TRUE, max = 100), XT)
  if(is.null(ls_prior_ab)){
    ls_prior_ab = dConfig$ab
  }
  GPlist = vector(mode='list',length=nPC)
  for(i in 1:nPC){
    GPlist[[i]] = newGPsep(X = XT,
                           Z = Vt[i, ],
                           d = rep(dConfig$start, ncol(XT)),
                           g = g,
                           dK = TRUE)
  }
  
  estLenscales = lapply(1:nPC, function(i) mleGPsep(GPlist[[i]],param='d',
                                                        tmin = dConfig$min, 
                                                        tmax = dConfig$max,
                                                        ab=ls_prior_ab,
                                                        maxit=200)$d)
  for(i in nPC){
    deleteGPsep(GPlist[[i]])
  }
  
  return = list()
  return$XT = estLenscales
  return$X = lapply(1:nPC, function(i) estLenscales[[i]][1:pX])
  return$T = lapply(1:nPC, function(i) estLenscales[[i]][(pX+1):(pX+pT)])
  return(return)
}

# FUNCTION: unit_xform
### Transform columns to lie in unit interval
# Parameters
# X (numeric) : matrix to be transformed
# Returns
# Xtrans (numeric) : transformed matrix
unit_xform = function(X,Xmin=NULL,Xrange=NULL){
  X = as.matrix(X)
  if(is.null(Xmin)){
    Xmin = apply(X,2,min)
  }
  if(is.null(Xrange)){
    Xmax = apply(X,2,max)
    Xrange = Xmax - Xmin
  }
  if(!isTRUE(all.equal(Xrange,rep(0,length(Xrange))))){
    Xtrans = t( (t(X) - Xmin) / Xrange )
  } else{
    # just return X if dummyx
    Xtrans = X
  }
  return(list(orig=X,
              trans=Xtrans,
              min=Xmin,
              range=Xrange))
}

transform_y = function(Ysim,YindSim,Yobs=NULL,YindObs=NULL,center=T,scale=T, scaletype='rowwise'){
  simList = list()
  obsList = list()
  dim = ncol(YindSim)
  
  # scale Ysim
  simList$orig = Ysim
  if(center){
    simList$mean = rowMeans(Ysim)
  } else{
    simList$mean = rep(0,nrow(Ysim))
  }
  if(scale){
    if(scaletype=='rowwise'){
      simList$sd = apply(Ysim,1,sd)
      # dont scale places with var 0
      simList$sd[simList$sd==0]=1
    } else if(scaletype=='scalar') {
      sd = sd(Ysim)
      simList$sd = rep(sd,nrow(Ysim))
    }
  } else{
    simList$sd = rep(1,nrow(Ysim))
  }
  simList$trans = t(scale(t(Ysim),simList$mean,simList$sd))
  
  # scale Yobs
  if(!is.null(Yobs)){
    obsList$orig = Yobs
    if(!isTRUE(all.equal(YindObs,YindSim))){
      # interpolation is needed
      if(dim==1){
        # 1d interpolation of sim mean onto yobs support
        obsList$mean = interp1(as.numeric(YindSim),simList$mean,as.numeric(YindObs))
        obsList$sd = interp1(as.numeric(YindSim),simList$sd,as.numeric(YindObs))
        obsList$trans = t(scale(t(Yobs),obsList$mean,obsList$sd))
      } else if(dim==2){
        # 2d interpolation of sim mean onto yobs support
        x = unique(YindSim[,1]); y = unique(YindSim[,2])
        xp = YindObs[,1]; yp = YindObs[,2]
        obsList$mean = interp2(x,y,matrix(simList$mean,length(y),length(x),byrow = T),xp,yp)
        if(scale){
          if(scaletype=='rowwise'){
            obsList$sd = interp2(x,y,matrix(simList$sd,length(y),length(x),byrow = T),xp,yp)
          } else if(scaletype=='scalar'){
            obsList$sd = rep(sd,nrow(Yobs))
          }
        } else{
          obsList$sd = rep(1,nrow(Yobs))
        }
        obsList$trans = t(scale(t(Yobs),obsList$mean,obsList$sd))
      }
    } else{
      obsList$trans = t(scale(t(Yobs),simList$mean,simList$sd)) 
      obsList$mean = simList$mean
      obsList$sd = simList$sd
    }
  }
  return(list(sim=simList,obs=obsList))
}

# FUNCTION: transform_xt
### Transforms inputs to lie on the unit hypercube
# Tobs are optional known calibration parameters for the observed data
# transforming these is useful for checking predictions at obs locations
transform_xt = function(Xsim,Tsim=NULL,Xobs=NULL,Tobs=NULL){
  pX = ncol(Xsim)
  if(!is.null(Tsim)){
    pT = ncol(Tsim)
  } else{
    pT = 0
  }
  simList = list()
  simList$X = unit_xform(Xsim)
  simList$T = NULL
  if(!is.null(Tsim)){
    simList$T = unit_xform(Tsim)
  }
  obsList = list()
  obsList$X = NULL
  obsList$T = NULL
  if(!is.null(Xobs)){
    obsList$X = unit_xform(Xobs,Xmin=simList$X$min,Xrange=simList$X$range)
  }
  if(!is.null(Tobs)){
    obsList$T = unit_xform(Tobs,Xmin=simList$T$min,Xrange=simList$T$range)
  }
  return(list(sim=simList,obs=obsList,pX=pX,pT=pT))
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

get_SC_inputs = function(estLS,XTdata,nPC){
  return = list()
  return$Xsim = lapply(1:nPC, function(i) sc_inputs(XTdata$sim$X$trans,estLS$X[[i]]))
  return$Tsim =  lapply(1:nPC, function(i) sc_inputs(XTdata$sim$T$trans,estLS$T[[i]]))
  return$XTsim = lapply(1:nPC, function(i) cbind(return$X[[i]],return$T[[i]]))
  
  if(!is.null(XTdata$obs$X)){
    return$Xobs = lapply(1:nPC, function(i) sc_inputs(XTdata$obs$X$trans,estLS$X[[i]]))
  }
  if(!is.null(XTdata$obs$T)){
    return$Tobs = lapply(1:nPC, function(i) sc_inputs(XTdata$obs$T$trans,estLS$T[[i]]))
  }
  return$pX = XTdata$pX
  return$pT = XTdata$pT
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

negll = function(theta,mvcData,lite=T,bias=T,sample=T){
  nPC = nrow(mvcData$obsBasis$Vt)
  n = nrow(mvcData$SCinputs$Xobs[[1]])
  nY = nrow(mvcData$Ydata$obs$orig)
  theta = matrix(theta,nrow=1)
  
  eta = NULL; delta = NULL;
  ll = numeric(n)
  
  # transform theta by estimated length-scales
  thetaSC = lapply(1:nPC, function(j) sc_inputs(theta,mvcData$estLS$T[[j]]))
  # make prediction matrix from Xobs and theta
  if(mvcData$SCinputs$pT>1){
    Xpred = lapply(1:nPC, function(j) cbind(mvcData$SCinputs$Xobs[[j]],t(replicate(n,thetaSC[[j]],simplify='matrix'))))
  } else{
    Xpred = lapply(1:nPC, function(j) cbind(mvcData$SCinputs$Xobs[[j]],t(t(replicate(n,thetaSC[[j]],simplify='matrix')))))
  }
  
  # Predict from emulator at [Xobs,theta]
  eta = aGPsep_SC_mv(X=mvcData$SCinputs$XTsim,
                     Z=lapply(1:nPC,function(j) mvcData$simBasis$Vt[j,]),
                     XX=Xpred,
                     Bobs = mvcData$obsBasis$B)
  
  if(sample){
    eta$wSamp = rmvnorm(1,as.numeric(eta$wMean),diag(as.numeric(eta$wVar)))
    dim(eta$wSamp) = dim(eta$wMean)
  } else{
    eta$wSamp = eta$wMean
  }

  # residuals
  yDisc = mvcData$Ydata$obs$trans - w_to_y(eta$wSamp,mvcData$obsBasis$B)
  ssq.hat = mean(yDisc^2)
  if(!lite){
    yDiscOrig = mvcData$Ydata$obs$orig - w_to_y(eta$wSamp,mvcData$obsBasis$B,mvcData$Ydata$obs$mean,mvcData$Ydata$obs$sd)
    ssq.hat.orig = mean(yDiscOrig^2)
  }
  if(bias){
    delta$basis = get_basis(yDisc,pctVar=.95,B=mvcData$obsBasis$D)
    delta$GPs = vector(mode='list',length=delta$basis$nPC)
    
    for(k in 1:delta$basis$nPC){
      d = darg(d=list(mle=T),X=mvcData$XTdata$obs$X$trans)
      g = garg(g=list(mle=T),y=delta$basis$Vt[k,])
      # Our responses have changed wrt to inputs, so we should not use d=1 anymore
      delta$GPs[[k]] <- newGPsep(X=mvcData$XTdata$obs$X$trans,
                                 Z=delta$basis$Vt[k,],
                                 d=d$start, g=g$start, dK=TRUE)
      cmle <- jmleGPsep(delta$GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      pred = predGPsep(delta$GPs[[k]], XX = mvcData$XTdata$obs$X$trans, lite=T)
      delta$wMean = rbind(delta$wMean, pred$mean)
      delta$wVar = rbind(delta$wVar, pred$s2)
      if(lite){
        deleteGPsep(delta$GPs[[k]])
      }
    }
    # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
    if(sample){
      delta$wSamp = rmvnorm(1,as.numeric(delta$wMean),diag(as.numeric(delta$wVar)))
      dim(delta$wSamp) = dim(delta$wMean)
    } else{
      delta$wSamp = delta$wMean
    }
    yDisc = yDisc - w_to_y(delta$wSamp,delta$basis$B)
    ssq.hat = mean(yDisc^2)
    if(!lite){
      newDiscOrig = yDiscOrig - w_to_y(delta$wSamp,delta$basis$B,sd=mvcData$Ydata$obs$sd)
      ssq.hat.orig = mean(newDiscOrig^2)
    }
    nY = nrow(yDisc)
    BD = cbind(mvcData$obsBasis$B,delta$basis$B)
  } else{
    BD = mvcData$obsBasis$B
  }
  
  if(sample){
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of yDisc
    ldetS = nY*log(ssq.hat)
    Sinv = (1/ssq.hat)*eye(nY)
    ll = apply(yDisc,2,function(d) -.5*(ldetS + d%*%Sinv%*%d))
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
    eyenY = eye(nY)
    Bhat = BDtBDinv%*%tBD%*%yDisc
    
    SigBhat = lapply(1:n, function(i) as.matrix(diag(c(eta$wVar[,i],delta$wVar[,i])) + ssq.hat*BDtBDinv))
    l_beta_hat = sapply(1:n, function(i) dmvnorm(x=as.numeric(Bhat[,i]),sigma = SigBhat[[i]], log = T))
    ll = sapply(1:n, function(i) -.5*((nY-rankBD)*log(2*pi*ssq.hat) + ldetBDtBD +
                                                  lambda*t(yDisc[,i])%*%(eyenY-BD%*%BDtBDinv%*%tBD)%*%yDisc[,i]) + 
                                           l_beta_hat[i])
  }
  if(lite){
    return(-sum(ll))
  } else{
    return = list()
    return$ssq.hat = ssq.hat
    return$ssq.hat.orig = ssq.hat.orig
    return$eta = eta
    return$delta = delta
    return(return)
  }
}



calib.init = function(tInit,mvcData,bias,sample=F){
  pT = mvcData$XTdata$pT
  # tInitOrig = t(t(tInit) * mvcData$XTdata$sim$T$range + mvcData$XTdata$sim$T$min)

  ll = sapply(1:nrow(tInit), function(j) -negll(tInit[j,,drop=F],mvcData,lite=T,bias=bias,sample=sample))
  
  ll_df = as.data.frame(matrix(cbind(ll,as.matrix(tInit)),nrow=length(ll),ncol=pT+1))
  names = c("ll")
  names = c(names,paste0("theta",1:pT))
  colnames(ll_df) = names
  ll_df = ll_df[order(ll_df$ll,decreasing=T),]
  
  return(ll_df)
}

mv.calib = function(mvcData,init,nrestarts=1,bias=F,sample=T,optimize=T,
                    opts=list("MAX_BB_EVAL" = 1000, "INITIAL_MESH_SIZE" = .1,
                              "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0)){

  if(optimize){
    # nrestarts cannot be greater than our set of initial values
    if(nrestarts>nrow(init)){
      cat('nrestarts>nrow(init), setting nrestarts=nrow(init)')
      nrestarts=nrow(init)
    }
    out = NULL
    for(i in 1:nrestarts){
      outi = snomadr(negll, n = mvcData$SCinputs$pT, bbin = rep(0, mvcData$SCinputs$pT), bbout = 0, 
                     x0 = as.numeric(init[i,2:ncol(init)]),
                     lb = rep(.01, mvcData$SCinputs$pT), 
                     ub = rep(.99, mvcData$SCinputs$pT),
                     opts = opts,
                     print.output=F,
                     mvcData = mvcData,
                     lite=T,
                     bias=bias,
                     sample=sample)
      if(is.null(out) || outi$objective < out$objective) out <- outi
    }
    return = negll(matrix(out$solution,nrow=1),mvcData,lite=F,bias=bias,sample=sample)
    return$snomadr = out
    return$theta.hat = out$solution
    return$theta.hat.orig = out$solution * mvcData$XTdata$sim$T$range + mvcData$XTdata$sim$T$min
  } else{
    return = negll(init,mvcData,lite=F,bias=bias,sample=sample)
    return$theta.hat = init
    return$theta.hat.orig = init * mvcData$XTdata$sim$T$range + mvcData$XTdata$sim$T$min
  }
  return(return)
}

mcmc = function(mvcData,tInit,bias=F,nsamples=100,nburn=10,prop.step=rep(.025,mvcData$XTdata$pT),verbose=F){
  ll_mcmc = function(theta){
    nPC = nrow(mvcData$obsBasis$Vt)
    n = nrow(mvcData$SCinputs$Xobs[[1]])
    nY = nrow(mvcData$Ydata$obs$orig)
    theta = matrix(theta,nrow=1)
    
    eta = NULL; delta = NULL;
    ll = numeric(n)
    
    # transform theta by estimated length-scales
    thetaSC = lapply(1:nPC, function(j) sc_inputs(theta,mvcData$estLS$T[[j]]))
    # make prediction matrix from Xobs and theta
    if(mvcData$SCinputs$pT>1){
      Xpred = lapply(1:nPC, function(j) cbind(mvcData$SCinputs$Xobs[[j]],t(replicate(n,thetaSC[[j]],simplify='matrix'))))
    } else{
      Xpred = lapply(1:nPC, function(j) cbind(mvcData$SCinputs$Xobs[[j]],t(t(replicate(n,thetaSC[[j]],simplify='matrix')))))
    }
    
    # Predict from emulator at [Xobs,theta]
    eta = aGPsep_SC_mv(X=mvcData$SCinputs$XTsim,
                       Z=lapply(1:nPC,function(j) mvcData$simBasis$Vt[j,]),
                       XX=Xpred,
                       Bobs = mvcData$obsBasis$B)
    eta$wSamp = rmvnorm(1,as.numeric(eta$wMean),diag(as.numeric(eta$wVar)))
    dim(eta$wSamp) = dim(eta$wMean)
    
    # residuals
    yDisc = mvcData$Ydata$obs$trans - w_to_y(eta$wSamp,mvcData$obsBasis$B)
    
    if(bias){
      delta$basis = get_basis(yDisc,pctVar=.95,B=mvcData$obsBasis$D)
      delta$GPs = vector(mode='list',length=delta$basis$nPC)
      
      for(k in 1:delta$basis$nPC){
        d = darg(d=list(mle=T),X=mvcData$XTdata$obs$X$trans)
        g = garg(g=list(mle=T),y=delta$basis$Vt[k,])
        # Our responses have changed wrt to inputs, so we should not use d=1 anymore
        delta$GPs[[k]] <- newGPsep(X=mvcData$XTdata$obs$X$trans,
                                   Z=delta$basis$Vt[k,],
                                   d=d$start, g=g$start, dK=TRUE)
        cmle <- jmleGPsep(delta$GPs[[k]],
                          drange=c(d$min, d$max),
                          grange=c(g$min, g$max),
                          dab=d$ab,
                          gab=g$ab)
        pred = predGPsep(delta$GPs[[k]], XX = mvcData$XTdata$obs$X$trans, lite=T)
        delta$wMean = rbind(delta$wMean, pred$mean)
        delta$wVar = rbind(delta$wVar, pred$s2)
        deleteGPsep(delta$GPs[[k]])
      }
      # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
      delta$wSamp = rmvnorm(1,as.numeric(delta$wMean),diag(as.numeric(delta$wVar)))
      dim(delta$wSamp) = dim(delta$wMean)
      # update yDisc
      yDisc = yDisc - w_to_y(delta$wSamp,delta$basis$B)
    } 
    ssq.hat = mean(yDisc^2)
    
    # compute simple llh - Sigma_w and Sigma_v are accounted for in calc of yDisc
    ldetS = nY*log(ssq.hat)
    Sinv = (1/ssq.hat)*eye(nY)
    ll = sum(apply(yDisc,2,function(d) -.5*(ldetS + d%*%Sinv%*%d))) - log(ssq.hat) + sum(dbeta(theta,2,2,log=T))
    
    return( list(ll=ll,ssq=ssq.hat) )
  }
  ptm = proc.time()
  t = tInit
  ssq = ll_mcmc(t)$ssq
  llProp = list()
  accept = matrix(0,nrow=1,ncol=pT)
  t.store = matrix(nrow=nsamples,ncol=mvcData$XTdata$pT);t.store[1,] = t
  ssq.store = numeric(nsamples);ssq.store[1] = ssq

  for(i in 2:nsamples){
    if(verbose){
      if(mod(i,nsamples/10)==0){cat(100*(i/nsamples),'% ')}
    }

    for(j in 1:pT){
      tProp = t # initialize tProp to t so that the parameter we are not updating is fixed at its previous value
      tProp[j] = t[j] + (prop.step[j] * runif(1,-0.5, 0.5))
      if(tProp[j]<0 | tProp[j]>1){
        llProp$ll = -Inf
      } else{
        llProp = ll_mcmc(tProp)
      }
      llt = ll_mcmc(t)
      # 4. Accept or reject t
      if(log(runif(1)) < (llProp$ll - llt$ll)){
        t[j] = tProp[j]
        ssq = llProp$ssq
        accept[j] = accept[j] + 1
      }
      t.store[i,j] = t[j]
    }
    # store ssq after both updates have been made
    ssq.store[i] = ssq
  }
  mcmc.time = ptm-proc.time()
  return(list(t.samp=t.store[(nburn+1):nsamples,],ssq.samp=ssq.store[nburn+1:nsamples],
              acpt.ratio=accept/nsamples,
              time=mcmc.time))
}

tune_step_sizes = function(mvcData,nburn,nlevels,target_acc=NULL){
  # t.default = .5
  # step.default = .2
  # step.sizes = matrix(nrow=nlevels,ncol=pT)
  # ex = seq(-(nlevels - 1)/2, (nlevels - 1)/2, length.out=nlevels)
  # for(i in 1:nlevels){
  #   for(j in 1:pT){
  #     if(ex[i] <= 0){
  #       base = 2.0
  #     } else{
  #       base = 20^(2/(nlevels-1))
  #     }
  #     step.sizes[i,j] = step.default * base^ex[i]
  #   }
  # }
  step.sizes = matrix(logseq(.01,2,n=nlevels),nrow=nlevels,ncol=pT,byrow = F)
  
  acc = matrix(nrow=nlevels,ncol=pT)
  for(i in 1:nlevels){
    acc[i,] = mcmc(mvcData,tInit = rep(.5,pT),bias=mvcData$bias,nsamples = nburn,prop.step = step.sizes[i,])$acpt.ratio
  }

  # Compute GLM for each parameter
  step.tuned = numeric(pT)
  if(is.null(target_acc)){
    target_logit = rep(log(1 / (exp(1) - 1)),pT)
  } else{
    target_logit = logit(target_acc)
  }
  for(j in 1:pT){
      y = cbind(acc[,j]*nburn,nburn-(acc[,j]*nburn))
      x = cbind(log(step.sizes[,j]))
      glm_model = glm(y~x,family=binomial(link='logit'))
      coefs = as.numeric(coef(glm_model))
      step.tuned[j] = exp((target_logit[j]-coefs[1])/coefs[2])
  }
  return(step.tuned)
}
