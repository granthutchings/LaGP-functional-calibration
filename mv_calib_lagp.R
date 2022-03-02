library(Matrix)
library(tidyverse)
library(laGP)
library(tgp)
library(parallel)
library(pracma)

# FUNCTION: ll_beta: log likelihood of beta(a,b) evaluated at x. Used as a prior on calibration parameters
ll_beta = function(x,a,b){
  (a-1)*log(x) + (b-1)*log(1-x) -lbeta(a,b)
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
mv.em.pred = function(theta,Xpred,estLS,XTdata,simVt){
  # transform and scale X
  nPC = nrow(simVt)
  XpredTrans = transform_xt(XTdata$sim$X$orig,XTdata$sim$T$orig,Xpred,Tobs=t(replicate(nrow(Xpred),theta,simplify='matrix')))
  XpredSC = get_SC_inputs(estLS,XpredTrans,nPC)
  emGPs = aGPsep_SC_mv(X=XpredSC$XTsim,
                       Z=lapply(1:nPC,function(i) simVt[i,]),
                       XX=lapply(1:nPC,function(i) cbind(XpredSC$Xobs[[i]],XpredSC$Tobs[[i]])),
                       returnYmean = F)
  return(list(GPs=emGPs,XpredSC=XpredSC))
}
yPred.nobias= function(emGPs,ssq.hat,B,BtBinv,nsamples,Yobs,ysamp=T){
  n = ncol(Yobs$orig)
  nY = nrow(Yobs$orig)
  YSamp = array(dim=c(nY,nsamples,n))
  if(ysamp){
    for(i in 1:n){
      YSamp[,,i] = t(mvtnorm::rmvnorm(nsamples,mean=B%*%emGPs$wMean[,i,drop=F],
                                      sigma = ssq.hat*eye(nY) + B%*%diag(emGPs$wVar[,i])%*%t(B))) * Yobs$sd + Yobs$mean
    }
  } else{
    for(i in 1:n){
      vSamp = t(mvtnorm::rmvnorm(nsamples,mean=emGPs$wMean[,i],
                                 sigma = ssq.hat*BtBinv + diag(emGPs$wVar[,i])))
      YSamp[,,i] = w_to_y(vSamp,B,Yobs$mean,Yobs$sd)
    }
  }
  YMean = apply(YSamp,c(1,3),mean)
  YVar = apply(YSamp,c(1,3),var)
  YInt = mclapply(1:n,function(i) t(apply(YSamp[,,i],1,quantile,c(.025,.975))))
  #return(list(vSamp=vSamp,YSamp=YSamp,YMean=YMean,YVar=YVar,YInt=YInt))
  return(list(YSamp=YSamp,YMean=YMean,YVar=YVar,YInt=YInt))
}

yPred.bias = function(theta,Xpred,estLS,XTdata,simVt,
                      dfit,Bobs,nsamples=1){
  nY = nrow(Bobs)
  n = nrow(Xpred)
  emSamp = array(dim=c(nY,nsamples,n))
  discSamp = array(dim=c(nY,nsamples,n))
  #YSamp = array(dim=c(nY,nsamples,n))
  emulator = mv.em.pred(theta,Xpred,estLS,XTdata,simVt)
  discrep = lapply(1:dfit$disc$basis$nPC, function(i) predGPsep(dfit$disc$GPs[[i]],Xpred))
  # make wMean and wVar for discrep
  wMean = NULL; wVar = NULL
  for(i in 1:dfit$disc$basis$nPC){
    wMean = rbind(wMean,discrep[[i]]$mean)
    wVar = rbind(wVar,diag(discrep[[i]]$Sigma))
  }
  
  for(j in 1:n){
    emSamp[,,j] = t(mvtnorm::rmvnorm(nsamples,mean=Bobs%*%emulator$GPs$wMean[,j,drop=F],
                                     sigma = Bobs%*%diag(emulator$GPs$wVar[,j])%*%t(Bobs))) *
                  Ydata$obs$sd + Ydata$obs$mean 
    discSamp[,,j] = t(mvtnorm::rmvnorm(nsamples,mean=dfit$disc$basis$B%*%wMean[,j,drop=F],
                                       sigma = dfit$disc$basis$B%*%diag(wVar[,j],nrow=length(wVar[,j]))%*%t(dfit$disc$basis$B) *
                                         dfit$ssq.hat*eye(nY) )) * Ydata$obs$sd
  }
  YSamp = emSamp + discSamp
  
  YMean = apply(YSamp,c(1,3),mean)
  YVar = apply(YSamp,c(1,3),var)
  YInt = mclapply(1:n,function(i) t(apply(YSamp[,,i],1,quantile,c(.025,.975))))
  emMean = apply(emSamp,c(1,3),mean)
  discMean = apply(discSamp,c(1,3),mean)
  return(list(YSamp=YSamp,YMean=YMean,YVar=YVar,YInt=YInt,
              emMean = emMean, discMean=discMean))
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
                        returnYmean=T, yMean = NULL, ySD = NULL, Bobs=NULL){
  nPC = length(Z)
  nY = nrow(XX[[1]])

  wMean = array(dim=c(nPC,nY))
  wVar = array(dim=c(nPC,nY))
  lagp_fit = vector(mode='list',length=nPC)
  for (i in 1:nPC) {
    lagp_fit[[i]] = aGPsep(X = X[[i]],
                           Z = Z[[i]],
                           XX = XX[[i]],
                           d = list(mle = FALSE, start = 1),
                           g = g,
                           start = start,
                           end = min(end,length(Z[[i]])-1),
                           method = 'nn',
                           verb=0)
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
get_basis = function(Y, nPC = 1, pctVar = NULL){

  retList = list()
  
  # Compute SVD
  nExp = ncol(Y)
  svdY = svd(Y)
  retList$B = svdY$u %*% diag(svdY$d) / sqrt(nExp)
  retList$Vt = t(svdY$v) * sqrt(nExp)

  # Return only the required number of basis vectors
  percentVar = zapsmall((svdY$d) / sum(svdY$d))
  if(!is.null(pctVar)){
    # percentage of variance explained specified
    nPC = which.max(cumsum(percentVar)>pctVar)
  }
  retList$nPC = nPC
  retList$pctVar = percentVar[1:retList$nPC]
  retList$B = as.matrix(retList$B[,1:retList$nPC,drop=F])
  retList$Vt = as.matrix(retList$Vt[1:retList$nPC,,drop=F])
  
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
  # THIS CALCULATION HAS CAUSED ME PROBLEMS
  # If emulation isn't working, check that obs basis is correct
  Bridge = 1e-6 * diag(rep(1,nPC)) # in case Bprod is ill-conditioned
  if(!is.null(sigY)){
    lamY = (1/sigY)*eye(dim(obsBasis$B)[1])
  } else{
    lamY = eye(dim(obsBasis$B)[1])
  }
  Bprod = t(obsBasis$B)%*%lamY%*%obsBasis$B
  obsBasis$Vt = solve(Bprod + Bridge) %*% t(obsBasis$B)%*%lamY%*%Yobs
  #obsBasis$Vt = pinv(obsBasis$B)%*%Yobs
  return(obsBasis)
}

# FUNCTION: mv_lengthscales
# TO DO: this function needs to be updated for only using a subset of the data
### Computes estimated length scale parameters for use in stretching and compressing the inputs space
# Parameters:
# X (numeric) : untransformed matrix of inputs
# Vt (numeric) : matrix 
mv_lengthscales = function(XTdata,Vt,g=1e-7){

  XT = cbind(XTdata$sim$X$trans,XTdata$sim$T$trans)
  pX = XTdata$pX
  pT = XTdata$pT
  nPC = nrow(Vt)
  #XT = t(t(XT))
  
  dConfig = darg(list(mle = TRUE, max = 100), XT)
  # can't instantiate GP's with mclapply because the GP's need unique integer identifiers
  # mclapply will try to give them the same identifier because they aren't communicating
  GPlist = vector(mode='list',length=nPC)
  for(i in 1:nPC){
    GPlist[[i]] = newGPsep(X = XT,
                           Z = Vt[i, ],
                           d = rep(dConfig$start, ncol(XT)),
                           g = g,
                           dK = TRUE)
  }
  
  estLenscales = mclapply(1:nPC, function(i) mleGPsep(GPlist[[i]],param='d',
                                                        tmin = dConfig$min, 
                                                        tmax = dConfig$max,
                                                        ab=dConfig$ab,
                                                        maxit=200)$d)
  for(i in nPC){
    deleteGPsep(GPlist[[i]])
  }
  
  return = list()
  return$XT = estLenscales
  return$X = mclapply(1:nPC, function(i) estLenscales[[i]][1:pX])
  return$T = mclapply(1:nPC, function(i) estLenscales[[i]][(pX+1):(pX+pT)])
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

transform_y = function(Ysim,YindSim,Yobs=NULL,YindObs=NULL,center=T,scale=T, scaletype='columnwise'){
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
    if(scaletype=='columnwise'){
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
          if(scaletype=='columnwise'){
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
  return$Xsim = mclapply(1:nPC, function(i) sc_inputs(XTdata$sim$X$trans,estLS$X[[i]]))
  return$Tsim =  mclapply(1:nPC, function(i) sc_inputs(XTdata$sim$T$trans,estLS$T[[i]]))
  return$XTsim = mclapply(1:nPC, function(i) cbind(return$X[[i]],return$T[[i]]))
  
  if(!is.null(XTdata$obs$X)){
    return$Xobs = mclapply(1:nPC, function(i) sc_inputs(XTdata$obs$X$trans,estLS$X[[i]]))
  }
  if(!is.null(XTdata$obs$T)){
    return$Tobs = mclapply(1:nPC, function(i) sc_inputs(XTdata$obs$T$trans,estLS$T[[i]]))
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

mv.calib.nobias = function(tCalib,SCinputs,estLS,Ydata,simVt,obsBasis,BtBinv,thetaHat=F){
  nPC = obsBasis$nPC
  n = ncol(Ydata$obs$orig)
  tCalibOrig = t(t(tCalib) * XTdata$sim$T$range + XTdata$sim$T$min)
  
  get_ll_ssq = function(theta){
    # transform theta by estimated length-scales
    thetaSC = mclapply(1:nPC, function(j) sc_inputs(theta,estLS$T[[j]]))
    # make prediction matrix from Xobs and theta
    if(SCinputs$pT>1){
      Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(replicate(n,thetaSC[[j]],simplify='matrix'))))
    } else{
      Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(t(replicate(n,thetaSC[[j]],simplify='matrix')))))
    }
    
    
    # Predict from emulator at [Xobs,theta]
    emGPs = aGPsep_SC_mv(X=SCinputs$XTsim,
                         Z=lapply(1:nPC,function(j) simVt[j,]),
                         XX=Xpred,
                         Bobs = obsBasis$B)
  
    # residuals
    wDisc = obsBasis$Vt - emGPs$wMean
    yDisc = Ydata$obs$trans - w_to_y(emGPs$wMean,obsBasis$B)
    yDiscNative = Ydata$obs$orig - w_to_y(emGPs$wMean,obsBasis$B,mean=Ydata$obs$mean,sd=Ydata$obs$sd)
    # error variance estimate
    ssq.hat = mean(yDisc^2)
    ssq.hat.native = mean(yDiscNative^2)
    # likelihood
    ll = numeric(n)
    for(i in 1:n){
      ll[i] = mvtnorm::dmvnorm(x=wDisc[,i],sigma=ssq.hat*BtBinv+diag(emGPs$wVar[,i]),log=T)
    }
    ll = sum(ll) #+ n*sum(ll_beta(theta,2,2))
    # return vector with first element likelihood and second element error variance estimate
    return(c(ll,ssq.hat,ssq.hat.native))
  }
  ll_df = matrix(nrow=nrow(tCalib),ncol=3+SCinputs$pT)
  ll_df = data.frame(cbind(matrix(unlist(mclapply(1:nrow(tCalib), function(i) get_ll_ssq(tCalib[i,,drop=F]))),ncol=3,byrow = T),tCalib))
  colnames(ll_df) = c('ll','ssq.std','ssq.native',paste0("theta",1:SCinputs$pT))
  if(thetaHat){
    theta.hat.id = which.max(ll_df$ll)
    theta.hat = as.matrix(tCalib[theta.hat.id,,drop=F],nrow=1)
    theta.hat.orig = as.matrix(tCalibOrig[theta.hat.id,,drop=F],nrow=1)
    ssq.hat = ll_df$ssq.std[theta.hat.id]
    ssq.hat.orig = ll_df$ssq.native[theta.hat.id]
    return(list(ll_df=ll_df,theta.hat=theta.hat,theta.hat.orig=theta.hat.orig,
                ssq.hat=ssq.hat,ssq.hat.orig=ssq.hat.orig))
  } else{
    return(ll_df)
  }
}

mv.calib.nobias_yllh = function(tCalib,SCinputs,estLS,Ydata,simVt,obsBasis,thetaHat=F){
  nPC = obsBasis$nPC
  n = ncol(Ydata$obs$orig)
  nY = nrow(Ydata$obs$orig)
  tCalibOrig = t(t(tCalib) * XTdata$sim$T$range + XTdata$sim$T$min)
  
  get_ll_ssq = function(theta){
    # transform theta by estimated length-scales
    thetaSC = mclapply(1:nPC, function(j) sc_inputs(theta,estLS$T[[j]]))
    # make prediction matrix from Xobs and theta
    Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(replicate(n,thetaSC[[j]],simplify='matrix'))))
    
    # Predict from emulator at [Xobs,theta]
    emGPs = aGPsep_SC_mv(X=SCinputs$XTsim,
                         Z=lapply(1:nPC,function(j) simVt[j,]),
                         XX=Xpred,
                         Bobs = obsBasis$B)
    
    # residuals
    wDisc = obsBasis$Vt - emGPs$wMean
    yDisc = Ydata$obs$trans - w_to_y(emGPs$wMean,obsBasis$B)
    yDiscNative = Ydata$obs$orig - w_to_y(emGPs$wMean,obsBasis$B,mean=Ydata$obs$mean,sd=Ydata$obs$sd)
    # error variance estimate
    ssq.hat = mean(yDisc^2)
    ssq.hat.native = mean(yDiscNative^2)
    #ssq.hat = apply(yDisc^2,2,mean)
    #ssq.hat.native = apply(yDiscNative^2,2,mean)
    # likelihood
    ll = numeric(n)
    for(i in 1:n){
      ll[i] = mvtnorm::dmvnorm(x=yDisc[,i],sigma=ssq.hat*eye(nY) + obsBasis$B%*%diag(emGPs$wVar[,i])%*%t(obsBasis$B),log=T)
      #ll[i] = mvtnorm::dmvnorm(x=yDisc[,i],sigma=ssq.hat[i]*eye(nY) + obsBasis$B%*%diag(emGPs$wVar[,i])%*%t(obsBasis$B),log=T)
      
    }
    ll = sum(ll) #+ n*sum(ll_beta(theta,2,2))
    # return vector with first element likelihood and second element error variance estimate
    return(c(ll,mean(ssq.hat),mean(ssq.hat.native)))
  }
  ll_df = matrix(nrow=nrow(tCalib),ncol=3+SCinputs$pT)
  ll_df = data.frame(cbind(matrix(unlist(mclapply(1:nrow(tCalib), function(i) get_ll_ssq(tCalib[i,]))),ncol=3,byrow = T),tCalib))
  colnames(ll_df) = c('ll','ssq.std','ssq.native',paste0("theta",1:pT))
  if(thetaHat){
    theta.hat.id = which.max(ll_df$ll)
    theta.hat = as.matrix(tCalib[theta.hat.id,,drop=F],nrow=1)
    theta.hat.orig = as.matrix(tCalibOrig[theta.hat.id,,drop=F],nrow=1)
    ssq.hat = ll_df$ssq.std[theta.hat.id]
    ssq.hat.orig = ll_df$ssq.native[theta.hat.id]
    return(list(ll_df=ll_df,theta.hat=theta.hat,theta.hat.orig=theta.hat.orig,
                ssq.hat=ssq.hat,ssq.hat.orig=ssq.hat.orig))
  } else{
    return(ll_df)
  }
}

mv.calib.bias = function(tCalib,SCinputs,estLS,Ydata,simVt,obsBasis,thetaHat=F){
  discBasis=NULL; emGPs=NULL; w = NULL; discGPs = NULL
  pT = XTdata$pT
  pX = XTdata$pX
  nPC = simBasis$nPC
  n = ncol(Ydata$obs$orig)
  tCalibOrig = t(t(tCalib) * XTdata$sim$T$range + XTdata$sim$T$min)

  get_disc = function(theta){
    thetaSC = lapply(1:nPC, function(j) sc_inputs(theta,estLS$T[[j]]))
    Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(replicate(n,thetaSC[[j]],simplify='matrix'))))
    emGPs = aGPsep_SC_mv(X=SCinputs$XTsim,
                         Z=lapply(1:nPC,function(j) simVt[j,]), XX=Xpred,
                         Bobs = obsBasis$B)
    
    Wdisc = obsBasis$Vt - emGPs$wMean
    Ydisc = Ydata$obs$trans - w_to_y(emGPs$wMean,obsBasis$B)
    YdiscNative = Ydata$obs$orig - w_to_y(emGPs$wMean,obsBasis$B,Ydata$obs$mean,Ydata$obs$sd)
    return(list(emGPs=emGPs,Wdisc=Wdisc,Ydisc=Ydisc,YdiscNative=YdiscNative))
  }
  fit_disc = function(theta,bias,clean){
    disc = NULL; disc$GPs = vector(mode='list',length=nPC)
    disc$basis = get_basis(bias$Ydisc,pctVar = .95)
    nY = nrow(bias$Ydisc)
    
    for(k in 1:disc$basis$nPC){
      d = darg(d=NULL,X=Xobs)
      g = garg(g=NULL,y=disc$basis$Vt[k,])
      # Our responses have changed wrt to inputs, so we should not use d=1 anymore
      disc$GPs[[k]] <- newGPsep(X=Xobs,
                                Z=disc$basis$Vt[k,],
                                d=d$start, g=g$start, dK=TRUE)
      cmle <- jmleGPsep(disc$GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      pred = predGPsep(disc$GPs[[k]], XX = Xobs, lite=T)
      disc$wMean = rbind(disc$wMean, pred$mean)
      # disc$wSigma[[k]] = pred$Sigma
      disc$wVar = rbind(disc$wVar, pred$s2)
      if(clean){
        deleteGPsep(disc$GPs[[k]])
      }
      # disc$GPs[[k]] = aGPsep(X=SCinputs$Xobs[[k]],
      #                        Z=disc$basis$Vt[k,],
      #                        XX=SCinputs$Xobs[[k]],
      #                        start=1,
      #                        end=min(50,length(disc$basis$Vt[k,])-1),
      #                        d = d,
      #                        g = list(start=g$start,mle=T),
      #                        verb=0)
      # disc$wMean = rbind(disc$wMean, disc$GPs[[k]]$mean)
      # disc$wVar = rbind(disc$wVar, disc$GPs[[k]]$var)
    }
    # new residuals (Y-emulator) - discrepancy estimate, this is mean zero
    newDisc = bias$Ydisc - w_to_y(disc$wMean,disc$basis$B)
    ssq.hat = mean(newDisc^2)
    newDiscNative = bias$YdiscNative - w_to_y(disc$wMean,disc$basis$B,sd=Ydata$obs$sd)
    ssq.hat.native = mean(newDiscNative^2)
    ll = 0
    for(i in 1:n){
      ll = ll + mvtnorm::dmvnorm(x=newDisc[,i],sigma=obsBasis$B%*%diag(bias$emGPs$wVar[,i])%*%t(obsBasis$B) +
                                   disc$basis$B%*%diag(disc$wVar[,i],nrow=disc$basis$nPC)%*%t(disc$basis$B) + 
                                   ssq.hat*eye(nY),log = T)
    }
    return(list(disc=disc,ll=ll,ssq.hat=ssq.hat,ssq.hat.native=ssq.hat.native))
  }
  
  ll = numeric(nrow(tCalib))
  for(i in 1:length(ll)){
    bias = get_disc(tCalib[i,,drop=F])
    ll[i] = fit_disc(tCalib[i,,drop=F],bias,clean=T)$ll
  }

  ll_df = as.data.frame(matrix(cbind(ll,tCalibOrig,as.matrix(tCalib)),nrow=length(ll),ncol=2*pT+1))
  names = c("ll")
  names = c(names,paste0("theta",1:pT),paste0("thetaScaled",1:pT))
  colnames(ll_df) = names
  
  if(thetaHat){
    theta.hat = as.matrix(tCalib[which.max(ll_df$ll),,drop=F],nrow=1)
    theta.hat.orig = as.matrix(tCalibOrig[which.max(ll_df$ll),,drop=F],nrow=1)
    bias = get_disc(theta.hat)
    dfit = fit_disc(theta.hat,bias,clean=F)
    return(list(ll_df=ll_df,theta.hat=theta.hat,theta.hat.orig=theta.hat.orig,bias=bias,dfit=dfit))
  } else{
    return(ll_df)
  }
}
