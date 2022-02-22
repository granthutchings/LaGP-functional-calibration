library(Matrix)
library(tidyverse)
library(laGP)
library(tgp)
library(parallel)
library(pracma)

ll_beta = function(x,a,b){
  (a-1)*log(x) + (b-1)*log(1-x) -lbeta(a,b)
}
# FUNCTION: multivariate aGPsep for use with stretched and compressed inputs only
## Parameters:
# X : list of length nPC which contains stretched and compressed training input matrices
# Z : list of length nPC which contains training response vectors
# XX: list of length nPC with contains stretched and compressed prediction input matrices

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
yPred.nobias= function(wMean,wVar,ssq.hat,B,BtBinv,nsamples,yMean,ySD){
  n = ncol(wMean)
  YSamp = array(dim=c(nrow(YindObs),nsamples,n))
  for(j in 1:n){
    vSamp = t(mvtnorm::rmvnorm(nsamples,mean=wMean[,j],
                                       sigma = ssq.hat*BtBinv + diag(wVar[,j])))
    YSamp[,,j] = w_to_y(vSamp,B,yMean,ySD)
  }
  YMean = apply(YSamp,c(1,3),mean)
  YVar = apply(YSamp,c(1,3),var)
  YInt = mclapply(1:n,function(i) t(apply(YSamp[,,i],1,quantile,c(.025,.975))))
  return(list(vSamp=vSamp,YSamp=YSamp,YMean=YMean,YVar=YVar,YInt=YInt))
}

aGPsep_SC_mv = function(X, Z, XX,
                        start=6, end=50, 
                        g=1/10000, method="nn",
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
                           method = method,
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
### Accepts a matrix of response values Y, and returns the basis decomposition from SVD
## Parameters:
# Y (numeric)       : matrix of response values to decompose into basis representation. Y should be
#                     formatted such that each row represents the multivariate response for a single experiment
# nPC (int)        : number of principal components desired from basis decomposition
# pctVar (numeric) : desired percentage of variability explained by PC's. Overrides nPC if specified
## Returns: List containings the following
# B (numeric)       : matrix of nPC basis vectors for Ytrans
# Vt (numeric)      : matrix of nPC response vectors for Ytrans
# nPC (integer)    : number of PC's used for the basis decomposition
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
get_obs_basis = function(simBasis,Yobs,YindSim,YindObs,sigY=NULL){
  obsBasis = list()
  obsBasis$nPC = simBasis$nPC
  if(!isTRUE(all.equal(YindObs,YindSim))){
    B = matrix(nrow=nrow(YindObs),ncol=nPC)
    if(ncol(YindSim)==1){
      # 1d interpolation
      for(i in 1:simBasis$nPC){
        B[,i] = interp1(as.numeric(YindSim),as.numeric(simBasis$B[,i]),as.numeric(YindObs))
      }
    } else if(ncol(YindSim)==2){
      # 2d interpolation
      y = unique(YindSim[,1]); x = unique(YindSim[,2])
      yp = YindObs[,1]; xp = YindObs[,2]
      for(i in 1:simBasis$nPC){
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

transform_y = function(Ysim,YindSim,Yobs=NULL,YindObs=NULL,center=T,scale=T, scaletype='scalar'){
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

mv.calib.nobias = function(thetaCalib,SCinputs,estLS,Ydata,simVt,obsBasis,BtBinv,clean=T,obsReps=1){
  nPC = obsBasis$nPC
  n = ncol(Ydata$obs$orig)

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
    # likelihood
    ll = numeric(n)
    for(i in 1:n){
      ll[i] = mvtnorm::dmvnorm(x=wDisc[,i],sigma=ssq.hat*BtBinv+diag(emGPs$wVar[,i]),log=T)
    }
    ll = sum(ll) + n*sum(ll_beta(theta,2,2))
    # return vector with first element likelihood and second element error variance estimate
    return(c(ll,ssq.hat,ssq.hat.native))
  }
  ll_df = matrix(nrow=nrow(thetaCalib),ncol=3+SCinputs$pT)
  ll_df = data.frame(cbind(matrix(unlist(mclapply(1:nrow(thetaCalib), function(i) get_ll_ssq(thetaCalib[i,]))),ncol=3,byrow = T),thetaCalib))
  colnames(ll_df) = c('ll','ssq.std','ssq.native',paste0("theta",1:pT))
  return(ll_df)
}

mv.calib.biasDtD = function(thetaCalib,SCinputs,estLS,Ydata,simVt,obsBasis,BtBinv,clean=T){
  nPC = simBasis$nPC
  n = ncol(Ydata$obs$orig)
  
  get_disc = function(theta){
    thetaSC = lapply(1:nPC, function(j) sc_inputs(theta,estLS$T[[j]]))
    Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(replicate(n*obsReps,thetaSC[[j]],simplify='matrix')))[1:n,])
    emGPs = aGPsep_SC_mv(X=SCinputs$XTsim,
                         Z=lapply(1:nPC,function(j) simVt[j,]), XX=Xpred,
                         Bobs = obsBasis$B)
    
    Wdisc = obsBasis$Vt - do.call(cbind, replicate(obsReps, emGPs$wMean, simplify=FALSE))
    Ydisc = Ydata$obs$trans - do.call(cbind, replicate(obsReps, emGPs$yMean, simplify=FALSE))
    ssq.hat = mean(Ydisc^2)
    return(list(emGPs=emGPs,Wdisc=Wdisc,Ydisc=Ydisc,ssq.hat=ssq.hat))
  }
  fit_disc = function(theta,bias,clean){
    discGPs = NULL
    # discBasis = get_basis(bias$Ydisc,nPC)
    discBasis = get_obs_basis(simBasis,bias$Ydisc,YindSim,YindObs)
    DtDinv = solve(t(discBasis$B)%*%discBasis$B)
    
    wMean = NULL; wSigma = NULL
    ll = numeric(nPC)
    for(k in 1:nPC){
      d = darg(d=NULL,X=SCinputs$Xobs[[k]])
      g = garg(g=NULL,y=discBasis$Vt[k,])
      # Our responses have changed wrt to inputs, so we should not use d=1 anymore
      discGPs$GPs[[k]] <- newGPsep(X=SCinputs$Xobs[[k]],
                                   Z=discBasis$Vt[k,],
                                   d=d$start, g=g$start, dK=TRUE)
      cmle <- jmleGPsep(discGPs$GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      ll[k] = llikGPsep(discGPs$GPs[[k]])
      pred = predGPsep(discGPs$GPs[[k]], XX = SCinputs$Xobs[[k]])
      wMean = rbind(wMean, pred$mean)
      wSigma[[k]] = pred$Sigma
      if(clean){
        deleteGPsep(discGPs$GPs[[k]])
      }
    }
    # # new residuals
    # wDisc = bias$Wdisc - wMean
    # # calculate log likelihood
    # llh = numeric(n)
    # for(i in 1:n){
    #   llh[i] = mvtnorm::dmvnorm(x=wDisc[,i],sigma=bias$ssq.hat*BtBinv+diag(bias$emGPs$wVar[,i]+diag(wVar[,i])),log=T)
    # }
    return(sum(ll) + nPC*sum(ll_beta(theta,2,2)))
  }
  
  ll = numeric(nrow(thetaCalib))
  for(i in 1:length(ll)){
    bias = get_disc(thetaCalib[i,])
    ll[i] = fit_disc(thetaCalib[i,],bias,clean=T)
  }
  
  ll_df = as.data.frame(matrix(cbind(ll,tCalibOrig,as.matrix(thetaCalib),ssq.hat),nrow=length(ll),ncol=2*pT+2))
  names = c("ll")
  names = c(names,paste0("theta",1:pT),paste0("thetaScaled",1:pT),'ssq')
  colnames(ll_df) = names
  return(ll_df)
}
mv.calib.bias = function(thetaCalib,SCinputs,estLS,Ydata,simVt,obsBasis,BtBinv,clean=T){
  discBasis=NULL; emGPs=NULL; w = NULL; discGPs = NULL
  set.seed(seed)
  pT = XTdata$pT
  pX = XTdata$pX
  nPC = simBasis$nPC
  n = ncol(Ydata$obs$orig)/obsReps

  #tCalibOrig = t(t(thetaCalib) * XTdata$sim$T$range + XTdata$sim$T$min)

  covMat = mclapply(1:nPC, function(k) exp(-0.5 * fields::rdist(SCinputs$Xobs[[k]])^2) + diag(1e-8, nrow(SCinputs$Xobs[[k]])))
  S = as.matrix(Matrix::bdiag(covMat))
  Sinv = chol2inv(chol(S))
  ldetS = determinant(S)$modulus

  # covMat = rep(list(BtBinv),n*obsReps)
  # covMat = as.matrix(Matrix::bdiag(covMat))
  # ldetB = n*obsReps*determinant(BtBinv)$modulus

  get_disc = function(theta){
    thetaSC = lapply(1:nPC, function(j) sc_inputs(theta,estLS$T[[j]]))
    Xpred = lapply(1:nPC, function(j) cbind(SCinputs$Xobs[[j]],t(replicate(n*obsReps,thetaSC[[j]],simplify='matrix')))[1:n,])
    emGPs = aGPsep_SC_mv(X=SCinputs$XTsim,
                         Z=lapply(1:nPC,function(j) simVt[j,]), XX=Xpred,
                         Bobs = obsBasis$B)

    Wdisc = obsBasis$Vt - do.call(cbind, replicate(obsReps, emGPs$wMean, simplify=FALSE))
    Ydisc = Ydata$obs$trans - do.call(cbind, replicate(obsReps, emGPs$yMean, simplify=FALSE))
    ssq.hat = mean(Ydisc^2)
    return(list(Wdisc=Wdisc,Ydisc=Ydisc,ssq.hat=ssq.hat))
  }
  fit_disc = function(theta,Ydisc,clean){
    discGPs = NULL
    discBasis = get_basis(Ydisc,nPC)
    wMean = NULL; wVar = NULL
    for(k in 1:nPC){
      d = darg(d=NULL,X=SCinputs$Xobs[[k]])
      g = garg(g=NULL,y=discBasis$Vt[k,])
      # Our responses have changed wrt to inputs, so we should not use d=1 anymore
      discGPs$GPs[[k]] <- newGPsep(X=SCinputs$Xobs[[k]],
                                   Z=discBasis$Vt[k,],
                                   d=d$start, g=g$start, dK=TRUE)
      cmle <- jmleGPsep(discGPs$GPs[[k]],
                        drange=c(d$min, d$max),
                        grange=c(g$min, g$max),
                        dab=d$ab,
                        gab=g$ab)
      pred = predGPsep(discGPs$GPs[[k]], XX = SCinputs$Xobs[[k]])
      wMean = rbind(wMean, pred$mean)
      wVar = rbind(wVar, pred$var)
      }
      if(clean){
        deleteGPsep(discGPs$GPs[[k]])
      }
    wVec = as.numeric(t(wMean))
    VtDiscVec <- as.numeric(t(discBasis$Vt))
    # calculate log likelihood
    diff = wVec - VtDiscVec
    N = length(diff)
    RSS = t(diff)%*%Sinv%*%diff
    return(-.5*(N*log(.5*RSS) + ldetS) + sum(ll_beta(theta,2,2)))
  }

  ll = numeric(nrow(thetaCalib))
  for(i in 1:length(ll)){
    Ydisc = get_disc(thetaCalib[i,])$Ydisc
    ll[i] = fit_disc(thetaCalib[i,],Ydisc,clean=T)
  }

  ll_df = as.data.frame(matrix(cbind(ll,tCalibOrig,as.matrix(thetaCalib),ssq.hat),nrow=length(ll),ncol=2*pT+2))
  names = c("ll")
  names = c(names,paste0("theta",1:pT),paste0("thetaScaled",1:pT),'ssq')
  colnames(ll_df) = names
  return(ll_df)
}
