simMultiWeight.mod = function (argvals, M, eFunType, ignoreDeg = NULL){
  p <- length(argvals)
  dimsSupp <- foreach::foreach(j = seq_len(p), .combine = "c") %do% 
    {
      length(argvals[[j]])
    }
  if (any(dimsSupp > 2)) 
    stop("Function simMultiWeight: method is not implemented for objects of dimension > 2!")
  if (p > 1) {
    if (isTRUE(do.call(all.equal, lapply(M, prod)))) {
      Mtotal <- prod(M[[1]])
    }
    else stop("Function simMultiWeight: basis dimensions must be equal!")
  }
  else {
    Mtotal <- prod(M[[1]])
  }
  basis <- vector("list", p)
  for (j in seq_len(p)) {
    if (dimsSupp[j] == 1) 
      basis[[j]] <- eFun(argvals[[j]][[1]], 
                         M = M[[j]], 
                         ignoreDeg = ignoreDeg[[j]], 
                         type = eFunType[[j]])
    else basis[[j]] <- tensorProduct(eFun(argvals[[j]][[1]], 
                                          M = M[[j]][1], 
                                          ignoreDeg = ignoreDeg[[j]][[1]], 
                                          type = eFunType[[j]][1]), 
                                     eFun(argvals[[j]][[2]], M = M[[j]][2], 
                                          ignoreDeg = ignoreDeg[[j]][[2]], 
                                          type = eFunType[[j]][2]))
  }
  return(multiFunData(basis))
}

# generate scores.
generate_score_fun = function(M.true.max, p, nSubs, N, is.indep = FALSE){
  scoreList = vector('list', p)
  
  etaMIJ = mvrnorm(nSubs*N, mu = rep(0, M.true.max),
                   Sigma = diag(rep(0.5, M.true.max)))
  
  for(j in 1:p){
    if(is.indep){
      rhoMI = 0
    }else{
      rhoMI = apply(mvrnorm(nSubs, mu=rep(0, M.true.max),
                            Sigma = diag(rep(0.1, M.true.max))),
                    2, function(obj){rep(obj, each=N)})
    }
    rhoMIJ = mvrnorm(nSubs*N, mu=rep(0, M.true.max), 
                     Sigma = diag(rep(j, M.true.max)))
    
    scoreList[[j]] = rhoMI + rhoMIJ + etaMIJ
  }
  return(scoreList)
}

# generate functional predictors
generate_func_predictor_fun = function(trueBasis, scoreRes, M.true.max, p, err = 0.){
  XFunDataTemp = vector('list', p)
  for (j in seq_len(p)) {
    Xbtemp = apply(trueBasis[[j]]@X, -1, function(v) { scoreRes[[j]] %*% v })
    dim(Xbtemp) = c(nrow(scoreRes[[j]]), nObsPoints(trueBasis[[j]]))
    XFunDataTemp[[j]] = funData(trueBasis[[j]]@argvals, Xbtemp)
  }
  XFunData = multiFunData(XFunDataTemp)
  if(err != 0)
    XFunData = addError(XFunData, sd=err)
  
  return(list(XFunData = XFunData))
}

# generate beta(t) and delta(t)
generate_beta_delta_fun = function(trueBasis, p, M.true.max){
  coefs = rep(2, M.true.max)
  betaFuns = vector("list", p)
  for (j in seq_len(p)) {
    Beta.temp = apply(trueBasis[[j]]@X, -1, function(obj) {coefs %*% obj})
    dim(Beta.temp) = c(1, nObsPoints(trueBasis[[j]]))
    betaFuns[[j]] = funData(trueBasis[[j]]@argvals, Beta.temp)
  }
  
  betaFuns = multiFunData(betaFuns)
  
  return(list(betaFuns=betaFuns, deltaFuns=-1*betaFuns))
}

# generate funcitonal random effects
generate_beta.i_fun = function(trueBasis, p, nSubs, N, M.true.max, is.dummy=FALSE, sdk=0.5){
  beta.i = vector('list', p)
  for(j in 1:p){
    theta.I = mvrnorm(nSubs, mu=rep(0, M.true.max),
                      Sigma = diag(c(sdk^(j-1), sdk^(j))))
    if(is.dummy){ 
      if(j == p){
        Bi.temp = array(0, dim = c(nSubs, dim(trueBasis[[j]]@X)[-1]))
      }else{
        Bi.temp = apply(trueBasis[[j]]@X, -1, function(obj){theta.I%*%obj})
      }
    }else{
      Bi.temp = apply(trueBasis[[j]]@X, -1, function(obj){theta.I%*%obj})
    }
    beta.i[[j]] = funData(trueBasis[[j]]@argvals, Bi.temp)
  }
  beta.I = if(p > 1){multiFunData(beta.i)}else{beta.i}
  return(beta.I)
}

# generate data function.
generate_data_fun = function(argvals, M, p, eFunType, ignoreDeg=NULL, nSubs=10, N=40,
                             is.dummy = FALSE, sdk=0.5, err=0., is.indep=FALSE, DataComb.fit=NULL){
  trueBasis = simMultiWeight.mod(argvals=argvals, M=M, 
                                 eFunType=eFunType, 
                                 ignoreDeg=ignoreDeg)
  M.true.max = max(unlist(M)) 
  scoreRes = generate_score_fun(M.true.max = M.true.max, p = p, nSubs = nSubs, N=N, is.indep=is.indep)
  XFunDataRes = generate_func_predictor_fun(trueBasis = trueBasis,
                                            scoreRes = scoreRes, M.true.max = M.true.max,
                                            p = p, err = err)
  
  alpha0 = 0.5
  alpha = 2
  gamma = 2
  if(is.null(DataComb.fit)){
    beta.delta.res = generate_beta_delta_fun(trueBasis=trueBasis, 
                                             p=p, M.true.max=M.true.max)
    betaFuns = beta.delta.res$betaFuns
    deltaFuns = beta.delta.res$deltaFuns
    
    
    beta.I = generate_beta.i_fun(trueBasis= trueBasis, p=p, sdk = sdk, 
                                 nSubs = nSubs, is.dummy = is.dummy,
                                 N=N, M.true.max = M.true.max)
    
    gammaI = rep(rnorm(nSubs, mean = gamma, sd=1), each=N)
    alphaI = rep(rnorm(nSubs, sd=1), each=N)
  }else{
    betaFuns = DataComb.fit$true.coef$beta.true.fun
    deltaFuns = DataComb.fit$true.coef$delta.true.fun
    beta.I = DataComb.fit$true.coef$beta.i.true.fun
    gammaI.uni = unique(DataComb.fit$true.coef$gammaI)
    alphaI.uni = unique(DataComb.fit$true.coef$alphaI)
    gammaI = rep(gammaI.uni[1:nSubs], each = N)
    alphaI = rep(alphaI.uni[1:nSubs], each = N)
  }
  
  Z = rep(runif(nSubs, min = 0, max = 3), each=N)
  B = rbinom(N*nSubs, 1, prob = 0.5)
  
  Xbeta = scalarProduct(XFunDataRes$XFunData, betaFuns)
  Xdelta = scalarProduct(XFunDataRes$XFunData, deltaFuns)
  
  XbetaI = as.numeric(sapply(seq_len(nSubs), function(m){
    scalarProduct(extractObs(beta.I, obs = m),
                  extractObs(XFunDataRes$XFunData, obs = (N*(m-1)+1):(N*m)))
  }))
  
  
  XdeltaB = Xdelta*B
  Zalpha = Z*alpha
  BgammaI = B*gammaI
  
  ETA = alpha0 + alphaI + Zalpha + BgammaI + Xbeta + XdeltaB + XbetaI
  prob= 1/(1+exp(-ETA))
  Y = rbinom(N*nSubs, 1, prob)
  subID = rep(1:nSubs, each=N)
  
  DataComb = list(Y=Y, Z=Z, B=B, subID = subID,
                  XfunData = XFunDataRes$XFunData)
  return(DataComb)
}

