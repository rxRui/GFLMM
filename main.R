pca.score.fun = function(XFunData, p, M = 6, pve=0.90){
  #Step 1: perform MFPCA to generate U and Uk.
  mComb = NULL
  pcaSepRes = vector('list', p)
  UKcomb = vector('list', p)
  xFunDimTemp = NULL
  lambda.values = vector('list', p+1)
  for(i in 1:p){
    XFunDataTemp = XFunData[[i]]
    xFunDim = length(dim(XFunDataTemp@X))
    xFunDimTemp = c(xFunDimTemp, xFunDim)
    if(xFunDim == 3){
      pcaSepRes[[i]] = MFPCA(mFData = multiFunData(XFunDataTemp), M=M,
                             uniExpansions = list(list(type = "splines2D", 
                                                       k = c(8,8))))
      lambda.values[[i+1]] = pcaSepRes[[i]]$values
      M1max = sum(cumsum(pcaSepRes[[i]]$values)/sum(pcaSepRes[[i]]$values) <= pve)+1
      
      UKcomb[[i]] = pcaSepRes[[i]]$scores[,1:M1max]
      mComb = c(mComb, M1max)
      
    } else if (xFunDim == 2){
      pcaSepRes[[i]] = MFPCA(mFData = multiFunData(XFunDataTemp), M=M,
                             uniExpansions = list(list(type = "uFPCA")))
      lambda.values[[i+1]] = pcaSepRes[[i]]$values
      M2max = sum(cumsum(pcaSepRes[[i]]$values)/sum(pcaSepRes[[i]]$values) <= pve)+1
      
      UKcomb[[i]] = pcaSepRes[[i]]$scores[,1:M2max]
      mComb = c(mComb, M2max)
    }else{
      stop('only allowed dimenion no larger than 2')
    }
  }
  if(var(xFunDimTemp)){
    if(p>2)
      stop('For multidimensional cases, the number of functional predictors $p$ is only allowed 2')
    pcaCombRes = MFPCA(mFData = XFunData, M=M,
                       uniExpansions = list(list(type = "splines2D", 
                                                 k = c(8,8)),
                                            list(type = "splines1D", k = 10)))
    Mmax = sum(cumsum(pcaCombRes$values)/sum(pcaCombRes$values) <= pve)+1
    lambda.values[[1]] = pcaCombRes$values
    
  }else{
    # when all is one-dimensional.
    pcaCombRes = MFPCA(mFData = XFunData, M=M,
                       uniExpansions = rep(list(list(type = "splines1D", k = 10)), p))
    Mmax = sum(cumsum(pcaCombRes$values)/sum(pcaCombRes$values) <= pve)+1
    lambda.values[[1]] = pcaCombRes$values
  }
  Ucomb = pcaCombRes$scores[,1:Mmax]
  mComb = c(Mmax, mComb)
  
  return(list(Ucomb = Ucomb, UKcomb = UKcomb, mComb=mComb, 
              pcaCombRes=pcaCombRes, pcaSepRes = pcaSepRes,
              lambda.values = lambda.values))
}

std.gflmm = function(fit.gflmm){
  mcall.unstd = fit.gflmm$mcall
  fit.unstd = fit.gflmm$fit
  unstd.data = eval(mcall.unstd$data)
  random.scalar = fit.gflmm$random.scalar
  random.function = fit.gflmm$random.function
  
  rand.fun.names = (all.vars(random.function))[-1]
  rand.scal.names = (all.vars(random.scalar))[-1]
  rand.names = c(rand.scal.names, rand.fun.names)
  
  W.half = diag(sqrt(unstd.data$weights)) 
  Y.std =  W.half %*% unstd.data$zz
  
  # fixed effects
  X.design = data.frame(model.matrix(fit.unstd, data=unstd.data)) 
  X.std = data.frame(W.half %*% as.matrix(X.design)) 
  X.std.names = names(X.std)
  fixed.function = as.formula(paste("y ~ 0 + ",paste(X.std.names, collapse="+"))) 
  fixed.std = data.frame(y=Y.std, X.std)
  fixed.names = c('y', X.std.names)
  
  # random effects
  group.names = names(fit.unstd@flist) # grouping variables
  extract.groups = data.frame(unstd.data[, (names(unstd.data) %in% group.names)])
  names(extract.groups) = group.names
  rand.vars.names = rand.names[!(rand.names %in% group.names)]
  rand.fun.names = rand.fun.names[!(rand.fun.names %in% group.names)]
  rand.scal.names = rand.scal.names[!(rand.scal.names %in% group.names)]
  
  std.rand.all = sapply(rand.vars.names, function(obj){
    as.matrix(unstd.data[, obj])*diag(W.half)
  })
  std.data.all = cbind(fixed.std, std.rand.all, extract.groups)
  names(std.data.all) = c(fixed.names, rand.vars.names, group.names)
  formula.update = update(fixed.function, random.function)
  if(!is.null(random.scalar)) formula.update = update(formula.update, random.scalar)
  
  mcall.std = mcall.unstd
  mcall.std$formula = formula.update
  mcall.std$data = std.data.all
  mcall.std$weights =  rep(1, nrow(unstd.data))
  
  
  return(list(mcall.std = mcall.std, 
              random.names = rand.vars.names,
              rand.fun.names = rand.fun.names,
              rand.scal.names = rand.scal.names,
              group.names = group.names,
              fixed.names = fixed.names))
}

svd.transform.fun = function(std.res, mComb, nSubs, nRepComb, k=NULL){
  mcall.nonsvd= mcall.std = std.res$mcall.std # save the standardized mcall.
  if(length(mComb)<1)
    stop('This method is built for testing functional random effects.')
  if(is.null(k))
    stop('The k has to be given!')
  
  std.random.fun = mcall.std$data[, std.res$rand.fun.names]
  std.random.scal = mcall.std$data[, std.res$rand.scal.names]
  std.fixed = mcall.std$data[, std.res$fixed.names]
  nonsvd.std = cbind(std.fixed, std.random.scal)
  nonsvd.std.col = ncol(nonsvd.std)
  group.data = mcall.std$data[, std.res$group.names]
  
  svd.random.fun = NULL
  svd.std.transf = NULL
  svd.random.fun.names = NULL
  To.test = NULL
  ID.uni = NULL
  idx.col = if(k == 1){1:mComb[1]}else{
    (sum(mComb[1:(k-1)])+1):sum(mComb[1:k]) }
  
  ID.uni.index = matrix(1:(nSubs*mComb[k]), nrow=mComb[k], ncol=nSubs)
  
  for(i in seq_len(nSubs)){
    ID.uni = c(ID.uni, ID.uni.index[,i], rep(0, nRepComb[i]-mComb[k]))
    
    idx.row = if(i == 1){ 1:nRepComb[1]}else{
      (sum(nRepComb[1:(i-1)])+1):(sum(nRepComb[1:i]))}
    
    svd.random = std.random.fun[idx.row, idx.col]
    nonsvd.random = std.random.fun[idx.row, -idx.col]
    
    svd.res = svd(as.matrix(svd.random)%*%t(as.matrix(svd.random)))
    svd.u.tra = t(svd.res$v)
    svd.d = (svd.res$d)[1:mComb[k]]
    
    To.test = c(To.test, sqrt(svd.d), rep(0, nRepComb[i]-mComb[k]))
    
    svd.std.transf.temp = svd.u.tra%*%as.matrix(nonsvd.std[idx.row, ])
    svd.std.transf = rbind(svd.std.transf, svd.std.transf.temp)
    if((ncol(nonsvd.random)) != 0){
      svd.random.fun.temp = svd.u.tra%*%as.matrix(nonsvd.random)
      svd.random.fun = rbind(svd.random.fun, svd.random.fun.temp)
    }
  }
  if((ncol(nonsvd.random)) != 0){svd.random.fun.names = std.res$rand.fun.names[-idx.col]}
  
  transf.data = data.frame(cbind(svd.std.transf, 
                                 svd.random.fun, 
                                 To.test = To.test, 
                                 ID.uni = ID.uni, 
                                 group.data))
  
  names(transf.data) = c(std.res$fixed.names, 
                         std.res$rand.scal.names,
                         svd.random.fun.names,
                         'To.test', 'ID.uni', 
                         std.res$group.names)
  
  fixed.formula = paste('y~0+', 
                        paste0(std.res$fixed.names[-1], 
                               collapse = '+'))
  if(length(std.res$rand.scal.names) != 0 ){
    random.formula = paste('.~. + (0 + To.test | ID.uni)',' + ',
                           paste('(0 +', std.res$rand.scal.names,'|', 
                                 std.res$group.names, ')', 
                                 sep=' ',collapse = '+'))
  }else{
    random.formula = paste('.~. + (0 + To.test | ID.uni)')
  }
  
  if((ncol(nonsvd.random)) != 0){
    random.formula = paste(random.formula, '+',
                           paste('(0 +', svd.random.fun.names,'|', 
                                 std.res$group.names, ')', sep=' ', collapse = '+'))
  }
  
  random.names = all.vars(as.formula(random.formula))
  random.names = random.names[!(random.names %in%c('ID.uni', std.res$group.names))][-1]
  
  svd.formula = update(as.formula(fixed.formula),
                       as.formula(random.formula))
  mcall.svd = mcall.nonsvd # update mcall.std
  mcall.svd$data = transf.data
  mcall.svd$formula = svd.formula
  
  
  
  return(list(mcall.svd = mcall.svd,
              random.names = random.names,
              fixed.names = std.res$fixed.names,
              random.formula = random.formula,
              fixed.formula = fixed.formula
  ))
}


test.fun = function(gflmm.res, mComb, nSubs, N, p, nRepComb=NULL, tol=1e-4, verbose=0L, is.joint = FALSE){
  if((length(mComb)>p) & !is.joint) mComb = mComb[-1]
  if(is.null(nRepComb))
    nRepComb = rep(N, nSubs)
  
  if(is.joint) 
    p = 1 # for joint testing only!!!
  if( (class(gflmm.res$fit)[[1]] == 'lmerMod') || (class(gflmm.res$fit)[[1]] == 'glmerMod')){
    std.res = std.gflmm(gflmm.res)
    test.res = lapply(seq_len(p), function(k){
      svd.std.res = svd.transform.fun(std.res = std.res,
                                      mComb = mComb,
                                      nSubs = nSubs,
                                      nRepComb = nRepComb, k = k)
      test.prob = fun.test.aRLRT(svd.std.res = svd.std.res)
      return(test.prob$aRLRT$p.value)
    })
  }else
    stop('only glmerMod or glm is accepted')
  
  return(unlist(test.res))
}

### Derived from Chen et al. (2019) glmmVCtest
gflmmPQL = function (fixed, random.function, random.scalar=NULL, data, family,
                     control= NULL, weights=NULL,
                     REML = TRUE, niter = 50, verbose = 0L){
  if (!requireNamespace("lme4", quietly = TRUE)){
    stop("package 'lme4' is essential")}
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if(is.null(control)) control = lme4::lmerControl(check.conv.singular = 'ignore')
  
  data.call = mcall = match.call()
  if(is.null(random.function)) stop("functional random effects are required")
  formula = update(fixed, random.function)
  if(!is.null(random.scalar)) formula = update(formula, random.scalar)
  formula.data = formula
  
  # construct updated dataset to new formula
  data.name = names(data.call)[-1L]
  keep.data.name = is.element(data.name, c('data', 'weights', 'subset'))
  for(i in data.name[!keep.data.name]) data.call[[i]] = NULL
  allvars = all.vars(formula.data)
  data.call$formula = as.formula(paste("~", paste(allvars, collapse = "+")))
  environment(data.call$formula) = environment(fixed)
  data.call$drop.unused.levels = TRUE
  data.call[[1L]] = quote(stats::model.frame)
  data.update = eval.parent(data.call)
  # get start value for further modeling by means of glm().
  offset = model.offset(data.update)
  if (is.null(offset)) offset = 0
  
  weights = model.weights(data.update)
  if (is.null(weights)) weights = rep(1, nrow(data.update))
  
  data.update$weights = weights
  fit.glm = glm(formula = fixed, family = family, data = data.update, 
                weights = weights)
  
  weights.prior = fit.glm$prior.weights
  eta = fit.glm$linear.predictors
  zz = eta + fit.glm$residuals - offset
  weights = fit.glm$weights
  
  mcall.name = names(mcall)[-1L]
  keep.formula.name = is.element(mcall.name, c('formula', 'data', 'REML', 
                                               'subset', 'weights',
                                               'offset'))
  for(i in mcall.name[!keep.formula.name]) mcall[[i]] = NULL
  mcall[[1L]] = quote(lme4::lmer)
  formula.update = update(formula, zz~.)
  mcall[["formula"]] = formula.update
  mcall$REML = REML 
  if(!is.null(control)) mcall$control = control
  data.update$zz = zz
  
  data.update$weights = sqrt(weights)  
  mcall$weights = quote(weights)
  mcall$data = quote(data.update)
  
  for (i in seq_len(niter)) {
    if (verbose)
      message(gettextf("iteration %d", i), domain = NA)
    etaold = eta
    fit = eval(mcall)
    eta = fitted(fit) + offset
    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2))
      break
    mu = family$linkinv(eta)
    mu.eta.val = family$mu.eta(eta)
    data.update$zz = eta + (fit.glm$y - mu)/mu.eta.val - offset
    weights = weights.prior * mu.eta.val^2/family$variance(mu)
    data.update$weights = sqrt(weights)
    mcall$weights = quote(weights)
    mcall$data = quote(data.update)
  }
  
  data.update$weights = sqrt(weights)
  mcall$data = data.update
  
  return(list(fit = fit, fit0=fit.glm,
              mcall = mcall, 
              random.scalar=random.scalar, 
              random.function=random.function))
}


gflmmLA.wei = function (fixed, random.scalar, random.function, data, family,
                        control= NULL, weights=NULL, offset = NULL, REML = TRUE, 
                        niter = 50, verbose = 0L){
  if (!requireNamespace("lme4", quietly = TRUE)){
    print('here')
    stop("package 'lme4' is essential")}
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if(is.null(control)) 
    control = lme4::lmerControl(check.conv.singular = 'ignore')
  
  data.call = mcall = match.call()
  if(is.null(random.function)) stop('functional random effect are required')
  formula = update(fixed, random.function)
  if(!is.null(random.scalar)) formula = update(formula, random.scalar)
  formula.data = formula
  
  # construct updated dataset to new formula
  data.name = names(data.call)[-1L]
  keep.data.name = is.element(data.name, c('data', 'subset', 'weights', 'offset'))
  for(i in data.name[!keep.data.name]) data.call[[i]] = NULL
  allvars = all.vars(formula.data)
  data.call$formula = as.formula(paste("~", paste(allvars, collapse = "+")))
  environment(data.call$formula) = environment(fixed)
  data.call$drop.unused.levels = TRUE
  data.call[[1L]] = quote(stats::model.frame)
  data.update = eval.parent(data.call)
  
  # get start value for further modeling by means of glm().
  offset = model.offset(data.update)
  if (is.null(offset)) offset = rep(0, nrow(data.update))
  data.update$offset = offset
  weights = model.weights(data.update)
  if (is.null(weights)) weights = rep(1, nrow(data.update))
  data.update$weights = weights
  
  glm.res = glm(formula = fixed, data = data.update, family = family, weights = weights)
  
  # random initiation of theta
  random.names = c(all.vars(random.scalar), all.vars(random.function))
  len.rand.names = length(random.names)
  num.con = sum(random.names == '.')
  initial.start = list(theta=rep(1, len.rand.names-2*num.con),
                       fixef=glm.res$coefficients)
  weights.prior = glm.res$prior.weights
  weights = 1/sqrt(weights.prior) 
  data.update$weights = weights
  glmer.res = glmer(formula = formula, data = data.update,
                    family = family, 
                    start = initial.start,
                    weights = weights,
                    control = glmerControl(check.conv.singular = 'ignore',
                                           optCtrl=list(maxfun = 1e+8)))
  
  weights.out = weights(glmer.res, type='working')
  mu = fitted(glmer.res)
  eta = family$linkfun(mu) + offset
  mu.eta.val = family$mu.eta(eta)
  data.update$zz = eta + (glm.res$y - mu)/mu.eta.val - offset
  data.update$weights = (weights.out)
  
  mcall.name = names(mcall)[-1L]
  keep.formula.name = is.element(mcall.name, c('formula', 'data', 'REML', 
                                               'subset', 'weights',
                                               'offset'))
  for(i in mcall.name[!keep.formula.name]) mcall[[i]] = NULL
  mcall[[1L]] = quote(lme4::lmer)
  formula.update = update(formula, zz~.)
  mcall[["formula"]] = formula.update
  mcall$REML = REML 
  if(!is.null(control)) mcall$control = control
  
  ### Ready for tuning
  data.update$weights = rep(1, nrow(data.update))
  mcall$weights = quote(weights)
  mcall$data = quote(data.update)
  fit = eval(mcall)
  etaold = eta
  
  dat = data.update[, names(fixef(fit))]
  eta = as.vector(as.matrix(dat)%*%fixef(fit))
  eta = eta/sigma(fit)
  
  mu = family$linkinv(eta)
  mu.eta.val = family$mu.eta(eta)
  weights.out = weights.prior * mu.eta.val^2/family$variance(mu)
  data.update.old = data.update
  data.update$zz = eta + (glm.res$y-mu)/mu.eta.val - offset
  data.update$weights = sqrt(weights.out)
  mcall$data = quote(data.update)
  mcall$weights = quote(weights)
  for( i in seq_len(niter)){
    if (verbose)
      message(gettextf("iteration %d", i), domain = NA)
    weights.out.old = weights.out
    fitold = fit
    etaold = eta
    fit = try(eval(mcall), silent = T)
    if(class(fit) =="try-error"){
      message('Warning: return the previous fit')
      break
    }
    dat = data.update[, names(fixef(fit))]
    eta = as.vector(as.matrix(dat)%*%fixef(fit))
    if (sum((eta - etaold)^2) < 1e-6*sum((eta)^2))
      break
    mu = family$linkinv(eta)
    mu.eta.val = family$mu.eta(eta)
    weights.out = weights.prior * mu.eta.val^2/family$variance(mu)
    data.update.old = data.update
    data.update$zz = (eta + (glm.res$y-mu)/mu.eta.val - offset)
    data.update$weights = sqrt(weights.out)
  }
  
  if(class(fit) =="try-error"){
    mcall$data = quote(data.update.old)
    mcall$weights = quote(weights)
    fit = eval(mcall)
    mcall$data = data.update.old
    mcall$weights = sqrt(weights.out.old)
  }else{
    data.update$weights = sqrt(weights.out)
    mcall$data = data.update
  }
  
  return(list(fit = fit,
              mcall = mcall, 
              random.scalar=random.scalar, 
              random.function=random.function))
}

### Derived from Chen et al. (2019) glmmVCtest
fun.test.aRLRT = function(svd.std.res){
  mcall.svd.origin = mcall.svd = svd.std.res$mcall.svd
  svd.data = mcall.svd$data
  fixed.names = svd.std.res$fixed.names
  fixed.formula = as.formula(paste(paste(fixed.names[1], '~'), 
                                   paste0(fixed.names[-1], collapse ='+')))
  
  fit.alt.std = try(eval(mcall.svd), silent=T)
  if('try-error' %in% class(fit.alt.std)){
    stop('Error in lmer model estimation under alternative. Consider simplifying
             or rescaling variables.')
  }
  
  if(length(svd.std.res$random.names)==1){ # no nuisance random effects
    mcall.null = mcall.svd # update from alt model fit
    null.formula = update(mcall.svd$formula, '.~. + (0 + X.Intercept.| ID.uni) - (0 + To.test | ID.uni)') #only nuis
    mcall.null$formula = null.formula
    fit.null.std = try(eval(mcall.null), silent = T) # fixed effects only
    if('try-error' %in% class(fit.null.std)){
      stop('Error in lmer model estimation under null. Consider simplifying
               or rescaling variables.')
    }
    fit.test.std = fit.alt.std # same as alternative model
  } else { # nuisance r.effect
    mcall.null = mcall.test = mcall.svd # update from alt model fit
    null.formula = update(mcall.svd$formula, '.~. - (0 + To.test | ID.uni)') #only nuis
    mcall.null$formula = null.formula
    
    test.formula = update(fixed.formula, '.~. + 0 + (0 + To.test | ID.uni)') # only test
    mcall.test$formula = test.formula
    
    fit.null.std = try(eval(mcall.null),silent=T)
    if('try-error' %in% class(fit.null.std)){
      stop('Error in lmer model estimation under null. Consider simplifying
                 or rescaling variables.')
    }
    fit.test.std = try(eval(mcall.test),silent=T)
    if('try-error' %in% class(fit.test.std)){
      stop('Error in lmer model estimation for testing variable. Consider rescaling.')
    }
  }
  # testing
  aRLRT = RLRsim::exactRLRT(m=fit.test.std, mA=fit.alt.std, m0=fit.null.std)
  return(list(aRLRT=aRLRT, fit.alt=fit.alt.std, fit.null=fit.null.std, fit.test=fit.test.std))
}
