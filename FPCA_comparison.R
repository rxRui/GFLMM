library(fgm)
library(lme4)
library(RLRsim)
library(fda)
library(mvtnorm)
library(MFPCA)
library(fda.usc)
options(warn = -1)
## Variables
# Omega - list of precision matrices, one per eigenfunction
# Sigma - list of covariance matrices, one per eigenfunction
# theta - list of functional principal component scores
# phi - list of eigenfunctions densely observed on a time grid
# y - list containing densely observed multivariate (p-dimensional) functional data

n_run = 100
p = 2
nSubs = 100
N = 50
nbasis = 6
time.grid.length = 101
Sigma = matrix(c(0.99, 0.25, 0.25, 0.99), 2, 2)
subID = rep(1:nSubs, each=N)

MISE = NULL
alpha_sigma = NULL
betaFunList1 = NULL
betaFunList2 = NULL


phi.basis=create.fourier.basis(rangeval=c(0,1), nbasis=21, period=1)
t = seq(0, 1, length.out = time.grid.length)
chosen.basis = c(2, 3, 6, 7, 10, 11, 16, 17, 20, 21)
phi.org = t(pblackict(phi.basis, t))
phi = phi.org[chosen.basis[1:nbasis],]

#'####### generate Beta(t)
nbasis.beta = 2
# beta_coef1 = runif(nbasis.beta, min = 0, max=5)
beta_coef1 = c(4, 0.5)
beta1 = as.numeric(beta_coef1%*%phi.org[chosen.basis[1:nbasis.beta],])
# plot(beta1)

beta_coef2 = c(0.5, -6)
beta2 = as.numeric(beta_coef2%*%phi.org[chosen.basis[1:nbasis.beta],])
# plot(beta2)

beta1Fun = funData(X=t(beta1), argvals = t)
beta2Fun = funData(X=t(beta2), argvals = t)

for(i in 1:n_run){
  theta = lapply(1:nbasis, function(b){
    t(rmvnorm(n = nSubs*N, sigma = Sigma^b))} )
  theta.reshaped = lapply( 1:p, function(j){
    t(sapply(1:nbasis, function(i) theta[[i]][j,]))
  })
  
  X = lapply(theta.reshaped, function(th) t(th)%*%phi)
  
  
  
  Xbeta =  apply(X[[1]], 1, function(obj) {
    int.simpson(fdataobj = fdata(mdata = obj*beta1, argvals = t))} ) + 
    apply(X[[2]], 1, function(obj) {
      int.simpson(fdataobj = fdata(mdata = obj*beta2, argvals = t))} )
  
  #'###### generate Beta_i (t)
  nbasis.beta.i = 2
  beta.i.coefs = lapply(1:nbasis.beta.i, function(b){
    t(rmvnorm(n = nSubs, sigma = Sigma))} )
  beta.i.coefs.reshaped = lapply( 1:p, function(j){
    t(sapply(1:nbasis.beta.i, function(i) beta.i.coefs[[i]][j,]))
  })
  
  beta.i = lapply(beta.i.coefs.reshaped, 
                  function(th) {t(th)%*%phi.org[chosen.basis[1:nbasis.beta.i],]})
  
  XbetaI = unlist(lapply(seq_len(nSubs), function(m){
    res = rowSums(do.call('cbind', 
                          lapply(1:p, function(j){
                            X.uni.temp = X[[j]][(N*(m-1)+1):(N*m), ]
                            x.fdata = fdata(mdata = X.uni.temp, argvals = t)
                            beta.i.fdata = fdata(mdata = beta.i[[j]][m, ], argvals = t)
                            int.res = int.simpson(x.fdata*beta.i.fdata)
                            return(int.res)}
                          )
    )) 
  }))
  
  #' generate alpha_i
  alpha.i = rep(rnorm(nSubs, sd=0.5), each=N)
  
  eta = Xbeta + XbetaI + alpha.i
  mu = 1/(1+exp(-eta))
  Y = rbinom(nSubs*N, size=1, prob = mu)
  
  

  ## Solve
  res = pfpca(X)
  mMax1 = min(which(res$FVE >= 0.95 * 100))
  mMax1
  
  score1 = NULL
  score2 = NULL
  for(l in seq_len(mMax1)){
    score1 = rbind(score1, res$theta[[l]][1,])
    score2 = rbind(score2, res$theta[[l]][2,])
  }
  scores1 = cbind(t(score1), t(score2))
  
  #' build glmer model for PFPCA
  test.data1 = data.frame(Y, alpha.i, scores1, subID, scores1)
  names(test.data1) = c( 'y', "I", paste('ascore', 1:mMax1, sep=''), 
                         paste('bscore', 1:mMax1, sep=''), 'subID',
                         paste('f1.score', 1:mMax1, sep=''),
                         paste('f2.score', 1:mMax1, sep=''))
  
  fixed.formula1 = as.formula(paste('y ~ 0', paste('ascore', 1:mMax1, sep='', collapse = '+'),
                                    paste('bscore', 1:mMax1, sep='', collapse = '+'),
                                    sep = '+'))
  random.fun.formula1 = as.formula(paste('~.+ (0 + 1 | subID)', paste('(0 + f1.score', 1:mMax1, '| subID)',
                                                                      sep='', collapse = ' +'),
                                         paste('(0 + f2.score', 1:mMax1, '| subID)', 
                                               sep='', collapse = ' +'), sep = '+'))
  all.formula1 = update(fixed.formula1, random.fun.formula1)
  glmer.model1 = glmer(formula = all.formula1, data = test.data1, 
                       family = binomial)
  
  
  beta1Fun1.mat = fixef(glmer.model1)[1:mMax1]%*%res$phi
  beta2Fun1.mat = fixef(glmer.model1)[-(1:mMax1)]%*%res$phi
  
  beta1Fun1 = funData(argvals = t, X=beta1Fun1.mat)
  beta2Fun1 = funData(argvals = t, X=beta2Fun1.mat)
  # plot(beta1Fun, col='black')
  # plot(beta1Fun1, add=T)
  # plot(beta2Fun, col='black')
  # plot(beta2Fun1, add=T)
  
  betaFun1 = list(beta1Fun1, beta2Fun1)
  betaFunList1 = append(betaFunList1, list(betaFun1))
  
  mise1 = norm((beta1Fun - beta1Fun1), square=T) + norm((beta2Fun - beta2Fun1), square=T) 
  sig1 = sqrt((summary(glmer.model1))$varcor$subID[1])
  
  
  
  # glm.model1 = glm(Y~-1+scores1, family = binomial)
  # 
  # 
  # beta1Fun1.mat = glm.model1$coefficients[1:mMax1]%*%res$phi
  # beta2Fun1.mat = glm.model1$coefficients[-(1:mMax1)]%*%res$phi
  # 
  # beta1Fun1 = funData(argvals = t, X=beta1Fun1.mat)
  # beta2Fun1 = funData(argvals = t, X=beta2Fun1.mat)
  # plot(beta1Fun1)
  # plot(funData(X=t(beta1), argvals = t), add=T, col='black')
  # plot(beta2Fun1)
  # plot(funData(X=t(beta2), argvals = t), add=T, col='black')
  
  
  
  # estimation based on MFPCA
  f1 = funData(argvals = t, X = X[[1]])
  f2 = funData(argvals = t, X = X[[2]])
  
  xFunData = multiFunData(f1, f2)
  # plot(xFunData)
  
  pcaRes = MFPCA(mFData = xFunData, M=20,
                 uniExpansions = list(list(type = "uFPCA"),
                                      list(type = "uFPCA")))
  mMax2 = min(which(cumsum(pcaRes$values)/sum(pcaRes$values) >= 0.95))
  mMax2
  scores2 = (pcaRes$scores[, 1:mMax2])
  
  
  #' build glmer model for MFPCA
  test.data2 = data.frame(Y, scores2, subID, scores2)
  names(test.data2) = c( 'y', paste('ascore', 1:mMax2, sep=''),
                         'subID',
                         paste('f1.score', 1:mMax2, sep=''))
  
  fixed.formula2 = as.formula(paste('y ~ 0', paste('ascore', 1:mMax2, sep='', collapse = '+'),
                                    sep = '+'))
  
  random.fun.formula2 = as.formula(paste('~.+ (0 + 1 | subID)', paste('(0 + f1.score', 1:mMax2, '| subID)',
                                                    sep='', collapse = ' +'), sep = '+'))
  all.formula2 = update(fixed.formula2, random.fun.formula2)
  glmer.model2 = glmer(formula = all.formula2, data = test.data2, family = binomial)
  
  # beta1
  beta1_MFPCA = pcaRes$functions[[1]][1:mMax2]
  beta1.basis = beta1_MFPCA*fixef(glmer.model2)
  beta1Fun2 = funData(X=t(colSums(beta1.basis@X)), argvals = beta1.basis@argvals)
  
  # beta2
  beta2_MFPCA = pcaRes$functions[[2]][1:mMax2]
  beta2.basis = beta2_MFPCA*fixef(glmer.model2)
  beta2Fun2 = funData(X=t(colSums(beta2.basis@X)), argvals = beta2.basis@argvals)
  
  # plot(beta1Fun2)
  # plot(funData(X=t(beta1), argvals = t), add=T, col='black')
  # plot(beta2Fun2)
  # plot(funData(X=t(beta2), argvals = t), add=T, col='black')
  # 
  betaFun2 = list(beta1Fun2, beta2Fun2)
  betaFunList2 = append(betaFunList2, list(betaFun2))
  
  
  mise2 = norm((beta1Fun - beta1Fun2), square=T) + norm((beta2Fun - beta2Fun2), square=T) 
  sig2 = sqrt((summary(glmer.model2))$varcor$subID[1])
  
  
  MISE = rbind(MISE, c(mise1, mise2))
  alpha_sigma = rbind(alpha_sigma, c(sig1, sig2))
  # glm.model2 = glm(Y~-1+scores2, family = binomial)
  # 
  # # beta1
  # beta1_MFPCA = pcaRes$functions[[1]][1:mMax2]
  # 
  # beta1.basis = beta1_MFPCA*glm.model2$coefficients
  # 
  # beta1Fun2 = funData(X=t(colSums(beta1.basis@X)), argvals = beta1.basis@argvals)
  # plot(beta1Fun2)
  # plot(funData(X=t(beta1), argvals = t), add=T, col='black')
  # 
  # # beta2
  # beta2_MFPCA = pcaRes$functions[[2]][1:mMax2]
  # 
  # beta2.basis = beta2_MFPCA*glm.model2$coefficients
  # 
  # beta2Fun2 = funData(X=t(colSums(beta2.basis@X)), argvals = beta2.basis@argvals)
  # plot(beta2Fun2)
  # plot(funData(X=t(beta2), argvals = t), add=T, col='black')
  
}

par(mfrow=c(2,1), mar=c(4,2,1,1))
colnames(MISE) = c('PFPCA', 'MFPCA')
boxplot(MISE, col=c('darkred', 'darkorange'), xlab='MISE')
colnames(alpha_sigma) = c('PFPCA', 'MFPCA')
boxplot(alpha_sigma, col=c('darkred', 'darkorange'), xlab='STD')


par(mfrow=c(2,2))
# PFPCA
coef_beta1Fun1 = NULL
plot(beta1Fun, col='black', lwd=4, xlab='PFPCA')
for(j in 1:n_run){
  plot(betaFunList1[[j]][[1]], add=T, col='gray')
  coef_beta1Fun1 = rbind(coef_beta1Fun1, betaFunList1[[j]][[1]]@X)
}
plot(beta1Fun, col='black', lwd=4, add=T)
lines(t, colMeans(coef_beta1Fun1), col='red', lty=2, lwd=4)

coef_beta2Fun1 = NULL
plot(beta2Fun, col='black', lwd=4, xlab='PFPCA')
for(j in 1:n_run){
  plot(betaFunList1[[j]][[2]], add=T, col='gray')
  coef_beta2Fun1 = rbind(coef_beta2Fun1, betaFunList1[[j]][[2]]@X)
}
plot(beta2Fun, col='black', lwd=4, add=T)
lines(t, colMeans(coef_beta2Fun1), col='red', lty=2, lwd=4)


# MFPCA
coef_beta1Fun2 = NULL
plot(beta1Fun, col='black', lwd=4, xlab="MFPCA")
for(j in 1:n_run){
  plot(betaFunList2[[j]][[1]], lwd=0.5, add=T, col='gray')
  coef_beta1Fun2 = rbind(coef_beta1Fun2, betaFunList2[[j]][[1]]@X)
}
plot(beta1Fun, col='black', lwd=4, add=T)
lines(t, colMeans(coef_beta1Fun2), col='darkorange', lty=2, lwd=4)



coef_beta2Fun2 = NULL
plot(beta2Fun, col='black', lwd=4, xlab="MFPCA")
for(j in 1:n_run){
  plot(betaFunList2[[j]][[2]], add=T, col='gray')
  coef_beta2Fun2 = rbind(coef_beta2Fun2, betaFunList2[[j]][[2]]@X)
  
}
plot(beta2Fun, col='black', lwd=4, add=T)
lines(t, colMeans(coef_beta2Fun2), col='darkorange', lty=2, lwd=4)

