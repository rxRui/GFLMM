setwd("./code/")
source('main.R')
source("data-generation.R")
library(RLRsim)
library(lme4)
library(nlme)
library(MASS)
library(foreach)
library(funData)
library(MFPCA)

# #mute for "Observed RLRT statistic is 0, no simulation performed." 
# # message("Observed RLRT statistic is 0, no simulation performed.")
# trace(exactRLRT, edit = T) 
options(warn = -1)

nSubs = 80
N = 50
p = 2
argvals = list(list(seq(0,1,0.01)),
               list(seq(0,1,0.01)))
M = list(2, 2)
eFunType = list( "PolyHigh", "PolyHigh")
ignoreDeg = list(1:2, 1)
DataComb = generate_data_fun(argvals=argvals, M=M, p=2, 
                             eFunType=eFunType, ignoreDeg=ignoreDeg,
                             nSub=nSubs, N=N, is.dummy = TRUE, sdk=0.5, 
                             err=0.5)

pcaRes = pca.score.fun(XFunData = DataComb$XfunData, p=p)
mComb = pcaRes$mComb
random.fun.data = do.call(cbind, pcaRes$UKcomb)
scores = pcaRes$Ucomb
Bscore = scores*DataComb$B

test.data = data.frame(1, DataComb$Y, DataComb$Z,
                       DataComb$B, Bscore, DataComb$subID,
                       scores, random.fun.data)
names(test.data) = c('inter', 'y', 'Z', 'B', paste('bscore', 1:mComb[1], sep=''),
                     'subID', paste('score', 1:mComb[1], sep=''),
                     paste('f1.score', 1:mComb[2], sep=''),
                     paste('f2.score', 1:mComb[3], sep=''))

fixed.formula = as.formula(paste('y ~ 0 +inter+ Z + B', paste('score', 1:mComb[1], sep='', collapse = '+'),
                                 paste('bscore', 1:mComb[1], sep='', collapse = '+'),
                                 sep = '+'))
random.fun.formula = as.formula(paste('~.',paste('(0 + f1.score', 1:mComb[2], '| subID)',
                                                 sep='', collapse = ' +'),
                                      paste('(0 + f2.score', 1:mComb[3], '| subID)', 
                                            sep='', collapse = ' +'), sep = '+'))

random.scal.formula = as.formula(paste('~.', '(0 + inter|subID)', '(0 + B|subID)', sep='+'))



gflmmPQL.res =  try(gflmmPQL(fixed = fixed.formula,
                             random.scalar = random.scal.formula,
                             random.function = random.fun.formula,
                             verbose = T,
                             family = binomial, data = test.data),
                    silent=T)

test.PQL.res = test.fun(gflmm.res = gflmmPQL.res, mComb = mComb, N = N, nSubs = nSubs, p=p)
print(test.PQL.res)

gflmmLA.wei.res = try(gflmmLA.wei(fixed = fixed.formula,
                                  random.scalar = random.scal.formula,
                                  random.function = random.fun.formula,
                                  family = binomial, data =test.data, 
                                  verbose = T),
                      silent=T)
test.LA.res = test.fun(gflmm.res = gflmmLA.wei.res, mComb = mComb, N = N, nSubs = nSubs, p=p)
print(test.LA.res)

