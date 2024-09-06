source('./main.R')
source("./generation.R")
library(RLRsim)
library(lme4)
library(MASS)
library(foreach)
library(MFPCA)

options(warn = -1)

nSubs = 80 # the number of subjects
N = 50 # the number of visits for each subject

## Here only shows the case of p=2, if one wants to consider different cases,
# please change the @param argvals, @param M, @param eFunType, and @ param ignoreDeg
# correspondingly.
# the number of functional predictors
p = 2 

## the range of argument @param t
argvals = list(list(seq(0,1,0.01)), list(seq(0,1,0.01))) 

## the number of true FPCs
M = list(2, 2)

## set eigenfunctions.
eFunType = list( "PolyHigh", "PolyHigh")
ignoreDeg = list(1:2, 1)

 
## set the second functional random effect to be insignificant.
# Determine whether the functional random effect corresponding 
# to the second functional predictor is significant.
# TRUE: non-significance, FALSE: significance.
is.dummy = TRUE 

## with measurement error.
# Determine whether the functional predictors are contaminated by 
# measurement errors. err = 0.5 with measurement errors, 
# err = 0 without measurement errors
err = 0.5  

## generate responses with Binomial distribution
fam = 'binomial'
adj = 1.0 # default setting for binomial.

## generate responses with Poisson distribution
# fam = 'poisson'
# adj = 0.20  # adjust factor for Poisson distribution



## generate pseudo-sample.
DataComb = generate_data_fun(argvals=argvals, M=M, p=2, 
                             eFunType=eFunType, ignoreDeg=ignoreDeg,
                             nSub=nSubs, N=N, is.dummy = is.dummy, err=err,
                             adj = adj, family = fam)

pcaRes = pca.score.fun(XFunData = DataComb$XfunData, p=p)
mComb = pcaRes$mComb
random.fun.data = do.call(cbind, pcaRes$UKcomb)
scores = pcaRes$Ucomb
Bscore = scores*DataComb$B

## prepare data and formula for further analysis.
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


## testing based on the candidate method
gflmmPQL.res =  try(gflmmPQL(fixed = fixed.formula,
                             random.scalar = random.scal.formula,
                             random.function = random.fun.formula,
                             family = fam, data = test.data),
                    silent=T)

test.PQL.res = test.fun(gflmm.res = gflmmPQL.res, mComb = mComb, N = N, nSubs = nSubs, p=p)


## testing based on the proposed method
gflmmLA.wei.res = try(gflmmLA.wei(fixed = fixed.formula,
                                  random.scalar = random.scal.formula,
                                  random.function = random.fun.formula,
                                  family = fam, data =test.data),
                      silent=T)
test.LA.res = test.fun(gflmm.res = gflmmLA.wei.res, mComb = mComb, N = N, nSubs = nSubs, p=p)

## print testing results
data.test = rbind(Candidate = test.PQL.res, Proposed = test.LA.res)
df = data.frame(data.test)
names(df) = c('beta1(t)', 'beta2(t)')
message( paste('Hypothesis testing, Default: Reject for beta1(t) (P<0.05) and ', 
             ifelse(is.dummy, 'accept for beta2(t) (P >= 0.05)', 'reject for beta2(t) (P< 0.05).'), 
             sep=''))
print(df)
