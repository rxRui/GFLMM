# GFLMM
Here are files to reproduce the functional random effects testing method for the GFLMM model proposed by Rui et al. (2024). For more details, please refer to the paper "Unveiling Schizophrenia: A Study with Generalized Functional Linear Mixed Model via the Investigation of Functional Random Effects". 

Before running the code, please make sure that the following required packages are installed and all .R files are placed at current workspace. Specifically, packages include:

library(MASS) 

library(foreach) 

library(MFPCA)

library(lme4)

library(RLRsim)

# main.R 
The main.R includes core functions related to the proposed testing method, which are:

pca.score.fun(): Perform joint and separate FPCA to get related FPCs;

gflmmPQL(): The modified function of glmmPQL.mod.R in glmmVCtest developed by Chen et al. (2019), which is used to get working responses as well as weight matrices;

gflmmLA.wei(): Our proposed method to get working responses and related weight matrices;

std.gflmm(): Obtain normalized working responses;

svd.transform.fun(): Transfer the model with multiple additive random terms to the desired one with single random term; 

fun.test.aRLRT(): Fit models that are used to feed into the exactRLRT() function in the RLRsim package;

test.fun(): Obtain P values corresponding to all functional random effects based on the fitted models from fun.test.aRLRT();

fit.fun(): Fit the data by means of  GFLM and GFLMM by using glmer() function in the lme4 package;

construct_beta_fun(): Extract estimated beta functions from the fit.fun();

get.est.coef.fun(): Extract all estimated terms from GFLM and GFLMM, including scalar-valued coefficients, functional coefficients, etc.

# generation.R
The generation.R includes functions to generate simulation data in which the functionality of each function is shown by its name:

simMultiWeight.mod(); 

generate_score_fun(); 

generate_func_predictor_fun(); 

generate_beta_delta_fun(); 

generate_beta.i_fun(); 

generate_data_fun().

# Demo.R 
Demo.R can be directly executed to obtain testing results based on the generated pseudo-sample, using our proposed method (gflmmLA.wei()) as well as a candidate method (gflmmPQL()). The Poisson-distributed case can also be accessed by simply changing the parameters from `fam=binomial` and `adj=1.0` to `fam=poisson` and `adj=0.25`. 
Here, the `adj` parameter is used to control the size of the Poisson-distributed response for numerical analysis.


# FPCA_comparison.R
The file FPCA_comparison.R is used to illustrate the robust performance of different FPCA methods for parameter estimation due to various representations of multivariate functional predictors $X$. In particular, FPCA methods in the R package \texttt{MFPCA} and \texttt{fgm} are used to produce different types of representations for the same $X$.

The following R packages are required:

library(fgm)

library(lme4)

library(RLRsim)

library(FDA)

library(mvtnorm)

library(MFPCA)

library(fda.usc)
