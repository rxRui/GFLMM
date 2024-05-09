# GFLMM
Here is the .R files to reproduce the results for the GFLMM and the testing method proposed by Rui et al. (2024). 
For more details, please refer to the paper "Unveiling Schizophrenia: A Study with Generalized Functional Linear Mixed Model via the Investigation of Functional Random Effects".
Before running the code, please make sure that the following required packages are installed, all .R files are placed at current work space. Specifically, packages include:
library(MASS)
library(foreach)
library(MFPCA)
library(lme4)
library(RLRsim)

# NOTE: The testing process for the real data is time-consuming due to the large size of data, especially because of the large number of time points (301). 
# It may take several hours depending on the used computer. Please wait. If one only concerned with the code reproducibility of our proposed method, or 
# the testing result on a pseudo-data, just run the Demo.R file directly.  
