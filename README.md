# CoxTools
Development of tools for clinical forecasting from drug screens using regularized Cox regression based on the Glmnet framework.
The package requires dplyr, reshape, doSNOW, glmnet and glmnetUtils.

# Installation:
``` r
devtools::install_github("AramNAndersen/CoxTools-dev",
                         ref="main")
``` 
# Example:
``` r
library(glmnet)
library(glmnetUtils)
library(reshape)
library(dplyr)
library(doSNOW)
library(CoxTools)

# Model C-index testing
list_Cox_testing <- Cox_forecasting_glmnet_CVA(X_data=X.AUC, 
                                               y_data=Y, 
                                               alpha=c(0,1,"CVA"), 
                                               lambda=c(exp(seq(-4,6, 0.1))),
                                               free_cores = 2,
                                               test.n= c(6,4), 
                                               nfolds = nrow(Y),
                                               iter=200,
                                               log_AUC=1:2,
                                               Patient.Z=1:2,
                                               Drug.Z =1:2,
                                               RCPC=0:4)
``` 
