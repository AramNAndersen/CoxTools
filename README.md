# CoxTools
Tools for clinical forecasting from drug screens using Cox regression.
The package requires dplyr, reshape, doSNOW, glmnet and glmnetUtils.

# Installation:
``` r
devtools::install_github("AramNAndersen/CoxTools",
                         ref="main",
                         auth_token = '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')
``` 
# Example:
``` r
library(glmnet)
library(glmnetUtils)
library(reshape)
library(dplyr)
library(doSNOW)
library("CoxTools")

# Model C-index testing
list_Cox_testing <- Cox_forecasting_glmnet_CVA(X_data=X.train.AUC, 
                                               y_data=y.train, 
                                               alpha=c(0,1,"CVA"), 
                                               lambda=c(exp(seq(-4,6, 0.1))),
                                               free_cores = 2,
                                               test.n= c(6,4), 
                                               nfolds = nrow(y.train.red),
                                               iter=200,
                                               log_AUC=1:2,
                                               Patient.Z=1:2,
                                               Drug.Z =1:2,
                                               RCPC=0:4)
                                                              
# Variable importance/predictivity scoring by drug withdrawal (with mode reduction)
Results_Cox_drug_naive <- Cox_forecasting_drug_withdrawal(X.train.AUC, 
                                                          y.train, 
                                                          Reduce="Naive",
                                                          alpha=0, 
                                                          lambda=c(exp(seq(-4,6, 0.1))),
                                                          free_cores = 2,
                                                          test.n= c(6,4), 
                                                          iter=200,
                                                          log_AUC=1,
                                                          Patient.Z=1,
                                                          Drug.Z =2,
                                                          RCPC=0)

# Variable survival association statistics
Results_Cox_betas <- Cox_bootstrapping(X.train.AUC, 
                                       y.train, 
                                       alpha=0, 
                                       lambda=c(exp(seq(-4,6, 0.1))),
                                       pre.CV=FALSE,
                                       lambda_opt = 0,
                                       free_cores = 2,
                                       iter=200,
                                       log_AUC=1,
                                       Patient.Z=1,
                                       Drug.Z =2,
                                       RCPC=0)
                                       
# Combination data model C-index testing
list_combined_data_Cox_testing <- Cox_forecasting_glmnet_combination(X.train.AUC, 
                                                                     y.train, 
                                                                     X.train.gen,
                                                                     alpha=0, 
                                                                     lambda=c(exp(seq(-4,6, 0.1))),
                                                                     free_cores = 2,
                                                                     test.n= c(6,4), 
                                                                     iter=200,
                                                                     log_AUC=c(1:2),
                                                                     Patient.Z=c(1:2),
                                                                     Drug.Z =c(1:2),
                                                                     RCPC=c(0,2))
``` 
