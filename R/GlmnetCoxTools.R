#'
#' Cox forecasting glmnet with alpha cross-validation
#'
#' C-index testing of Cox partial likelihood glmnet models for survival forecasting from ex vivo drug screens. With CVA.
#'
#' @param X_data Input drug screen data
#' @param y_data Survival data
#' @param alpha Penalty type. 0 for L2 (Ridge), 1 for L1 (Lasso), or any number between 0 and 1 for elastic net penalty. Set "CVA" for automatic selection of optimal penalty type.
#' @param lambda Vector of regularization penalties.
#' @param free_cores Number of free cores.
#' @param test.n Vector of test set balance. First integer is the number of survivals, and second integer is the number of deaths.
#' @param nfolds Cross-validation fold.
#' @param iter Number of test iterations.
#' @param log_AUC log2 transform of AUC. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Patient.Z Patient-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Drug.Z Drug-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param RCPC Removal of confounding principal components.
#' @importFrom magrittr %>%
#'
#' @return
#' A list containing results from cross-validation alpha optimization and C-index test results
#'
#'
#' @author
#' Aram N. Andersen \email{aram.n.andersen@@gmail.com}
#'
#' @export
Cox_forecasting_glmnet_CVA <- function(X_data,
                                       y_data,
                                       alpha=0,
                                       lambda=c(exp(seq(-4,6, 0.1))),
                                       free_cores = 2,
                                       test.n= c(6,4),
                                       nfolds = nrow(y.train.red),
                                       iter=200,
                                       log_AUC=c(1:2),
                                       Patient.Z=c(1:2),
                                       Drug.Z =c(1:2),
                                       RCPC=c(0,1,2,3,4)){
  cva_fit_cv_best_alldata <- c()
  cva_fit_plots_alldata <- list()
  cva_fit_coefs.alldata <- c()
  df_C.index.alldata <- c()

  for(al in alpha){
    for(Transform in c("log2(AUC)", "AUC")[log_AUC]){
      for(pt.st in c(TRUE, FALSE)[Patient.Z] ){
        for(cpc in RCPC){
          for(drug.st in c(TRUE, FALSE)[Drug.Z]){
            i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", al)
            set.seed(1)
            X <- X_data
            l=min(data.frame(replace(X, X == 0, 1)))
            if(Transform=="log2(AUC)"){
              if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
            }else{
              X <- data.matrix(X)
            }
            if(pt.st){
              X <- (X - rowMeans(X))/matrixStats::rowSds(X)
            }
            if(cpc>0){
              x_sd <- matrixStats::colSds(X)
              x_mean <- colMeans(X)
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
              svdz <- svd(t(X))
              X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
              X <- X*x_sd + x_mean
              X <- t(X)
              rm(svdz, x_mean, x_sd)
            }

            if(drug.st){
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
            }

            Y <- y_data

            cores <- parallel::detectCores()-free_cores
            cluster.cores<-makeCluster(cores)
            cva <- cva.glmnet(x=X, y=Y, alpha = seq(0, 1, len = 20)^3, family = "cox",
                              standardize = FALSE, type.measure= "deviance", cv.type	="min", nfolds=nfolds,
                              outerParallel=cluster.cores)
            stopCluster(cluster.cores)
            closeAllConnections()
            rm(cluster.cores)
            rm(cores)

            names(cva$modlist) <- round(cva$alpha, digits = 6)
            cva_fit_cv <- cva$modlist %>% lapply(function(x) do.call("cbind", x[c(1:6)]) %>% as_tibble())
            cva_fit_cv <- cva_fit_cv %>%
              bind_rows(.id = "alpha") %>%
              mutate(alpha = as.numeric(alpha))

            # Collect cvm plots drug screen data
            p=cva_fit_cv %>%
              mutate(f_alpha = factor(alpha), cv_min = min(cvm)) %>%
              ggplot(aes(x = lambda, y = cvm, colour = f_alpha)) +
              scale_x_log10() +
              geom_line() +
              geom_hline(aes(yintercept = cv_min)) +
              labs(colour = expression(alpha)) +
              facet_wrap(~f_alpha, scales = "free_y")+
              theme_minimal()+
              theme()+
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            cva_fit_plots_alldata[[i]] <- p

            # Collect top cva results drug screen data
            cva_fit_cv$rank <- rank(cva_fit_cv$cvm)
            cva_fit_cv$class <- "Ridge"
            cva_fit_cv$class[cva_fit_cv$alpha !=0] <- "ElasticNet"
            cva_fit_cv$class[cva_fit_cv$alpha ==1] <- "Lasso"
            cva_fit_cv <- data.frame(cva_fit_cv)
            rownames(cva_fit_cv) <- paste0(cva_fit_cv$class,cva_fit_cv$alpha, cva_fit_cv$rank)
            cva_fit_cv_best <- cva_fit_cv[(cva_fit_cv %>% mutate(id = paste0(class, alpha)) %>% group_by(id) %>% summarise(rank.min = min(rank)) %>% mutate(best.hit = paste0(id, rank.min)))$best.hit,]
            cva_fit_cv_best$dataset.id <- i
            cva_fit_cv_best_alldata <- rbind(cva_fit_cv_best_alldata, cva_fit_cv_best)

            # Collect coefficients for all models (drug screen data)
            coefs.alpha <- c()
            for(j in 1:nrow(cva_fit_cv_best)){
              model <- glmnet(X, Y, family = "cox", alpha = cva_fit_cv_best$alpha[j], standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              coefs  <- data.frame(as.matrix(coef(model, s = cva_fit_cv_best$lambda[j])))
              coefs <- data.matrix(coefs)
              colnames(coefs) <- paste0("alpha", cva_fit_cv_best$alpha[j])
              coefs.alpha <- cbind(coefs.alpha, coefs)
              rm(model, coefs)
            }
            coefs.alpha<-data.frame(coefs.alpha)
            coefs.alpha$dataset.id <- i
            cva_fit_coefs.alldata <- rbind(cva_fit_coefs.alldata, coefs.alpha)

            # Compute out-of-sample prediction c-index (drug screen data)
            if(al == "CVA"){
              a = cva_fit_cv_best$alpha[which.min(cva_fit_cv_best$cvm)]
            }else{
              a = as.numeric(al)
            }
            cores <- parallel::detectCores()-free_cores
            cluster.cores<-makeCluster(cores)
            registerDoSNOW(cluster.cores)
            pb <- txtProgressBar(max=iter, style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
              set.seed(b)
              setTxtProgressBar(pb,b)
              X.0 <- X[Y[,2]==0,]
              X.1 <- X[Y[,2]==1,]
              Y.0 <- Y[Y[,2]==0,]
              Y.1 <- Y[Y[,2]==1,]

              ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
              ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

              train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
              test <- rbind(X.0[ind.0,], X.1[ind.1,])

              Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
              Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

              model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              nfolds=nrow(train)
              model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

              c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
              ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

              rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

              return(c(c,ct))
            }
            stopCluster(cluster.cores)
            closeAllConnections()
            rm(cluster.cores)
            rm(cores)
            gc()

            df_C_index <- do.call(rbind, list_C_index)
            df_C_index <- data.frame(df_C_index)
            colnames(df_C_index) <- c("C_index_test", "C_index_train")
            df_C_index$Iteration <- 1:iter
            df_C_index$ID <- i

            df_C.index.alldata <- rbind(df_C.index.alldata, df_C_index)

            cat("\nAnalysis completed for: ", i,"\n")
            rm(X,Y, cva_fit_cv_best, coefs.alpha, df_C_index, cva, cva_fit_cv, list_C_index, opts, pb, p)

          }
        }
      }
    }
  }
  return(list(CVA_results = cva_fit_cv_best_alldata,
              CVA_plots = cva_fit_plots_alldata,
              CVA_coefficients= cva_fit_coefs.alldata,
              C_index_results = df_C.index.alldata))
}
#'
#'
#' Cox forecasting glmnet
#'
#' C-index testing of Cox partial likelihood glmnet models for survival forecasting from ex vivo drug screens.
#'
#' @param X_data Input drug screen data
#' @param y_data Survival data
#' @param alpha Penalty type. 0 for L2 (Ridge), 1 for L1 (Lasso), or any number between 0 and 1 for elastic net penalty.
#' @param lambda Vector of regularization penalties.
#' @param free_cores Number of free cores.
#' @param test.n Vector of test set balance. First integer is the number of survivals, and second integer is the number of deaths.
#' @param iter Number of test iterations.
#' @param log_AUC log2 transform of AUC. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Patient.Z Patient-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Drug.Z Drug-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param RCPC Removal of confounding principal components.
#' @importFrom magrittr %>%
#'
#' @return
#' A list containing C-index test results
#'
#'
#' @author
#' Aram N. Andersen \email{aram.n.andersen@@gmail.com}
#'
#' @export
Cox_forecasting_glmnet <- function(X_data,
                                   y_data,
                                   alpha=0,
                                   lambda=c(exp(seq(-4,6, 0.1))),
                                   free_cores = 2,
                                   test.n= c(6,4),
                                   nfolds=nrow(y_data),
                                   iter=200,
                                   log_AUC=c(1:2),
                                   Patient.Z=c(1:2),
                                   Drug.Z =c(1:2),
                                   RCPC=c(0,1,2,3,4)){

  df_C.index.alldata <- c()

  for(a in alpha){
    for(Transform in c("log2(AUC)", "AUC")[log_AUC]){
      for(pt.st in c(TRUE, FALSE)[Patient.Z] ){
        for(cpc in RCPC){
          for(drug.st in c(TRUE, FALSE)[Drug.Z]){
            i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
            set.seed(1)
            X <- X_data
            l=min(data.frame(replace(X, X == 0, 1)))
            if(Transform=="log2(AUC)"){
              if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
            }else{
              X <- data.matrix(X)
            }
            if(pt.st){
              X <- (X - rowMeans(X))/matrixStats::rowSds(X)
            }
            if(cpc>0){
              x_sd <- matrixStats::colSds(X)
              x_mean <- colMeans(X)
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
              svdz <- svd(t(X))
              X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
              X <- X*x_sd + x_mean
              X <- t(X)
              rm(svdz, x_mean, x_sd)
            }

            if(drug.st){
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
            }

            Y <- y_data

            cores <- parallel::detectCores()-free_cores
            cluster.cores<-makeCluster(cores)
            registerDoSNOW(cluster.cores)
            pb <- txtProgressBar(max=iter, style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
              set.seed(b)
              setTxtProgressBar(pb,b)
              X.0 <- X[Y[,2]==0,]
              X.1 <- X[Y[,2]==1,]
              Y.0 <- Y[Y[,2]==0,]
              Y.1 <- Y[Y[,2]==1,]

              ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
              ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

              train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
              test <- rbind(X.0[ind.0,], X.1[ind.1,])

              Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
              Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

              model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              nfolds=nrow(train)
              model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

              c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
              ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

              rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

              return(c(c,ct))
            }
            stopCluster(cluster.cores)
            closeAllConnections()
            rm(cluster.cores)
            rm(cores)
            gc()

            df_C_index <- do.call(rbind, list_C_index)
            df_C_index <- data.frame(df_C_index)
            colnames(df_C_index) <- c("C_index_test", "C_index_train")
            df_C_index$Iteration <- 1:iter
            df_C_index$ID <- i

            df_C.index.alldata <- rbind(df_C.index.alldata, df_C_index)

            cat("\nAnalysis completed for: ", i,"\n")
            rm(X,Y,  df_C_index, list_C_index, opts, pb, p)

          }
        }
      }
    }
  }
  return(list(C_index_results = df_C.index.alldata))
}
#'
#' Cox forecasting glmnet, combination data
#'
#' C-index testing of Cox partial likelihood glmnet models for survival forecasting from ex vivo drug screens and genetic/clinical data.
#'
#' @param X_data Input drug screen data
#' @param y_data Survival data
#' @param Z_data Input other clinical and genetic covariates.
#' @param alpha Penalty type. 0 for L2 (Ridge), 1 for L1 (Lasso), or any number between 0 and 1 for elastic net penalty.
#' @param lambda Vector of regularization penalties.
#' @param free_cores Number of free cores.
#' @param test.n Vector of test set balance. First integer is the number of survivals, and second integer is the number of deaths.
#' @param iter Number of test iterations.
#' @param log_AUC log2 transform of AUC. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Patient.Z Patient-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param Drug.Z Drug-wise standardization. 1 for yes, 2 for no, or c(1,2) for both.
#' @param RCPC Removal of confounding principal components.
#' @importFrom magrittr %>%
#'
#' @return
#' A list containing C-index test results for X and Z individually and combined.
#'
#'
#' @author
#' Aram N. Andersen \email{aram.n.andersen@@gmail.com}
#'
#' @export
Cox_forecasting_glmnet_combination <- function(X_data,
                                               y_data,
                                               Z_data,
                                               alpha=0,
                                               lambda=c(exp(seq(-4,6, 0.1))),
                                               free_cores = 2,
                                               test.n= c(6,4),
                                               iter=200,
                                               log_AUC=c(1:2),
                                               Patient.Z=c(1:2),
                                               Drug.Z =c(1:2),
                                               RCPC=c(0,1,2,3,4)){

  df_C.index.alldata <- c()

  # Clinical/genetic data
  for(a in alpha){
    i=paste0("Clinical/genetic", "/Penalty_", a)
    X <- Z_data
    Y <- y_data

    cores <- parallel::detectCores()-free_cores
    cluster.cores<-makeCluster(cores)
    registerDoSNOW(cluster.cores)
    pb <- txtProgressBar(max=iter, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
      set.seed(b)
      setTxtProgressBar(pb,b)
      X.0 <- X[Y[,2]==0,]
      X.1 <- X[Y[,2]==1,]
      Y.0 <- Y[Y[,2]==0,]
      Y.1 <- Y[Y[,2]==1,]

      ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
      ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

      train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
      test <- rbind(X.0[ind.0,], X.1[ind.1,])

      Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
      Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

      model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
      nfolds=nrow(train)
      model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

      c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
      ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

      rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

      return(c(c,ct))
    }
    stopCluster(cluster.cores)
    closeAllConnections()
    rm(cluster.cores)
    rm(cores)
    gc()

    df_C_index <- do.call(rbind, list_C_index)
    df_C_index <- data.frame(df_C_index)
    colnames(df_C_index) <- c("C_index_test", "C_index_train")
    df_C_index$Iteration <- 1:iter
    df_C_index$ID <- i
    df_C_index$Data <- "Clinical/genetic"

    df_C.index.alldata <- rbind(df_C.index.alldata, df_C_index)

    cat("\nAnalysis completed for: ", i,"\n")
    rm(X,Y,  df_C_index, list_C_index, opts, pb, p)


  }

  # Combined data
  for(a in alpha){
    for(Transform in c("log2(AUC)", "AUC")[log_AUC]){
      for(pt.st in c(TRUE, FALSE)[Patient.Z] ){
        for(cpc in RCPC){
          for(drug.st in c(TRUE, FALSE)[Drug.Z]){
            i=paste0("Clinical/genetic/", Transform,c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
            set.seed(1)
            X <- X_data
            l=min(data.frame(replace(X, X == 0, 1)))
            if(Transform=="log2(AUC)"){
              if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
            }else{
              X <- data.matrix(X)
            }
            if(pt.st){
              X <- (X - rowMeans(X))/matrixStats::rowSds(X)
            }
            if(cpc>0){
              x_sd <- matrixStats::colSds(X)
              x_mean <- colMeans(X)
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
              svdz <- svd(t(X))
              X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
              X <- X*x_sd + x_mean
              X <- t(X)
              rm(svdz, x_mean, x_sd)
            }

            if(drug.st){
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
            }
            X <- cbind(X, Z_data)
            Y <- y_data

            cores <- parallel::detectCores()-free_cores
            cluster.cores<-makeCluster(cores)
            registerDoSNOW(cluster.cores)
            pb <- txtProgressBar(max=iter, style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
              set.seed(b)
              setTxtProgressBar(pb,b)
              X.0 <- X[Y[,2]==0,]
              X.1 <- X[Y[,2]==1,]
              Y.0 <- Y[Y[,2]==0,]
              Y.1 <- Y[Y[,2]==1,]

              ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
              ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

              train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
              test <- rbind(X.0[ind.0,], X.1[ind.1,])

              Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
              Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

              model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              nfolds=nrow(train)
              model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

              c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
              ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

              rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

              return(c(c,ct))
            }
            stopCluster(cluster.cores)
            closeAllConnections()
            rm(cluster.cores)
            rm(cores)
            gc()

            df_C_index <- do.call(rbind, list_C_index)
            df_C_index <- data.frame(df_C_index)
            colnames(df_C_index) <- c("C_index_test", "C_index_train")
            df_C_index$Iteration <- 1:iter
            df_C_index$ID <- i
            df_C_index$Data <- "Combined"

            df_C.index.alldata <- rbind(df_C.index.alldata, df_C_index)

            cat("\nAnalysis completed for: ", i,"\n")
            rm(X,Y,  df_C_index, list_C_index, opts, pb, p)

          }
        }
      }
    }
  }

  # Drug screen data
  for(a in alpha){
    for(Transform in c("log2(AUC)", "AUC")[log_AUC]){
      for(pt.st in c(TRUE, FALSE)[Patient.Z] ){
        for(cpc in RCPC){
          for(drug.st in c(TRUE, FALSE)[Drug.Z]){
            i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
            set.seed(1)
            X <- X_data
            l=min(data.frame(replace(X, X == 0, 1)))
            if(Transform=="log2(AUC)"){
              if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
            }else{
              X <- data.matrix(X)
            }
            if(pt.st){
              X <- (X - rowMeans(X))/matrixStats::rowSds(X)
            }
            if(cpc>0){
              x_sd <- matrixStats::colSds(X)
              x_mean <- colMeans(X)
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
              svdz <- svd(t(X))
              X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
              X <- X*x_sd + x_mean
              X <- t(X)
              rm(svdz, x_mean, x_sd)
            }

            if(drug.st){
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
            }

            Y <- y_data

            cores <- parallel::detectCores()-free_cores
            cluster.cores<-makeCluster(cores)
            registerDoSNOW(cluster.cores)
            pb <- txtProgressBar(max=iter, style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
              set.seed(b)
              setTxtProgressBar(pb,b)
              X.0 <- X[Y[,2]==0,]
              X.1 <- X[Y[,2]==1,]
              Y.0 <- Y[Y[,2]==0,]
              Y.1 <- Y[Y[,2]==1,]

              ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
              ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

              train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
              test <- rbind(X.0[ind.0,], X.1[ind.1,])

              Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
              Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

              model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              nfolds=nrow(train)
              model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

              c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
              ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

              rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

              return(c(c,ct))
            }
            stopCluster(cluster.cores)
            closeAllConnections()
            rm(cluster.cores)
            rm(cores)
            gc()

            df_C_index <- do.call(rbind, list_C_index)
            df_C_index <- data.frame(df_C_index)
            colnames(df_C_index) <- c("C_index_test", "C_index_train")
            df_C_index$Iteration <- 1:iter
            df_C_index$ID <- i
            df_C_index$Data <- "Drug screen"

            df_C.index.alldata <- rbind(df_C.index.alldata, df_C_index)

            cat("\nAnalysis completed for: ", i,"\n")
            rm(X,Y,  df_C_index, list_C_index, opts, pb, p)

          }
        }
      }
    }
  }


  return(list(C_index_results = df_C.index.alldata))
}
#'
#' Cox glmnet model bootstrapping for survival association hypothesis testing
#'
#' Survival association statistics.
#'
#' @param X_data Input drug screen data
#' @param y_data Survival data
#' @param alpha Penalty type. 0 for L2 (Ridge), 1 for L1 (Lasso), or any number between 0 and 1 for elastic net penalty.
#' @param lambda Vector of regularization penalties if pre.CV = FALSE.
#' @param pre.CV Logical. pre.CV=TRUE if optimal lambda is provided.
#' @param lambda_opt Optimal lambda if pre.CV=TRUE.
#' @param free_cores Number of free cores.
#' @param test.n Vector of test set balance. First integer is the number of survivals, and second integer is the number of deaths.
#' @param iter Number of test iterations.
#' @param log_AUC log2 transform of AUC. 1 for yes, 2 for no.
#' @param Patient.Z Patient-wise standardization. 1 for yes, 2 for no.
#' @param Drug.Z Drug-wise standardization. 1 for yes, 2 for no.
#' @param RCPC Removal of confounding principal components.
#' @importFrom magrittr %>%
#'
#' @return
#' A list containing coefficients and associated statistics.
#'
#' @author
#' Aram N. Andersen \email{aram.n.andersen@@gmail.com}
#'
#' @export
Cox_bootstrapping <- function(X_data,
                              y_data,
                              alpha=0,
                              lambda=c(exp(seq(-4,6, 0.1))),
                              pre.CV=TRUE,
                              lambda_opt = 0,
                              free_cores = 2,
                              iter=200,
                              log_AUC=c(1:2),
                              Patient.Z=c(1:2),
                              Drug.Z =c(1:2),
                              RCPC=c(0,1,2,3,4)){
  # Data preparation
  a = alpha
  Transform = c("log2(AUC)", "AUC")[log_AUC[1]]
  pt.st = c(TRUE, FALSE)[Patient.Z[1]]
  cpc = RCPC[1]
  drug.st = c(TRUE, FALSE)[Drug.Z[1]]
  i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
  set.seed(1)
  X <- X_data
  l=min(data.frame(replace(X, X == 0, 1)))
  if(Transform=="log2(AUC)"){
    if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
  }else{
    X <- data.matrix(X)
  }
  if(pt.st){
    X <- (X - rowMeans(X))/matrixStats::rowSds(X)
  }
  if(cpc>0){
    x_sd <- matrixStats::colSds(X)
    x_mean <- colMeans(X)
    X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
    if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
    svdz <- svd(t(X))
    X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
    X <- X*x_sd + x_mean
    X <- t(X)
    rm(svdz, x_mean, x_sd)
  }
  if(drug.st){
    X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
    if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
  }
  Y <- y_data
  colnames(X) <- colnames(X_data)
  X_data <- X

  # Prediction function
  Boot_cox <- function(X,
                       Y,
                       test.n=test.n,
                       pre.CV=pre.CV,
                       lambda_opt = lambda_opt,
                       a=a,lambda=lambda, iter=iter, i=i){

    cores <- parallel::detectCores()-free_cores
    cluster.cores<-makeCluster(cores)
    registerDoSNOW(cluster.cores)
    pb <- txtProgressBar(max=iter, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_coefs <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
      set.seed(b)
      setTxtProgressBar(pb,b)

      ids <- sample(nrow(X), nrow(X), replace=TRUE)
      X_boot <- X[ids,]
      Y_boot <- Y[ids,]

      if(pre.CV == FALSE){
        model.loop <- glmnet(X_boot, Y_boot, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
        nfolds=nrow(X_boot)
        model.loop.cv <- cv.glmnet(X_boot, Y_boot, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
        coefs  <- data.frame(as.matrix(coef(model.loop, s = model.loop.cv$lambda.min)))
      }else{
        model.loop <- glmnet(X_boot, Y_boot, family = "cox", alpha = a, standardize = FALSE,lambda=lambda_opt,  type.measure = "deviance")
        coefs  <- data.frame(as.matrix(coef(model.loop, s =lambda_opt)))
      }
      coefs$boot <- b
      coefs$drug <- rownames(coefs)

      return(coefs)
    }
    stopCluster(cluster.cores)
    closeAllConnections()
    rm(cluster.cores)
    rm(cores)
    gc()

    return(list_coefs)
  }


  list_coefficients <- Boot_cox(X,
                                Y,
                                test.n=test.n,
                                pre.CV=pre.CV,
                                lambda_opt = lambda_opt,
                                a=a,lambda=lambda, iter=iter, i=i)

  df_coefficients <- bind_rows(list_coefficients)
  colnames(df_coefficients)[1] <- "Beta"

  df_coefficients.statistics <- df_coefficients %>% group_by(drug) %>%
    summarise(Mean = mean(Beta),
              StErr = sd(Beta))
  df_coefficients.statistics <- data.frame(df_coefficients.statistics)
  df_coefficients.statistics$Lower_95CI = df_coefficients.statistics$Mean - df_coefficients.statistics$StErr*1.96
  df_coefficients.statistics$Upper_95CI = df_coefficients.statistics$Mean + df_coefficients.statistics$StErr*1.96
  df_coefficients.statistics$t.stat <- df_coefficients.statistics$Mean/df_coefficients.statistics$StErr
  df_coefficients.statistics$Rank <- rank(df_coefficients.statistics$t.stat)

  return(list(Betas=df_coefficients,
              Betas_statistics=df_coefficients.statistics))

}
#'
#' Cox forecasting glmnet, drug withdrawal testing
#'
#' C-index variable importance testing of Cox partial likelihood glmnet model for survival forecasting from ex vivo drug screens.
#'
#' @param X_data Input drug screen data
#' @param y_data Survival data
#' @param Reduce Variable removal for dimension reduction. "None" for only testing the variables. "Naive" reduce the model by removing one and one variable starting from highest C-index gain until optimal model is found. "Iterative" will re-test variables as the model is reduced (UNDER DEVELOPMENT).
#' @param alpha Penalty type. 0 for L2 (Ridge), 1 for L1 (Lasso), or any number between 0 and 1 for elastic net penalty.
#' @param lambda Vector of regularization penalties.
#' @param free_cores Number of free cores.
#' @param test.n Vector of test set balance. First integer is the number of survivals, and second integer is the number of deaths.
#' @param iter Number of test iterations.
#' @param log_AUC log2 transform of AUC. 1 for yes, 2 for no.
#' @param Patient.Z Patient-wise standardization. 1 for yes, 2 for no.
#' @param Drug.Z Drug-wise standardization. 1 for yes, 2 for no.
#' @param RCPC Removal of confounding principal components.
#' @importFrom magrittr %>%
#'
#' @return
#' A list containing C-index test results from the variable withdrawal routine.
#'   Contains the test results from the parent data, and the C-index differential upon drug withdrawal.
#'   Optimization results from removing non-predictive variables.
#'
#' @author
#' Aram N. Andersen \email{aram.n.andersen@@gmail.com}
#'
#' @export
Cox_forecasting_drug_withdrawal <- function(X_data,
                                            y_data,
                                            Reduce=c("None","Naive", "Iterative"),
                                            alpha=0,
                                            lambda=c(exp(seq(-4,6, 0.1))),
                                            free_cores = 2,
                                            test.n= c(6,4),
                                            iter=200,
                                            log_AUC=c(1:2),
                                            Patient.Z=c(1:2),
                                            Drug.Z =c(1:2),
                                            RCPC=c(0,1,2,3,4),
                                            path=path){
  
  Drug_WD_test <- function(X_data,
                           y_data,
                           Reduce=c("None","Naive"),
                           alpha=0,
                           lambda=c(exp(seq(-4,6, 0.1))),
                           free_cores = 2,
                           test.n= c(6,4),
                           iter=200,
                           log_AUC=c(1:2),
                           Patient.Z=c(1:2),
                           Drug.Z =c(1:2),
                           RCPC=c(0,1,2,3,4)){
    # Data preparation
    a = alpha
    Transform = c("log2(AUC)", "AUC")[log_AUC[1]]
    pt.st = c(TRUE, FALSE)[Patient.Z[1]]
    cpc = RCPC[1]
    drug.st = c(TRUE, FALSE)[Drug.Z[1]]
    i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
    set.seed(1)
    X <- X_data
    l=min(data.frame(replace(X, X == 0, 1)))
    if(Transform=="log2(AUC)"){
      if(length(which(X <= 0))){
        X <- data.matrix(-log2(X + l))
      }else{
        X <- data.matrix(-log2(X))
      }
    }else{
      X <- data.matrix(X)
    }
    if(pt.st){
      X <- (X - rowMeans(X))/matrixStats::rowSds(X)
    }
    if(cpc>0){
      x_sd <- matrixStats::colSds(X)
      x_mean <- colMeans(X)
      X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
      if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
      svdz <- svd(t(X))
      X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
      X <- X*x_sd + x_mean
      X <- t(X)
      rm(svdz, x_mean, x_sd)
    }
    if(drug.st){
      X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
      if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
    }
    Y <- y_data
    colnames(X) <- colnames(X_data)
    X_data <- X
    
    # Prediction function
    Pred_cox <- function(X,
                         Y,
                         test.n=test.n,
                         a=a,lambda=lambda, iter=iter, i=i){
      
      cores <- parallel::detectCores()-free_cores
      cluster.cores<-makeCluster(cores)
      registerDoSNOW(cluster.cores)
      pb <- txtProgressBar(max=iter, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
        set.seed(b)
        setTxtProgressBar(pb,b)
        X.0 <- X[Y[,2]==0,]
        X.1 <- X[Y[,2]==1,]
        Y.0 <- Y[Y[,2]==0,]
        Y.1 <- Y[Y[,2]==1,]
        
        ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
        ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)
        
        train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
        test <- rbind(X.0[ind.0,], X.1[ind.1,])
        
        Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
        Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])
        
        model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
        nfolds=nrow(train)
        model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
        
        c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
        ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))
        
        rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
        
        return(c(c,ct))
      }
      stopCluster(cluster.cores)
      closeAllConnections()
      rm(cluster.cores)
      rm(cores)
      gc()
      
      df_C_index <- do.call(rbind, list_C_index)
      df_C_index <- data.frame(df_C_index)
      colnames(df_C_index) <- c("C_index_test", "C_index_train")
      df_C_index$ID <- i
      
      rm(X,Y, list_C_index, opts, pb)
      
      
      return(df_C_index)
    }
    
    # Run on full data
    df_parent <- Pred_cox(X, Y, test.n, a,lambda, iter, i)
    
    # Prepare drug withdrawal data
    mat.C.index.test <- data.frame(array(data=NA, dim=c(nrow(df_parent), ncol(X_data))))
    colnames(mat.C.index.test) <- colnames(X_data)
    mat.C.index.train <- mat.C.index.test
    
    cat("\n", "Parent prediction completed...")
    
    for(dr in colnames(X_data)){
      X_data1 <- X_data[,colnames(X_data) != dr]
      df_wd <- Pred_cox(X_data1, Y, test.n, a,lambda, iter, i)
      mat.C.index.test[,dr] <- df_wd$C_index_test
      mat.C.index.train[,dr] <- df_wd$C_index_train
      cat("\n",dr, "withdrawal prediction completed...")
      rm(df_wd, X_data1)
    }
    
    C_loss_test <- colMeans((mat.C.index.test - df_parent$C_index_test))
    C_loss_train <- colMeans((mat.C.index.train - df_parent$C_index_train))
    
    df_loss <- data.frame(C_loss_test=C_loss_test,
                          C_loss_train=C_loss_train,
                          Var=names(C_loss_train))
    
    
    
    Optimization <- list()
    if("Naive" %in% Reduce ){
      cat("\nInitiating naive model reduction:\n")
      
      dr=df_loss$Var[which.max(df_loss$C_loss_test)]
      X_data1 <- X_data[,colnames(X_data) != dr]
      df_parent1 <- df_parent
      df_parent1$C_index_test <- mat.C.index.test[,dr]
      df_parent1$C_index_train <- mat.C.index.train[,dr]
      loss <- mean(df_parent1$C_index_test - df_parent$C_index_test)
      df_loss1 <- df_loss[df_loss$Var != dr,]
      list_naive_reduction <- list()
      list_naive_reduction[[1]] <- df_parent
      list_naive_reduction[[2]] <- df_parent1
      list_naive_reduction[[1]]$Data <- "WD0_full"
      list_naive_reduction[[2]]$Data <- paste0("WD1_", dr)
      j=3
      while(loss>=0 & ncol(X_data1)>2){
        dr=df_loss1$Var[which.max(df_loss1$C_loss_test)]
        X_data1 <- X_data1[,colnames(X_data1) != dr]
        
        df_wd <- Pred_cox(X_data1, Y, test.n, a,lambda, iter, i)
        loss <- mean(df_wd$C_index_test - df_parent1$C_index_test)
        
        list_naive_reduction[[j]] <- df_wd
        list_naive_reduction[[j]]$Data <- paste0("WD",j-1,"_", dr)
        
        df_loss1 <- df_loss1[df_loss1$Var != dr,]
        df_parent1 <- df_wd
        j=j+1
      }
      
      df_naive_reduction <- bind_rows(list_naive_reduction)
      df_naive_reduction$Data <- factor(df_naive_reduction$Data, levels=unique(df_naive_reduction$Data))
      
      Optimization[["Naive_reduction"]] <- df_naive_reduction
      
      cat("\nNaive model reduction completed...\n")
    }
    
    return(list(df_parent = df_parent,
                WD_initial = list(df_loss=df_loss,
                                  mat.C.index.test=mat.C.index.test,
                                  mat.C.index.train=mat.C.index.train),
                WD_optimization = Optimization))
  }
  
  if(Reduce != "Iterative"){
    list_res_initial <- Drug_WD_test(X_data, 
                                     y_data, 
                                     Reduce=Reduce,
                                     alpha=alpha, 
                                     lambda=lambda,
                                     free_cores = free_cores,
                                     test.n= test.n, 
                                     iter=iter,
                                     log_AUC=log_AUC,
                                     Patient.Z=Patient.Z,
                                     Drug.Z =Drug.Z,
                                     RCPC=RCPC)
    if(Reduce == "Naive"){
      df_loss <- list_res_initial$WD_initial$df_loss
      df_naive_reduction <- list_res_initial$WD_optimization$Naive_reduction
      df_naive_reduction_sum <- df_naive_reduction %>% group_by(ID, Data) %>% summarise(Mean=mean(C_index_test), Median=median(C_index_test))
      remove_stop <- which.max(df_naive_reduction_sum$Mean)
      remove_drugs <-gsub(".*_","", df_naive_reduction_sum$Data[2:remove_stop])
      remaining_drugs <- colnames(X)
      remaining_drugs <- remaining_drugs[!(remaining_drugs %in% remove_drugs)]
      X_data_red <- X_data[,which(colnames(X_data) %in% remaining_drugs)]
      list_res_initial$WD_optimization$Naive_reduction_summary <- df_naive_reduction_sum
      list_res_initial$WD_optimization$remove_drugs <- remove_drugs
      list_res_initial$WD_optimization$remaining_drugs <- remaining_drugs
      list_res_initial$WD_optimization$X_reduced <- X_data_red
      return(list_res_initial)
    }else{
      return(list_res_initial)
    }
    save(list_res_initial, file=path)
  }else{
    list_res_initial <- Drug_WD_test(X_data, 
                                     y_data, 
                                     Reduce="Naive",
                                     alpha=alpha, 
                                     lambda=lambda,
                                     free_cores = free_cores,
                                     test.n= test.n, 
                                     iter=iter,
                                     log_AUC=log_AUC,
                                     Patient.Z=Patient.Z,
                                     Drug.Z =Drug.Z,
                                     RCPC=RCPC)
    df_loss <- list_res_initial$WD_initial$df_loss
    df_naive_reduction <- list_res_initial$WD_optimization$Naive_reduction
    df_naive_reduction_sum <- df_naive_reduction %>% group_by(ID, Data) %>% summarise(Mean=mean(C_index_test), Median=median(C_index_test))
    remove_stop <- which.max(df_naive_reduction_sum$Mean)
    remove_drugs <-gsub(".*_","", df_naive_reduction_sum$Data[2:remove_stop])
    remaining_drugs <- colnames(X)
    remaining_drugs <- remaining_drugs[!(remaining_drugs %in% remove_drugs)]
    X_data_red <- X_data[,which(colnames(X_data) %in% remaining_drugs)]
    list_res_initial$WD_optimization$Naive_reduction_summary <- df_naive_reduction_sum
    list_res_initial$WD_optimization$X_reduced <- X_data_red
    
    list_res_iterations <- list()
    list_res_iterations[["WD_reduction_0"]] <- list_res_initial
    c = 1
    while(ncol(X_data_red)<ncol(X_data) & ncol(X_data_red)>2){
      cat("\nIterative model reduction:",c,"\n")
      X_data <- X_data_red
      list_res_iter <- Drug_WD_test(X_data, 
                                    y_data, 
                                    Reduce="Naive",
                                    alpha=alpha, 
                                    lambda=lambda,
                                    free_cores = free_cores,
                                    test.n= test.n, 
                                    iter=iter,
                                    log_AUC=log_AUC,
                                    Patient.Z=Patient.Z,
                                    Drug.Z =Drug.Z,
                                    RCPC=RCPC)
      
      X_data <- X_data_red
      df_loss <- list_res_iter$WD_initial$df_loss
      df_naive_reduction <- list_res_iter$WD_optimization$Naive_reduction
      df_naive_reduction_sum <- df_naive_reduction %>% group_by(ID, Data) %>% summarise(Mean=mean(C_index_test), Median=median(C_index_test))
      remove_stop <- which.max(df_naive_reduction_sum$Mean)
      remove_drugs <-gsub(".*_","", df_naive_reduction_sum$Data[2:remove_stop])
      remaining_drugs <- colnames(X)
      remaining_drugs <- remaining_drugs[!(remaining_drugs %in% remove_drugs)]
      X_data_red <- X_data[,which(colnames(X_data) %in% remaining_drugs)]
      list_res_iter$WD_optimization$Naive_reduction_summary <- df_naive_reduction_sum
      list_res_iter$WD_optimization$X_reduced <- X_data_red
      list_res_iterations[[paste0("WD_reduction_", c)]] <- list_res_iter
      save(list_res_iterations, file=path)
      c = c + 1
    }
    return(list_res_iterations)
  }
}
#'
