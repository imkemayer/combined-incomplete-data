library(reshape2)
library(MASS)
library(dplyr)
library(misaem)  # for glm (linear and logistic) with missing data
library(grf)     # for generalized random forests, option for incomplete data
library(devtools)
library(pracma)
library(genRCT)  # for Calibration Weighting (Dong's implementation)
library(nleqslv) # for searchZeros function required in genRCT package
library(micemd)  # for mice on multilevel data
library(purrr)   # for the map function used in data.frame handling
source_url("https://rmisstastic.netlify.app/how-to/generate/amputation.R")
source("estimators.R")

library(norm)


# Dong's transformation to mispecify
transform_X_star <-function(DF){
  DF_transformed <- DF
  DF_transformed$X1 <- exp(DF$X1 / 3)
  DF_transformed$X2 <- DF$X2 / (1+exp(DF$X1)) + 10
  DF_transformed$X3 <- DF$X1 * DF$X3 / 25 + 0.6
  DF_transformed$X4 <- DF$X1 + DF$X4 + 20
  matrix_transformed <- as.matrix(DF_transformed[, c("X1", "X2", "X3", "X4")])
  matrix_transformed <- scale(matrix_transformed, center = c(1,1,1,1))
  return(matrix_transformed)
}

# simulation with continuous outcome
simulate_continuous <- function(n = 1000, m = 49000, p = 4, mu = rep(1, p), Sigma = diag(p), 
                                bs = c(-0.5, -0.3, -0.5, -0.4), bs0 = -2.5, 
                                beta = c(27.4, 13.7, 13.7, 13.7), b0 = - 100, sigma = 1, snr=NULL,
                                misRCT = "correct", misoutcome = "correct", Nested = FALSE,
                                na_rct = NULL, # missing values in RCT
                                na_rwe = NULL, # missing values in RWE
                                na_df = NULL, # missing values in both, needs to be specified for CIS=T assumption
                                link="linear") {
  
  # Target population generation 
  covariates <- mvrnorm(n = 50*n, mu, Sigma, tol = 1e-06, empirical = FALSE) # 50*n is roughly the initial population size necessary to have the n
  DF <- as.data.frame(covariates)
  names(DF) <- paste("X", 1:p, sep = "")
  covariates_names <- names(DF)
  DF_orig <- DF
  if (link=="nonlinear"){
    tmp <- DF
    for (j in 1:p){
      if (mod(j,3)==0){
        tmp[,j] <- tmp[,j] * tmp[,(j-2)]
      }
      if (mod(j,p)==1){
        tmp[,j] <- sin(tmp[,j])*sign(tmp[,j])+1.5
      }
      if (mod(j,2)==0){
        tmp[,j] <- sign(tmp[,j])*mean(tmp[,j], na.rm=T)
      }
    }
    DF <- tmp
  }
  tau = beta[1]*mean(DF$X1)
  
  DF_orig_full <- DF_orig
  DF_full <- DF
  
  df_nas_idx <- NULL
  if (!is.null(na_df)){
    self.mask <- NULL
    mechanism <- na_df[['mechanism']]
    if (na_df[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_df[['cis']]){
      if (length(c(na_df[['prop_miss']]))>1){
        df_nas_idx <- matrix(FALSE, nrow=nrow(DF_orig), ncol=length(covariates_names))
        for (prop.miss in unique(na_df[['prop_miss']])){
          df_nas <- produce_NA(data=DF_orig[,covariates_names],
                               mechanism = mechanism,
                               self.mask = self.mask,
                               perc.missing = prop.miss,
                               idx.incomplete = as.numeric(na_df[['prop_miss']]==prop.miss))
          df_nas_idx <- apply(pmax(df_nas_idx,df_nas[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        df_nas <- produce_NA(data=DF_orig[,covariates_names],
                             mechanism = mechanism,
                             self.mask = self.mask,
                             perc.missing = na_df[['prop_miss']],
                             idx.incomplete = na_df[['idx_incomplete']])
        df_nas_idx <- df_nas[['idx_newNA']]
      }
      DF[,covariates_names][df_nas_idx] <- 0
    }
  }
  covariates <- as.matrix(DF)
  
  if (!is.null(na_df)){
    if (na_df[['cio']]){
      tau = beta[1]*mean(DF$X1)
    }
  }
  
  # RCT probability to sample according to model
  if (misRCT == "correct"){
    etas <- as.vector(covariates %*% bs + bs0)
    
  } else if (misRCT == "exponential") {
    
    # RCT misspecification with exp on all covariates
    etas <- as.vector(exp(covariates) %*% bs + bs0 + 3) # 3 was found manually to keep same proportion m and n
    
  } else if (misRCT == "Partial_X2only"){
    
    # partial misspecification with only X2 affected
    DF_mis <- DF
    DF_mis$X2 <- exp(DF_mis$X2)
    mis_covariatesX2 <- as.matrix(DF_mis[, c("X1", "X2", "X3", "X4")])
    etas <- as.vector(mis_covariatesX2 %*% bs + bs0 + 0.1)
    
  } else if (misRCT == "Partial_X1only"){
    
    # partial misspecification with only X1 affected
    DF_mis <- DF
    DF_mis$X1 <- exp(DF_mis$X1)
    mis_covariatesX1 <- as.matrix(DF_mis[, c("X1", "X2", "X3", "X4")])
    etas <- as.vector(mis_covariatesX1 %*% bs + bs0 + 0.1)
    
  } else if(misRCT == "dong"){
    
    miscovariatesDong <- transform_X_star(DF)
    etas <- as.vector (miscovariatesDong %*% bs + bs0)
    
  } else if (misRCT == "strongbias"){
    bs = c(-1.5, -0.3, -0.5, -0.4)
    etas <- as.vector (covariates %*% bs + bs0)
  }  else {
    
    print("Error in RCT specification arguments.")
    break
    
  }
  
  ps = 1 / (1 + exp(-etas))
  DF$ps <- ps
  
 
  
  # from probability to RCT indicator
  RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
  DF$V <- RCT_indicator 
  
  # random treatment assignment within the RCT
  DF$A <- ifelse(DF$V == 1, rbinom(nrow(DF), 1, 0.5), NA)
  
  # keep only interesting variables
  DF <- DF[, c(covariates_names, "A", "V")]
  
  if (!Nested) {
    
    # drop other data
    DF_rct <- DF[DF$V == 1,] 
    
    # untransformed covariates (if link=="linear", then this is the same as DF_rct)
    DF_orig_rct <- data.frame(cbind(DF_orig[DF$V==1, ], DF$A[DF$V==1 ],DF$V[DF$V==1 ]))
    
    # untransformed covariates without missing values
    DF_orig_full_rct <- data.frame(cbind(DF_orig_full[DF$V==1, ], DF$A[DF$V==1 ],DF$V[DF$V==1 ]))
    
    names(DF_orig_rct) <- names(DF_orig_full_rct) <- names(DF_rct)
    
    # generate new observational data
    covariates_rwe <- mvrnorm(n = m, mu, Sigma, tol = 1e-06, empirical = FALSE) 
    DF_rwe <- as.data.frame(covariates_rwe)
    names(DF_rwe) <- paste("X", 1:p, sep = "")
    DF_orig_rwe <- DF_rwe
    DF_orig_full_rwe <- DF_rwe
    if (link=="nonlinear"){
      tmp <- DF_rwe
      for (j in 1:p){
        if (mod(j,3)==0){
          tmp[,j] <- tmp[,j] * tmp[,(j-2)]
        }
        if (mod(j,p)==1){
          tmp[,j] <- sin(tmp[,j])*sign(tmp[,j])+1.5
        }
        if (mod(j,2)==0){
          tmp[,j] <- sign(tmp[,j])*mean(tmp[,j], na.rm=T)
        }
      }
      DF_rwe <- tmp
    }
    DF_rwe$A <- rep(NA, m)
    DF_rwe$V <- rep(0, m)
    
    
    DF_full_rct <- DF_full[DF$V == 1, ]
    DF_full_rct$A <- DF_rct$A
    DF_full_rct$V <- DF_rct$V
    
    DF_full_rwe <- DF_rwe[, names(DF_rct)]
    DF_full_rwe$A <- rep(NA,m)
    DF_full_rwe$V <- rep(0,m)
    
    DF_orig_full_rct <- DF_orig_full[DF$V == 1, ]
    DF_orig_full_rct$A <- DF_orig_rct$A
    DF_orig_full_rct$V <- DF_orig_rct$V
    
    DF_orig_rwe <- data.frame(cbind(DF_orig_rwe, rep(NA, m),rep(0, m)))
    names(DF_orig_rwe) <- names(DF_rwe)
    DF_orig_full_rwe <- DF_orig_rwe
    
    if (is.null(df_nas_idx)){
      df_nas_idx_rct <- matrix(FALSE, nrow=n, ncol=length(covariates_names))
    } else {
      df_nas_idx_rct <- df_nas_idx[which(DF$V==1),]
    }
    df_nas_idx_rwe <- matrix(FALSE, nrow=m, ncol=length(covariates_names))
    
    
  } else {
    
    #here we need to drop values such that the final data set contains m observational values and n RCT values.
    DF_rct <- DF[DF$V == 1,]
    DF_rwe <- DF[DF$V == 0,]
    
    DF_orig_rct <- data.frame(cbind(DF_orig[DF$V==1, ], DF$A[DF$V==1 ],DF$V[DF$V==1 ]))
    DF_orig_rwe <- data.frame(cbind(DF_orig[DF$V==0, ], DF$A[DF$V==0 ],DF$V[DF$V==0 ]))
    names(DF_orig_rct) <- names(DF_rct)
    names(DF_orig_rwe) <- names(DF_rwe)
    
    DF_full_rct <- DF_full[DF$V == 1,]
    DF_full_rwe <- DF_full[DF$V == 0,]
    DF_orig_full_rct <- DF_orig_full[DF$V == 1,]
    DF_orig_full_rwe <- DF_orig_full[DF$V == 0,]
    
    df_nas_idx_rct <- df_nas_idx[which(DF$V==1),]
    df_nas_idx_rwe <- df_nas_idx[which(DF$V==0),]
  }
  
  rct_nas <- rwe_nas <- NULL
  if (!is.null(na_rct)){
    self.mask <- NULL
    mechanism <- na_rct[['mechanism']]
    if (na_rct[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_rct[['cio']]){
      if (length(c(na_rct[['prop_miss']]))>1){
        rct_nas_idx <- matrix(FALSE, nrow=nrow(DF_orig_rct), ncol=length(covariates_names))
        for (prop.miss in unique(na_rct[['prop_miss']])){
          rct_nas <- produce_NA(data=DF_orig_rct[,covariates_names],
                               mechanism = mechanism,
                               self.mask = self.mask,
                               perc.missing = prop.miss,
                               idx.incomplete = as.numeric(na_rct[['prop_miss']]==prop.miss))
          rct_nas_idx <- apply(pmax(rct_nas_idx, rct_nas[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        rct_nas <- produce_NA(data=DF_orig_rct[,covariates_names],
                              mechanism = mechanism,
                              self.mask = self.mask,
                              perc.missing = na_rct[['prop_miss']],
                              idx.incomplete = na_rct[['idx_incomplete']])
        rct_nas_idx <- rct_nas[['idx_newNA']]
      }
      DF_rct[,covariates_names][rct_nas_idx] <- 0
    }
  }
  if (!is.null(na_rwe)){
    self.mask <- NULL
    mechanism <- na_rwe[['mechanism']]
    if (na_rwe[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_rwe[['cio']]){
      if (length(c(na_rwe[['prop_miss']]))>1){
        rwe_nas_idx <- matrix(FALSE, nrow=nrow(DF_orig_rwe), ncol=length(covariates_names))
        for (prop.miss in unique(na_rwe[['prop_miss']])){
          rwe_nas <- produce_NA(data=DF_orig_rwe[,covariates_names],
                                mechanism = mechanism,
                                self.mask = self.mask,
                                perc.missing = prop.miss,
                                idx.incomplete = as.numeric(na_rwe[['prop_miss']]==prop.miss))
          rwe_nas_idx <- apply(pmax(rwe_nas_idx, rwe_nas[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        rwe_nas <- produce_NA(data=DF_orig_rwe[,covariates_names],
                             mechanism = mechanism,
                             self.mask = self.mask,
                             perc.missing = na_rwe[['prop_miss']],
                             idx.incomplete = na_rwe[['idx_incomplete']])
        rwe_nas_idx <- rwe_nas[['idx_newNA']]
      }
      DF_rwe[,covariates_names][rwe_nas_idx] <- 0
    }
  }
  
  if (!is.null(na_rct)){
    if (na_rct[['cio']]){
      tau = beta[1]*mean(DF$X1)
    }
  }
  
  # stack RCT and RWE
  DF <- rbind(DF_rct, DF_rwe)
  DF_full <- rbind(DF_full_rct, DF_full_rwe)
  DF_orig <- rbind(DF_orig_rct, DF_orig_rwe)
  DF_orig_full <- rbind(DF_orig_full_rct, DF_orig_full_rwe)
  
  
  
  # reset row number
  rownames(DF) <- 1:nrow(DF)
  rownames(DF_full) <- 1:nrow(DF_full)
  rownames(DF_orig) <- 1:nrow(DF_orig)
  rownames(DF_orig_full) <- 1:nrow(DF_orig_full)
  
  if (!is.null(na_df)){
    if (!na_df[['cio']] & is.null(rwe_nas) & is.null(rct_nas)) {
      DF[,covariates_names] <- DF_full[,covariates_names]
    }
  }
  
  # compute Y  
  ## if the argument `snr` is provided, we compute a corresponding sigma for the outcome model
  if (!is.null(snr)){
    sigma <- sqrt(t(beta)%*%var(DF_full[,1:4])%*%beta/snr)
  }
  if (misoutcome == "correct"){
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF$X1 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error
  } 
  else if (misoutcome == "+a"){
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + DF$X1 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error + beta[1]*(DF$A == 1)
  }
  else if (misoutcome == "wrong")
  {
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF$X1*DF$X2 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error
  } 
  else if (misoutcome == "dong") {
    miscovariatesDong <- transform_X_star(DF)
    DF_Xstar <- as.data.frame(miscovariatesDong)
    names(DF_Xstar) <- paste("X", 1:p, sep = "")
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF_Xstar$X1 + beta[2]*DF_Xstar$X2 + beta[3]*DF_Xstar$X3 + beta[4]*DF_Xstar$X4 + error
  }
  else {
    print("Parameters misoutcome is badly specified")
    break
  }
  
  DF_full$Y <- DF$Y
  DF_orig$Y <- DF$Y
  DF_orig_full$Y <- DF$Y
  
  if (!is.null(na_df)){
    self.mask <- NULL
    mechanism <- na_df[['mechanism']]
    if (na_df[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_df[['cis']]){
      DF[, covariates_names][rbind(df_nas_idx_rct, df_nas_idx_rwe)] <- NA
      DF_orig[, covariates_names][rbind(df_nas_idx_rct, df_nas_idx_rwe)] <- NA
    } else {
      if (length(c(na_df[['prop_miss']]))>1){
        idx_tmp <- matrix(FALSE, nrow=nrow(DF_orig), ncol=length(covariates_names))
        for (prop.miss in unique(na_df[['prop_miss']])){
          tmp <- produce_NA(data=DF_orig[,covariates_names],
                            mechanism = mechanism,
                            self.mask = self.mask,
                            perc.missing = prop.miss,
                            idx.incomplete = as.numeric(na_df[['prop_miss']]==prop.miss))
          idx_tmp <- apply(pmax(idx_tmp, tmp[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        idx_tmp <- produce_NA(data=DF_orig[,covariates_names],
                              mechanism = mechanism,
                              self.mask = self.mask,
                              perc.missing = na_df[['prop_miss']],
                              idx.incomplete = na_df[['idx_incomplete']])[['idx_newNA']]
      }
      DF[,covariates_names][idx_tmp] <- NA
      DF_orig[, covariates_names][idx_tmp] <- NA
    }
  }
  
  if (!is.null(na_rct)){
    self.mask <- NULL
    mechanism <- na_rct[['mechanism']]
    if (na_rct[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_rct[['cio']]){
      DF[DF$V == 1, covariates_names][rct_nas_idx] <- NA
      DF_orig[DF$V == 1, covariates_names][rct_nas_idx] <- NA
    } else {
      if (length(c(na_rct[['prop_miss']]))>1){
        idx_tmp <- matrix(FALSE, nrow=sum(DF$V == 1), ncol=length(covariates_names))
        for (prop.miss in unique(na_rct[['prop_miss']])){
          tmp <- produce_NA(data=DF_orig[DF$V == 1,covariates_names],
                            mechanism = mechanism,
                            self.mask = self.mask,
                            perc.missing = prop.miss,
                            idx.incomplete = as.numeric(na_rct[['prop_miss']]==prop.miss))
          idx_tmp <- apply(pmax(idx_tmp, tmp[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        idx_tmp <- produce_NA(data=DF_orig[DF$V == 1,covariates_names],
                              mechanism = mechanism,
                              self.mask = self.mask,
                              perc.missing = na_rct[['prop_miss']],
                              idx.incomplete = na_rct[['idx_incomplete']])[['idx_newNA']]
      }
      DF[DF$V == 1,covariates_names][idx_tmp] <- NA
      DF_orig[DF$V == 1,covariates_names][idx_tmp] <- NA
    }
  }
  if (!is.null(na_rwe)){
    self.mask <- NULL
    mechanism <- na_rwe[['mechanism']]
    if (na_rwe[['mechanism']] == "MNAR_selfmask") { 
      self.mask <- "upper"
      mechanism <- "MNAR" 
    }
    if (na_rwe[['cio']]){
      DF[DF$V == 0, covariates_names][rwe_nas_idx] <- NA
      DF_orig[DF$V == 0, covariates_names][rwe_nas_idx] <- NA
    } else {
      if (length(c(na_rwe[['prop_miss']]))>1){
        idx_tmp <- matrix(FALSE, nrow=sum(DF$V == 0), ncol=length(covariates_names))
        for (prop.miss in unique(na_rwe[['prop_miss']])){
          tmp <- produce_NA(data=DF_orig[DF$V == 0,covariates_names],
                            mechanism = mechanism,
                            self.mask = self.mask,
                            perc.missing = prop.miss,
                            idx.incomplete = as.numeric(na_rwe[['prop_miss']]==prop.miss))
          idx_tmp <- apply(pmax(idx_tmp, tmp[['idx_newNA']]), c(1,2), as.logical)
        }
      } else {
        idx_tmp <- produce_NA(data=DF_orig[DF$V == 0,covariates_names],
                              mechanism = mechanism,
                              self.mask = self.mask,
                              perc.missing = na_rwe[['prop_miss']],
                              idx.incomplete = na_rwe[['idx_incomplete']])[['idx_newNA']]
      }
      DF[DF$V == 0,covariates_names][idx_tmp] <- NA
      DF_orig[DF$V == 0,covariates_names][idx_tmp] <- NA
    }
  }
  if (link == "nonlinear"){
    DF <- DF_orig
  }
  return(list("DF"=DF,"DF_full"=DF_orig_full, "tau"=tau))
}

compute_estimators_mi <- function(DF, method="glm", nb_mi=5, 
                                  outcome_name="Y", treatment_name="A",
                                  tau=NULL, strategy="within",
                                  nb_strat=1, bin="quantile",
                                  complete_cases=F, verbose=F, verbose_intern=F) {
  # strategy is one of "within", "joint-fixed-effect-wY", "joint-fixed-effect-woY",
  # "multi-level-woY", "joint-naive-wY", "joint-naive-woY"
  
  nb_mi_max <- max(c(nb_mi))
  
  if (strategy=="within"){
    DF_rct <- DF[which(DF$V==1),]
    DF_rwe <- DF[which(DF$V==0),]
    
    m0 <- mice(DF_rct, maxit=0, print=F)
    predMatrix <- m0$pred
    predMatrix[,c("V")] <- 0
    if (verbose) cat("\n\tMI for RCT, ")
    DF_rct_mice <- mice.par(DF_rct,predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
    if (verbose) cat(" end.", end="\n")
  
    m0 <- mice(DF_rwe, maxit=0, print=F)
    predMatrix <- m0$pred
    predMatrix[,c("V", treatment_name, outcome_name)] <- 0
    if (verbose) cat("\n\tMI for RWE, ")
    DF_rwe_mice <- mice.par(DF_rwe, predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
    if (verbose) cat(" end.", end="\n")
    
    nb_mi <- c(0, nb_mi)
    results <- c()
    for (meth in c(method)){
      res_mi <- c()
      if (verbose) cat("\tUse RCT imputation ")
      for (nb_mi_idx in 1:(length(nb_mi)-1)){
        for (i_rct_imp in (nb_mi[nb_mi_idx]+1):nb_mi[nb_mi_idx+1]){
          if (verbose) cat(i_rct_imp, end=", ")
          for (i_rwe_imp in (nb_mi[nb_mi_idx]+1):nb_mi[nb_mi_idx+1]){
            DF_imp <- rbind(complete(DF_rct_mice,i_rct_imp), complete(DF_rwe_mice,i_rwe_imp))
            rct_ate <- mean(DF_imp[which(DF_imp$A == 1 & DF_imp$V == 1), "Y"]) - mean(DF_imp[which(DF_imp$A == 0  & DF_imp$V == 1), "Y"])
            
            res_tmp <- compute_all(DF_imp, outcome_name=outcome_name, treatment_name=treatment_name,
                                   method=meth, nb_strat=1)
            
            res_mi <- rbind(res_mi, 
                            data.frame("i_rct_imp"=i_rct_imp,
                                       "i_rwe_imp"=i_rwe_imp,
                                       "RCT" = rct_ate,
                                       "IPSW" = res_tmp$ipsw_hat,
                                       "IPSW.norm" = res_tmp$ipsw.norm_hat,
                                       "G.formula" = res_tmp$gformula_hat,
                                       "AIPSW" = res_tmp$aipsw_hat,
                                       "CW" = res_tmp$cw_hat))
          }
        }
        
        results <- rbind(results, data.frame("RCT" = mean(res_mi[,"RCT"]),
                                             "IPSW" = mean(res_mi[,"IPSW"]),
                                             "IPSW norm" = mean(res_mi[,"IPSW.norm"]),
                                             "G-formula" = mean(res_mi[,"G.formula"]),
                                             "AIPSW" = mean(res_mi[,"AIPSW"]),
                                             "CW" = mean(res_mi[,"CW"]),
                                             "tau" = tau,
                                             "method" = meth,
                                             "nb_mi"=nb_mi[nb_mi_idx+1]))
      }
    }
    if (verbose) cat(" end.", end="\n")
  } else {
    if (strategy == "joint-fixed-effect-wY"){
      
      m0 <- mice(DF, maxit=0, print=F)
      predMatrix <- m0$pred
      predMatrix[,c(treatment_name)] <- 0
      if (verbose) cat("\n\tMI on joint dataset with source indicator (fixed effects model), ")
      DF_mice <- mice.par(DF,predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
      if (verbose) cat(" end.", end="\n")
      
    }
    if (strategy == "joint-fixed-effect-woY"){
      
      m0 <- mice(DF, maxit=0, print=F)
      predMatrix <- m0$pred
      predMatrix[,c(outcome_name, treatment_name)] <- 0
      if (verbose) cat("\n\tMI on joint dataset with source indicator (fixed effects model), ")
      DF_mice <- mice.par(DF,predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
      if (verbose) cat(" end.", end="\n")
    }
    if (strategy == "multi-level-woY"){
      ind.clust<-which(colnames(DF)=="V") #index for the cluster variable
      #initialisation of the argument predictorMatrix
      m0 <- mice(DF, maxit=0, print=F)
      predMatrix <- m0$pred
      predMatrix[ind.clust,ind.clust]<-0
      predMatrix[-ind.clust,ind.clust]<- -2
      predMatrix[predMatrix==1]<-2
      predMatrix[,c(treatment_name, outcome_name)] <- 0
      #initialisation of the argument method
      micemd_method<-find.defaultMethod(DF,ind.clust)
      micemd_method[treatment_name] <- ""
      micemd_method[outcome_name] <- ""
      if (verbose) cat("\n\tMI on joint dataset with 2-level models, ")
      #multiple imputation by chained equations (parallel calculation)
      DF_mice <- mice.par(DF,predictorMatrix = predMatrix,
                         method=micemd_method,m = nb_mi_max,printFlag=verbose_intern)
      if (verbose) cat(" end.", end="\n")
    }
    if (strategy == "multi-level-wY"){
      print("TBA")
    }
    if (strategy == "joint-naive-wY"){
      
      m0 <- mice(DF, maxit=0, print=F)
      predMatrix <- m0$pred
      predMatrix[,c(treatment_name,"V")] <- 0
      if (verbose) cat("\n\tMI on joint dataset without source indicator, ")
      DF_mice <- mice.par(DF,predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
      if (verbose) cat(" end.", end="\n")
      
    }
    if (strategy == "joint-naive-woY"){
      
      m0 <- mice(DF, maxit=0, print=F)
      predMatrix <- m0$pred
      predMatrix[,c("V",outcome_name, treatment_name)] <- 0
      if (verbose) cat("\n\tMI on joint dataset without source indicator, ")
      DF_mice <- mice.par(DF,predictorMatrix = predMatrix, m=nb_mi_max, printFlag=verbose_intern)
      if (verbose) cat(" end.", end="\n")
      
    }
    
    results <- c()
    if (verbose) cat("\n\tUse joint imputations ")
    for (nb in nb_mi){
      for (meth in methods){
        rct_ate_tmp <- DF_mice %>%
          mice::complete(1:nb, mild = TRUE) %>%
          map(filter, V==1) %>%  
          map(function(x) mean(x[which(x[,treatment_name] == 1), outcome_name]) - mean(x[which(x[,treatment_name] == 0), outcome_name])) 
        rct_ate <- Reduce("+", rct_ate_tmp) / length(rct_ate_tmp)
        
        res_tmp <- DF_mice %>%
          mice::complete(1:nb, mild = TRUE) %>%
          map(compute_all, outcome_name=outcome_name, treatment_name=treatment_name,
              method=meth, nb_strat=1) 
        res_mi <- data.frame(t(apply(sapply(data.frame(do.call(rbind,res_tmp)), unlist), 2, mean)))
        
        results <- rbind(results, data.frame("RCT" = rct_ate,
                                             "IPSW" = res_mi$ipsw_hat,
                                             "IPSW norm" = res_mi$ipsw.norm_hat,
                                             "G-formula" = res_mi$gformula_hat,
                                             "AIPSW" = res_mi$aipsw_hat,
                                             "CW" = res_mi$cw_hat,
                                             "tau" = tau,
                                             "method" = meth,
                                             "nb_mi"=nb))
      }
    }
    if (verbose) cat(" end.", end="\n")
  }
  return(results)
}

# Function that launches rep times the simulation and returns a data.frame with results
compute_estimators_and_store <- function(rep, misoutcome = "correct", misRCT = "correct", N = 50000, m = 49000, n = 1000, bs0=-2.5,
                                         p = 4, mu = rep(1, p), Sigma = diag(p), sigma=1, snr=NULL,
                                         na_df=NULL, na_rct=NULL, na_rwe=NULL, link="linear", method="glm", nb_strat=10, bin="quantile",
                                         complete_cases=FALSE, full_data=FALSE,
                                         do_mi=FALSE, nb_mi=5, strategy="within",
                                         verbose=FALSE, verbose_intern=FALSE,
                                         seed=NULL){
  
  # rct_ate <- c()
  # ipsw <- c()
  # ipsw_norm <- c()
  # strat_10 <- c()
  # gformula <- c()
  # cw <- c()
  # aipsw <- c()
  # tau <- c()
  
  results <- c()
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (verbose) cat("Iteration ")
  
  for (i in 1:rep){
    if (verbose) cat(i, end=", ")
    temp <- simulate_continuous(misoutcome = misoutcome, misRCT = misRCT, m = m, n = n, 
                                p=p, mu=mu, Sigma=Sigma, sigma=sigma, snr=snr,
                                na_df=na_df, na_rct = na_rct, na_rwe = na_rwe, link=link, bs0=bs0)
    if (!full_data){
      DF <- temp$DF
    } else {
      DF <- temp$DF_full
    }
    if (verbose) cat(paste0(" (mean(X1)=", round(mean(DF$X1, na.rm=T),2), ")",
                            " (var(Y)=", round(var(DF$Y, na.rm=T),2), ")"), end=", ")
    
    # results <- c()
    
    if (do_mi){
      res_mi <- compute_estimators_mi(DF, nb_mi=nb_mi, method=method, strategy=strategy,
                                      outcome_name="Y", treatment_name="A",
                                      tau=temp$tau,
                                      nb_strat=nb_strat, bin=bin,
                                      complete_cases=complete_cases, verbose=verbose_intern)
      if (complete_cases){
        complete_rows <- apply(dplyr::select(DF, -c("Y", "A")), 1, function(x) all(!is.na(x)))
        DF.tmp<-DF[which(complete_rows),]
        n_effective <- sum(DF.tmp$V==1)
        m_effective <- sum(DF.tmp$V==0)
      } else {
        n_effective <- sum(DF$V==1)
        m_effective <- sum(DF$V==0)
      }
      res_mi$n_effective <- n_effective
      res_mi$m_effective <- m_effective

      results <- rbind(results, res_mi)
    } else {
      for (meth in c(method)){
      
        
      # naive estimator
        rct_ate <- mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"])
        
        res_tmp <- compute_all(DF, outcome_name="Y", treatment_name="A",
                                method=meth, nb_strat=nb_strat, bin=bin,
                                complete_cases=complete_cases, verbose=verbose_intern)
          
        #ispw
        #ipsw  <- c(ipsw, res_tmp$ipsw_hat)
        #ipsw_norm <- c(ipsw_norm, res_tmp$ipsw.norm_hat)
        ipsw <- res_tmp$ipsw_hat
        ipsw_norm <- res_tmp$ipsw.norm_hat
        
        #strat
        #strat_10 <- c(strat_10, res_tmp$strat_hat)
        
        #gformula
        #gformula <- c(gformula, res_tmp$gformula_hat)
        gformula <- res_tmp$gformula_hat
        
        #aipsw
        #aipsw <- c(aipsw, res_tmp$aipsw_hat)
        aipsw <- res_tmp$aipsw_hat
        
        #cw
        cw <- res_tmp$cw_hat
        
        # tau
        #tau <- c(tau, temp$tau)
        tau <- temp$tau
        
        if (complete_cases){
          complete_rows <- apply(dplyr::select(DF, -c("Y", "A")), 1, function(x) all(!is.na(x)))
          DF.tmp<-DF[which(complete_rows),]
          n_effective <- sum(DF.tmp$V==1)
          m_effective <- sum(DF.tmp$V==0)
        } else {
          n_effective <- sum(DF$V==1)
          m_effective <- sum(DF$V==0)
        }
        results <- rbind(results, data.frame("RCT" = rct_ate,
                              "IPSW" = ipsw,
                              "IPSW norm" = ipsw_norm,
                              #"Stratification n=10" = strat_10,
                              "G-formula" = gformula,
                              "AIPSW" = aipsw,
                              "CW" = cw,
                              "tau" = tau,
                              "method" = meth,
                              "n_effective"=n_effective,
                              "m_effective"=m_effective))
      }
    }
    if (verbose) cat(paste0(i,"-end."), end="\n")
  }
  return(results)
}


