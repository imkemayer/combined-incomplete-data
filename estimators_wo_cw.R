################################################################################
## Estimators for RCT data
################################################################################

risk_ratio <- function(DF, outcome_name="Y", treatment_name="A"){
  
  treated_idx <- which(DF[,treatment_name] == 1)
  control_idx <- which(DF[,treatment_name] == 0)
  
  IE <- sum(DF[treated_idx, outcome_name])
  IN <- sum(1 - DF[treated_idx, outcome_name])  
  CE <- sum(DF[control_idx, outcome_name]) 
  CN <- sum(1 - DF[control_idx, outcome_name])
  
  RR <- (IE*(CE+CN)) / (CE*(IE+IN))
  
  CI <- 1.96*sqrt(IN/(IE*(IE+IN)) +  CN/(CE*(CE+CN)))
  
  return(c(risk_placebo = 100*CE /(CE+CN), risk_TXA = 100*IE/(IE+IN), RR = RR, lower_ci = RR - CI , upper_ci = RR + CI))
}



# to have the equivalence in ATE
difference_in_means <- function(DF, outcome_name="Y", treatment_name="A") {
  
  treated_idx <- which(DF[,treatment_name] == 1)
  control_idx <- which(DF[,treatment_name] == 0)
  
  # Filter treatment / control observations, pulls outcome variable as a vector
  y1 <- DF[treated_idx, outcome_name] # Outcome in treatment grp
  y0 <- DF[control_idx, outcome_name] # Outcome in control group
  
  n1 <- length(treated_idx)     # Number of obs in treatment
  n0 <- length(control_idx) # Number of obs in control
  
  # Difference in means is ATE
  tauhat <- mean(y1) - mean(y0)
  
  # 95% Confidence intervals
  se_hat <- sqrt( var(y0)/(n0-1) + var(y1)/(n1-1) )
  lower_ci <- tauhat - 1.96 * se_hat
  upper_ci <- tauhat + 1.96 * se_hat
  
  return(c(ATE = tauhat, lower_ci = lower_ci, upper_ci = upper_ci))
}

difference_in_condmeans_ols <- function(DF, outcome_name="Y", treatment_name="A") {
  complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
  DF<-DF[which(complete_rows),]
  n <- dim(DF)[1]
  treated_idx <- which(DF[,treatment_name] == 1)
  control_idx <- which(DF[,treatment_name] == 0)
  
  # Filter treatment / control observations, pulls outcome variable as a vector
  y1 <- DF[treated_idx, outcome_name] # Outcome in treatment grp
  y0 <- DF[control_idx, outcome_name] # Outcome in control group
  
  n1 <- sum(DF[,treatment_name])     # Number of obs in treatment
  n0 <- sum(1 - DF[,treatment_name]) # Number of obs in control
  
  dataset_spread <- DF
  df0 <- dataset_spread[control_idx,] %>%
            dplyr::select(-c(treatment_name)) %>%
            mutate(outcome=eval(parse(text=outcome_name))) %>%
            dplyr::select(-c(outcome_name))
  df1 <- dataset_spread[treated_idx,] %>%
            dplyr::select(-c(treatment_name)) %>%
            mutate(outcome=eval(parse(text=outcome_name))) %>%
            dplyr::select(-c(outcome_name))
  
  mu0 <- lm(outcome ~ ., data = df0)
  mu1 <- lm(outcome ~ ., data = df1)
  
  # Difference in predicted means is ATE
  tauhat <- mean(predict(mu1, newdata = dplyr::select(dataset_spread, -c(treatment_name, outcome_name))) - predict(mu0, newdata = dplyr::select(dataset_spread, -c(treatment_name, outcome_name))), na.rm = TRUE)
  
  # 95% Confidence intervals
  beta0 <- mu0$coefficients[-1]
  beta0[which(is.na(beta0))] <- 0
  beta1 <- mu1$coefficients[-1]
  beta1[which(is.na(beta1))] <- 0
  se_hat <- sqrt( n*(var(y0)/(n0) + var(y1)/(n1)) - t(beta0+beta1)%*%var(dplyr::select(dataset_spread, -c(treatment_name,outcome_name)))%*%(beta0+beta1))/sqrt(n)
  lower_ci <- tauhat - 1.96 * se_hat
  upper_ci <- tauhat + 1.96 * se_hat
  
  return(c(ATE = tauhat, lower_ci = lower_ci, upper_ci = upper_ci))
}

compute_mean_diff_RCT <- function(DF, outcome_name="Y", treatment_name="A"){
  RCT_ATE <- mean(DF[DF[,treatment_name] == 1 & DF$V == 1, outcome_name]) - mean(DF[DF[,treatment_name] == 0  & DF$V == 1, outcome_name])  
  return(RCT_ATE)
}

################################################################################
## Nuisance functions and 
## Estimators for generalization from RCT to observational data
################################################################################

sampling_propensities <- function(DF, method="glm", seed=100){
  complete <- !any(is.na(DF))
  
  # Logistic regression
  if (method=="glm"){
    # Regular logistic regression in case of complete data
    if (complete){
      p.fit <- glm(V~., data = DF,
                   family = "binomial")
      pi_s_hat <- predict(p.fit, type = "response")
    # EM for logistic regression in case of incomplete data
    } else {
      p.fit <- miss.glm(V~., data = DF, print_iter=F, seed = seed)
      pi_s_hat <- predict(p.fit, newdata = DF[, !names(DF) %in% c("V")])
    } 
  # Random forest regression
  } else if (method=="grf"){
    na.action <- options()$na.action
    options(na.action = 'na.pass')
    if (is.data.frame(DF)){
      X.m = model.matrix(~.-1, data = DF[, !names(DF) %in% c("V")])
    } else { # X can also be a matrix
      X.m = model.matrix(~.-1, data=data.frame(DF[, !names(DF) %in% c("V")]))
    }
    options(na.action = na.action)
    
    forest.W = regression_forest(X.m, DF$V, tune.parameters = "all")
    pi_s_hat = predict(forest.W)$predictions
  } else {
    print("Incorrect 'method' parameter, must be either 'glm' or 'grf'.")
    return(NULL)
  }
  return(pi_s_hat)
}

outcome_regressions <- function(DF, method="glm", outcome_name="Y", treatment_name="A", seed=100){
  complete <- !any(is.na(DF[, !names(DF) %in% c(outcome_name, treatment_name)]))
  temp <- DF
  temp$outcome <- temp[,outcome_name]
  temp$treatment <- temp[,treatment_name]
  temp <- dplyr::select(temp, -c(outcome_name,treatment_name))
  
  binary_y <- FALSE
  if (length(unique(temp$outcome))==2) binary_y <- TRUE
  if (method=="glm") {
    # Regular logistic regression in case of complete data
    if (complete) {
      if (binary_y){
        mu_1 <- glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 1, !names(temp) %in% c("V", "treatment")], family="binomial")
        mu_0 <- glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 0, !names(temp) %in% c("V", "treatment")], family="binomial")
        
        mu_1_hat <- predict(mu_1, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")], type="response")
        mu_0_hat <- predict(mu_0, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")], type="response")
        rct_mu_11_hat <- predict(mu_1, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")], type="response")
        rct_mu_10_hat <- predict(mu_0, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")], type="response")
      } else {
        mu_1 <- glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 1, !names(temp) %in% c("V", "treatment")])
        mu_0 <- glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 0, !names(temp) %in% c("V", "treatment")])
        
        mu_1_hat <- predict(mu_1, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")])
        mu_0_hat <- predict(mu_0, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")])
        rct_mu_11_hat <- predict(mu_1, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")])
        rct_mu_10_hat <- predict(mu_0, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")])
      }
    } else {
      if (binary_y){
        mu_1 <- miss.glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 1, !names(temp) %in% c("V", "treatment")], print_iter=F, seed=100)
        mu_0 <- miss.glm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 0, !names(temp) %in% c("V", "treatment")], print_iter=F, seed=100)
      } else {
        mu_1 <- miss.lm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 1, !names(temp) %in% c("V", "treatment")])
        mu_0 <- miss.lm(outcome ~., data = temp[temp$V == 1 & temp$treatment == 0, !names(temp) %in% c("V", "treatment")])
      } 
      mu_1_hat <- predict(mu_1, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")])
      mu_0_hat <- predict(mu_0, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "treatment", "outcome")])
      rct_mu_11_hat <- predict(mu_1, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")])
      rct_mu_10_hat <- predict(mu_0, newdata = temp[which(temp$V == 1), !names(temp) %in% c("V", "treatment", "outcome")])
      
    }
  # Random forest regression
  } else if (method=="grf") {
    X <- temp[, !names(temp) %in% c("V", "treatment", "outcome")]
    # the two lines below are necessary if there are missing values in X
    na.action <- options()$na.action
    if (!complete){
      options(na.action = 'na.pass')
    } 
    if (is.data.frame(X)){
      X.m = model.matrix(~.-1, data = X)
    } else { # X can also be a matrix
      X.m = model.matrix(~.-1, data=data.frame(X))
    }
    options(na.action = na.action)
    
    idx_0 <- which(temp$V == 1 & temp$treatment == 0)
    mu_0 = regression_forest(X.m[idx_0, ], temp[idx_0, "outcome"], tune.parameters = "all")
    
    idx_1 <- which(temp$V == 1 & temp$treatment == 1)
    mu_1 = regression_forest(X.m[idx_1, ], temp[idx_1, "outcome"], tune.parameters = "all")
    
    mu_1_hat <- predict(mu_1, newdata = X.m[which(temp$V==0),])$predictions
    mu_0_hat <- predict(mu_0, newdata = X.m[which(temp$V==0),])$predictions
    rct_mu_11_hat <- predict(mu_1, newdata = X.m[which(temp$V == 1),])$predictions
    rct_mu_10_hat <- predict(mu_0, newdata = X.m[which(temp$V == 1),])$predictions
  } else {
    print("Incorrect 'method' parameter, must be either 'glm' or 'grf'.")
    return()
  }
  return(list("mu_1_hat"=mu_1_hat,"mu_0_hat"=mu_0_hat,
              "rct_mu_11_hat"=rct_mu_11_hat, "rct_mu_10_hat"=rct_mu_10_hat))
}

compute_ipsw <- function(DF, outcome_name="Y", treatment_name="A",
                         method="glm",
                         complete_cases=FALSE,
                         covariates="all"){
  ## simply delete this paragraph when removal of the covariates parameter
  if (covariates == "all"){
    temp <- DF
  } else if (covariates == "X1"){
    temp <- DF[, c("X1", "V", treatment_name, outcome_name)]
  } else if (covariates == "-X1"){
    temp <- DF[, !(names(DF) %in%  c("X1"))]
  } else {
    print("Covariates parameters must be all, X1, -X1.")
    break
  }
  DF <- temp
  
  # Either only keep complete observations
  if (complete_cases){
    complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
    temp<-DF[which(complete_rows),]
  # Or keep all observations and handle NA later
  } else{
    temp <- DF
  }
  temp$outcome <- temp[,outcome_name]
  temp$treatment <- temp[,treatment_name]
  temp <- dplyr::select(temp, -c(outcome_name, treatment_name))
  
  m <- nrow(temp[temp$V ==0, ])
  
  # Estimation of P(V = 1 | X)
  # p <-- P(V = 1 | X) 
  pi_s_hat <- sampling_propensities(temp[,!names(temp) %in% c("outcome", "treatment")], method=method)
  
  # Store odds
  temp$odds <- ((1 - pi_s_hat)/pi_s_hat)
  
  # Keep only RCT for the rest of the computations
  temp <- temp[which(temp$V == 1),]
  
  tau_ipsw_hat <- (2/m)*with(temp, sum(odds*treatment*outcome - odds*(1-treatment)*outcome))
  tau_ipsw_norm_hat <- with(temp, sum(odds*treatment*outcome/sum(odds*treatment) - odds*(1-treatment)*outcome/sum(odds*(1-treatment))))
  return(c(tau_ipsw_hat, tau_ipsw_norm_hat))
}

compute_gformula <- function(DF, outcome_name="Y", treatment_name="A",
                             method="glm",
                             complete_cases=FALSE) {
  # Either only keep complete observations
  if (complete_cases){
    complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
    temp<-DF[which(complete_rows),]
  # Or keep all observations and handle NA later
  } else{
    temp <- DF
  }
  
  mus_hat <- outcome_regressions(temp, method=method, outcome_name=outcome_name, treatment_name=treatment_name)
  
  tau_gformula_hat <- mean(mus_hat[['mu_1_hat']]) - mean(mus_hat[['mu_0_hat']])
  return(tau_gformula_hat)
}  

compute_aipsw <- function(DF, outcome_name="Y", treatment_name="A",
                          method="glm", normalized=FALSE,
                          complete_cases=FALSE) {
  # Either only keep complete observations
  if (complete_cases){
    complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
    temp<-DF[which(complete_rows),]
    # temp <- na.omit(DF)
  # Or keep all observations and handle NA later
  } else{
    temp <- DF
  }
    
  m <- nrow(temp[temp$V ==0, ])
    
  # G-formula part
  mus_hat <- outcome_regressions(temp, method=method, outcome_name=outcome_name, treatment_name=treatment_name)
  tau_gformula <- mean(mus_hat[['mu_1_hat']]) - mean(mus_hat[['mu_0_hat']])
  
  # IPSW part
  pi_s_hat <- sampling_propensities(temp[,!names(temp) %in% c(outcome_name, treatment_name)], method=method)
  
  # Store odds
  temp$odds <- ((1 - pi_s_hat)/pi_s_hat)
    
  rct <- temp[which(temp$V == 1),]
  rct$mu_11 <-  mus_hat[['rct_mu_11_hat']]     
  rct$mu_10 <-  mus_hat[['rct_mu_10_hat']]
  rct$outcome <- rct[,outcome_name]
  rct$treatment <- rct[,treatment_name]
  rct <- dplyr::select(rct, -c(outcome_name, treatment_name))
  
  if (normalized == FALSE){
    tau_ipsw <- (2/m)*with(rct, sum(odds*treatment*(outcome - mu_11) - odds*(1-treatment)*(outcome - mu_10)))  
  } else {
    tau_ipsw <- with(rct, sum(odds*treatment*(outcome - mu_11)/sum(odds*treatment) - odds*(1-treatment)*(outcome - mu_10)/sum(odds*(1-treatment))))
  }
    
  tau_aipsw_hat <- tau_ipsw + tau_gformula
  return(tau_aipsw_hat)
}

# STRATIFICATION
compute_stratification <- function(DF, nb_strat = 10, bin = "quantile",
                                   outcome_name="Y", treatment_name="A",
                                   method="glm",
                                   complete_cases=FALSE,
                                   pi_strat=NULL){
  # Only keep complete observations
  if (complete_cases){
    complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
    temp<-DF[which(complete_rows),]
    # keep all observations and handle NA later
  } else{
    temp <- DF
  }
  
  temp$outcome <- temp[,outcome_name]
  temp$treatment <- temp[,treatment_name]
  temp <- dplyr::select(temp, -c(outcome_name,treatment_name))
  
  # logit : sampling score
  # temp$V <- as.numeric(temp$V)
  if (is.null(pi_strat)){
    pi_strat <- sampling_propensities(temp[,!names(temp) %in% c("outcome", "treatment")], method=method)
    
    # following lines are equivalent to predict
    #pi_s_coefficient <- pi_s_reg$coefficients
    #X <- as.matrix(temp[, !names(temp) %in% c("Y", "A", "S")])
    #pi_strat <- as.numeric(expit(cbind(1, X) %*% pi_s_coefficient))
  }  
  temp$pi_strat <- (1-pi_strat)/pi_strat
  
  # decompose in strata 
  if (bin == "quantile"){
    temp$strata <- as.numeric(cut_number(temp$pi_strat, nb_strat))
  } 
  else if (bin == "regular"){
    temp$strata <- as.numeric(cut(temp$pi_strat, nb_strat))
  }
  else {
    print("Wrong `bin` argument, must be either 'quantile' or 'regular'")
    break
  }
  
  
  rct <- temp[temp$V ==1,]
  m <- nrow(temp[temp$V == 0,])
  tau_strat <- 0
  
  
  for (s in unique(rct$strata)){
    # compute strata ate
    strata_ate <- mean(rct[rct$strata == s & rct$treatment == 1, "outcome"]) - mean(rct[rct$strata == s & rct$treatment == 0, "outcome"])
    weigth <- nrow(temp[temp$V ==0 & temp$strata == s, ]) / m
    tau_strat <- tau_strat + weigth*strata_ate
    
  }
  
  return(tau_strat)
}


compute_all <- function(DF, outcome_name="Y", treatment_name="A",
                        method="glm", nb_strat=10, bin="quantile",
                        complete_cases=FALSE,
                        verbose=FALSE) {
  # Only keep complete observations
  if (complete_cases){
    complete_rows <- apply(dplyr::select(DF, -c(outcome_name, treatment_name)), 1, function(x) all(!is.na(x)))
    temp<-DF[which(complete_rows),]
    # keep all observations and handle NA later
  } else{
    temp <- DF
  }
  
  m <- nrow(temp[temp$V ==0, ])
  
  # G-formula part
  if (verbose) cat("G-formula start")
  mus_hat <- outcome_regressions(temp, method=method, outcome_name=outcome_name, treatment_name=treatment_name)
  tau_gformula_hat <- mean(mus_hat[['mu_1_hat']]) - mean(mus_hat[['mu_0_hat']])
  
  # IPSW part
  if (verbose) {cat(", G-formula end.",sep = "\n"); cat("IPSW start")}
  pi_s_hat <- sampling_propensities(temp[,!names(temp) %in% c(outcome_name, treatment_name)], method=method)
  
  tau_strat_hat <- compute_stratification(temp, nb_strat = nb_strat, bin = bin,
                                     outcome_name=outcome_name, treatment_name=treatment_name,
                                     method=method,
                                     pi_strat=pi_s_hat)
  # Store odds
  temp$odds <- ((1 - pi_s_hat)/pi_s_hat)
  
  # Keep only RCT for the rest of the calculus
  rct <- temp[temp$V == 1,]
  
  # IPSW estimator
  rct$outcome <- rct[,outcome_name]
  rct$treatment <- rct[,treatment_name]
  rct <- dplyr::select(rct, -c(outcome_name, treatment_name))
  tau_ipsw_hat <- (2/m)*with(rct, sum(odds*treatment*outcome - odds*(1-treatment)*outcome))
  tau_ipsw_norm_hat <- with(rct, sum(odds*treatment*outcome/sum(odds*treatment) - odds*(1-treatment)*outcome/sum(odds*(1-treatment))))
  
  if (verbose) {cat(", IPSW end.",sep = "\n"); cat("AIPSW start")}
  rct$mu_11 <-  mus_hat[['rct_mu_11_hat']]     
  rct$mu_10 <-  mus_hat[['rct_mu_10_hat']]
  
  
  tau_ipsw <- with(rct, sum(odds*treatment*(outcome - mu_11)/sum(odds*treatment) - odds*(1-treatment)*(outcome - mu_10)/sum(odds*(1-treatment))))
  
  #tau_ipsw <- (2/m)*with(rct, sum(odds*treatment*(outcome - mu_11) - odds*(1-treatment)*(outcome - mu_10)))  
  
  tau_aipsw_hat <- tau_ipsw + tau_gformula_hat
  if (verbose) cat(", AIPSW end.",sep = "\n")
  
  return(list("ipsw_hat" = tau_ipsw_hat, 
              "ipsw.norm_hat" = tau_ipsw_norm_hat,
              "gformula_hat" = tau_gformula_hat,
              "aipsw_hat" = tau_aipsw_hat,
              "strat_hat" = tau_strat_hat))
}

  