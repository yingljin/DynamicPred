library(cubature)
# This script writes functions for 
# calculating conditional expectation
# of fPCA scores 
# of an out-of-sample observation

#### Numeric integration: denominator of E(xi|Y) ####
# xi: scores. The value to integrate over
# tao: estimated covariance matrix(eigenvalues) of score from fPCA
# f0: mean function from fPCA (up to mth bin)
# phi: basis function of fPCA (up to mth bin)
# df_out: newly observed data (up to mth bin)

  
joint_pdf <- function(xi, tao, f0, phi, df_out){
  
  # numbers
  ns <- table(df_out$bin) # number of observations
  hs <- df_out %>% group_by(bin) %>% summarize_at("Y", sum) %>% 
    select(Y) # number of success
  nf <- ns-hs # number of failure
  
  # probability of conditional distribution of Y
  eta = f0+phi%*%xi
  pr = exp(eta)/(1+exp(eta))
  # the part relevant to P(Y|xi)
  p1 <- prod(dbinom(x=unlist(hs), size = ns, prob=pr))
  p2 <- exp(-t(xi)%*%solve(tao)%*%xi/2)
  
  return(p1*p2)
  
}

#test 
# joint_pdf(xi, tao, f0, phi, df_out)
# 
# adaptIntegrate(joint_pdf, lower = rep(-500, K), upper = rep(500, K), 
#           tao=tao, f0=f0, phi=phi, df_out=df_out)
  
# takes a really long time to run!

#### Laplace approximation with Bayes package ####

library(LaplacesDemon)
library(mvtnorm)

## model
Model <- function(parm, Data){
  xi <- parm[Data$pos.xi]
  
  # log-prior
  xi.prior <- dmvnorm(xi, mean = rep(0, Data$J), sigma=Data$tao, log = TRUE)
  
  # log-posterior likelihood
  eta <- Data$f0+Data$X %*% xi
  p <- exp(eta)/(1+exp(eta))
  LL <- sum(dbinom(x=Data$y, size = Data$n, prob=p, log = TRUE)) # log likelihood of Y|xi
  LP <- LL+sum(xi.prior) # joint log likelihood of (Y, xi)
  
  # output
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LL,
                   yhat=Data$y, parm=parm)
  return(Modelout)
}


# data: 
out_pred_laplace <- function(fpca_fit, df_new = df %>% filter(id==1 & sind_inx<=195) %>% select(-eta_i)){
  
  # put data into correct format
  ns <- as.vector(table(df_new$bin)) # number of observations
  hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% select(Y) %>% unlist()# number of success
  nf <- ns-hs # number of failure
  max_bin <- length(unique(df_new$bin)) # assume new skipped bins  
  df_new2 <- data.frame(bin = unique(df_new$bin), ns, hs, nf)  
  
  # into a list
  tao <- diag(fpca_fit$evalues)
  f0 <- fpca_fit$mu[1:max_bin]
  N <- nrow(df_new2) # "sample size" which is in fact number of observed bins in this case for a new subject
  n <- df_new2$ns # number of experiements at each bin
  y <- df_new2$hs # outcome, number of 1
  X <- fpca_fit$efunctions[1:max_bin, ] # "covariates", in this case PC functions
  J <- ncol(fpca_fit$efunctions) # number of parameters/scores
  mon.names <- "LP"
  parm.names <- as.parm.names(list(xi1=0, xi2=0, xi3=0, xi4=0))
  pos.xi <- grep("xi", parm.names)
  PGF <- function(Data) {
    xi <- rnorm(Data$J)
    return(xi)
  }
  
  MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
                 parm.names=parm.names, pos.xi=pos.xi, y=y, n=n, tao=tao, f0=f0)
  
  
  # fit laplace approximation
  Fit <- LaplaceApproximation(Model, parm = rep(0, 4), Data=MyData)
  score <- Fit$Summary1[, "Mode"]
  
  # prediction
  eta_pred_out <- fpca_fit$mu+fpca_fit$efunctions%*%score
  
  return(list(eta_pred = eta_pred_out,
               score_out = score))
}



# fit LaplacesDemon model
# Fit <- LaplacesDemon(Model, Data=MyData, rep(0, 4),
#                      Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#                      Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

# test

# out_pred_laplace(fpca_fit,  df %>% filter(id==2 & sind_inx<=395) %>% select(-eta_i))


