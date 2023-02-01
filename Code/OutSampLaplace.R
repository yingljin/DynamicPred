
# This script trys to approximate the integrals involved in posterior mean of scores
# using laplace approximation

library(mvtnorm)
library(rootSolve)


##### Marginal distribution of outcome #####

# first, approximate the denominator
# We need: 
# joint log-likelihood of outcome and score
# MLE of score
# second derivative of joint log-likelihood wrt scores

# subject 1 with observations up to t=195

df_new <- df %>% filter(id==1 & sind_inx<=195) %>% select(-eta_i)

ns <- as.vector(table(df_new$bin)) # number of observations
hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% select(Y) %>% unlist()# number of success
nf <- ns-hs # number of failure
max_bin <- length(unique(df_new$bin)) # assume new skipped bins  

# additional information needed from FPCA 
tao <- diag(fpca_fit$evalues)
f0 <- fpca_fit$mu[1:max_bin]
Phi <- fpca_fit$efunctions[1:max_bin, ]
K <- ncol(Phi)

# function for joint likelihood

joint_llh <- function(xi){
  eta <- f0+ Phi%*%xi
  p <- exp(eta)/(1+exp(eta))
  pos_y_llh <- sum(dbinom(x=hs, size = ns, prob=p, log=T)) # posterior log-likelihood of outcome
  pri_xi_llh <- dmvnorm(xi, mean = rep(0, K), sigma=tao, log = T) # prior log-likelihood of scores
  return(pos_y_llh+pri_xi_llh)
}

# MLE
maximum <- optim(par=rep(0, K), fn = joint_llh, control = list(fnscale=-1))
mle_score <-maximum$par
mle_max <- maximum$value

# 2nd derivative
joint_llh_div2 <- function(xi = mle_score){
  eta <- f0+ Phi%*%xi
  div <- 0
  for(i in 1:max_bin){
    div <- div+ns[i]*exp(eta[i])/(1+exp(eta[i]))^2*crossprod(t(Phi[i, ]))
  }
  
  return(div)
}

# approximate marginal likelihood of Y
den <- exp(mle_max+log((2*pi)^K/det(joint_llh_div2()))/2)


##### joint expectation of scores #####

# then, approximate the numerator
# We need: log score, posterior log likelihood of observation, prior marginal log likelihood of score
# but the target function can be negative, thus cannot be approximated by Laplace's method

# target function of laplace approximation
tfunc <- function(xi=mle_score){
  return(log(xi)+joint_llh(xi))
}

# MLE
maximum2 <- optim(par=mle_score, fn = tfunc, control = list(fnscale=-1))
mle_score <-maximum$par
mle_max <- maximum$value

##### joint likelihood of scores and outcome #####




    