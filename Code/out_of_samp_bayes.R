library(cubature)
# This script writes functions for 
# calculating conditional expectation
# of fPCA scores 
# of an out-of-sample observation

# first, only subject one, time up to 395

#### Denominator of E(xi|Y) ####
# xi: scores. The value to integrate over
# tao: estimated covariance matrix(eigenvalues) of score from fPCA
# f0: mean function from fPCA (up to mth bin)
# phi: basis function of fPCA (up to mth bin)
# Y: newly observed data (up to mth bin)

# just for mess around
tao <- diag(fpca_fit$evalues)
f0 <- fpca_fit$mu[1:max_bin]
phi <- fpca_fit$efunctions[1:max_bin, ] 
df_out <- df %>% filter(id == 1 & sind_inx <= 395) %>% select(-eta_i)
max_bin <- which(mid==395)
xi <- fpca_fit$scores[1, ]
xi <- matrix(xi, nrow = 4, ncol = 1)
k <- 4
  
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
joint_pdf(xi, tao, f0, phi, df_out)

adaptIntegrate(joint_pdf, lower = rep(-500, K), upper = rep(500, K), 
          tao=tao, f0=f0, phi=phi, df_out=df_out)
  
  