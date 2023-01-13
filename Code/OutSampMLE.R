
library(rootSolve)

# This script writes functions for 
# calculating MLE solution of fPCA scores 
# based on conditonal distribution of score (conditional on observations)
# of an out-of-sample observation
# just for mess around


#### MLE of out-of-sample score ####

# First derivative of Log-likelihood
# df_new: new observations, including at least columns Y and bin
# fpca_fit: fpca object


llh_div <- function(xi, df_new, fpca_fit){
  
  # observations
  ns <- as.vector(table(df_new$bin)) # number of observations
  hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% 
    select(Y) %>% unlist()# number of success
  nf <- ns-hs # number of failure
  max_bin <- length(unique(df_new$bin)) # assume new skipped bins 
  
  # fpca objects
  phi_m <- fpca_fit$efunctions[1:max_bin, ]
  f0_m <- fpca_fit$mu[1:max_bin]
  eta_s <- f0_m + phi_m %*% xi
  tao <- diag(fpca_fit$evalues)
  
  # derivative
  val1 <- hs-ns*(exp(eta_s)/(1+exp(eta_s)))
  div <- colSums(apply(phi_m, 2, function(x)x*val1))-xi%*%solve(tao)
  return(div)
}

# test
# llh_div(rep(0, 4), df_new=df_new, fpca_fit = fpca_fit)

#### Prediction with MLE scores #####

out_samp_dyn_pred <- function(df_new, fpca_fit){
  score_out <- multiroot(f=llh_div, start = rep(0, 4), df_new=df_new, fpca_fit=fpca_fit)$root
  eta_pred_out <- fpca_fit$mu+fpca_fit$efunctions%*%score_out
  
  return(list(eta_pred = eta_pred_out,
              score_out = score_out))
  
}


# test
# out_samp_dyn_pred(df_new, fpca_fit)



