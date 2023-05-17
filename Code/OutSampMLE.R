
library(rootSolve)

# This script writes functions for 
# calculating MLE solution of fPCA scores 
# based on conditonal distribution of score (conditional on observations)
# of an out-of-sample observation
# just for mess around

#### Likelihood function ####
joint_pdf <- function(xi, df_new, fpca_fit){
  
  # observations 
  ns <- as.vector(table(df_new$bin)) # number of observations
  hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% 
    select(Y) %>% unlist()# number of success
  nf <- ns-hs # number of failure
  max_bin <- length(unique(df_new$bin)) # assume no skipped bins 
  
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
  max_bin <- length(unique(df_new$bin)) # assume no skipped bins 
  
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

out_samp_dyn_pred <- function(df_new, fpca_fit, K){
  score_out <- multiroot(f=llh_div, start = rep(0, K), df_new=df_new, fpca_fit=fpca_fit)$root
  eta_pred_out <- fpca_fit$mu+fpca_fit$efunctions%*%score_out
  
  
  return(list(eta_pred = eta_pred_out,
              score_out = score_out))
}

# test
# multiroot(f=llh_div, start = rep(0, 4), df_new=df_i %>% filter(sind_inx<=800), fpca_fit=fpca_mod)

#### Estimate standard error of MLE ####

# # Jacabian matrix of first dirivative of llh
# try <- out_samp_dyn_pred(df_new = df %>% filter(id==1 & sind_inx<=195) %>%
#                     select(-eta_i),
#                   fpca_fit = fpca_fit)$score_out
# 
# I_mat <- -gradient(f = llh_div, x = try,
#          df_new = df %>% filter(id==1 & sind_inx<=195) %>%
#            select(-eta_i),
#          fpca_fit=fpca_fit)
# 
# var <- diag(solve(I_mat))
# sqrt(diag(fpca_fit$efunctions %*% solve(I_mat) %*% t(fpca_fit$efunctions)))

# I_mat <- function(xi, df_new, fpca_fit){
#   
#   # observations
#   ns <- as.vector(table(df_new$bin)) # number of observations
#   hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% 
#     select(Y) %>% unlist()# number of success
#   max_bin <- length(unique(df_new$bin)) # assume new skipped bins 
#   
#   # fpca objects
#   phi_m <- fpca_fit$efunctions[1:max_bin, ]
#   f0_m <- fpca_fit$mu[1:max_bin]
#   eta_s <- f0_m + phi_m %*% xi
#   tao <- diag(fpca_fit$evalues)
#   
#   # derivative
#   val1 <- ns*(exp(eta_s)/(1+exp(eta_s))^2)
#   div2 <- lapply(1:max_bin,  function(x)val1[x]*phi_m[x, ] %*% t(phi_m[x, ]))
#   div2 <- array(unlist(div2), dim = c(K, K, max_bin))
#   div2 <- apply(div2, 1:2, sum)+solve(tao)
#   
#   return(div2)
# }
# 
# var <- diag(solve(I_mat))
# sqrt(var)
# # # test
# # out_samp_dyn_pred(df_new, fpca_fit)
# 
# I_mat(try, df %>% filter(id==1 & sind_inx<=195) %>%
#         select(-eta_i), fpca_fit)


