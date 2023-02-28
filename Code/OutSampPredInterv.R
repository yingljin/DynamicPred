# This Script calculates variance of prediction latent function


#### EB Guassian ####

# same as in-sample

pred_std_eb <- function(fpca_fit, tm){
  m <- which(mid==tm)  # maximum observation time on bin label scale
  tao <- diag(fpca_fit$evalues) ## egienvalues, also var-cov matrix of scores
  Phi <- fpca_fit$efunctions[1:m, ] ## eigenfunctions on the entire binned grid
  
  H <-  tao %*% t(Phi) %*% solve(Phi%*%tao%*%t(Phi)+fpca_fit$sigma2*diag(m))
  eta_pred_var <- fpca_fit$efunctions %*% (tao - H %*% Phi %*% tao) %*% t(fpca_fit$efunctions) + fpca_fit$sigma2
  
  std <- sqrt(diag(eta_pred_var))
  
  return(std)
  
}

# pred_std(fpca_fit, 195)[61:65]


#### MLE ####


