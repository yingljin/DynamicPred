
##### In-sample prediction on binned grid #####
# fpca_fit: fpca model fit on estimated latent Gaussian function
# df_part: data frame including partial observations, with at least columns "bin", "eta_hat"


# return: predictions on the full binned grid

in_samp_dyn_pred <- function(fpca_fit, df_part){
  # values for prediction
  tm <- length(unique(df_part$bin)) # maximum observation time on bin label scale
  tao <- diag(fpca_fit$evalues) ## egienvalues, also var-cov matrix of scores
  Phi <- fpca_fit$efunctions[1:tm, ] ## eigenfunctions on the entire binned grid
  
  
  # re-estimate score
  H <-  tao %*% t(Phi) %*% solve(Phi%*%tao%*%t(Phi)+fpca_fit$sigma2*diag(tm))
  score_hat <- H %*% (df_part$eta_hat-fpca_fit$mu[1:tm])
  
  # predict latent function
  eta_pred <- fpca_fit$mu+fpca_fit$efunctions%*%score_hat
  
  # predict std
  eta_pred_var <- fpca_fit$efunctions %*% (tao - H %*% Phi %*% tao) %*% t(fpca_fit$efunctions) + fpca_fit$sigma2
  
  return(list(eta_pred = eta_pred, 
              eta_pred_var = eta_pred_var,
              eta_std_var = sqrt(diag(eta_pred_var)),
              score_hat = score_hat))
}

# return: predicted latent function, its variance/std, PC scores

# try <- in_samp_dyn_pred(fpca_fit, 60, df_pred_unique, 4)
# 
# try$eta_std_var
# 
# plot(1:length(try$eta_std_var), try$eta_std_var)
# fpca_fit$diag.var
# dim(fpca_fit$VarMats[[2]])
