# This script codes the computaiton of evaluation of predictive performance

library(ROCR)
library(GLMMadaptive)


obs_t <- c(195, 395, 595)

ise_fgfpca <- rep(NA, 3)
ise_adglmm <- rep(NA, 3)

auc_fgfpca <- rep(NA, 3)
auc_adglmm <- rep(NA, 3)

comp_t_fgfpca <- rep(NA, 3)
comp_t_adglmm <- rep(NA, 3)

adglmm <- mixed_model(Y ~ sind, random = ~ sind | id, data =  df,
                      family = binomial())

##### Midpoint ISE ##### 

# true eta and true y
# head(df_true_eta)
# table(df_true_eta$id)
# unique(df_true_eta$mid)

eta_true_mat <- df_true_eta %>%
  pivot_wider(names_from = id, values_from = true) %>%
  column_to_rownames("mid")

y_mat <- df %>% select(id, Y, sind_inx) %>% 
  pivot_wider(names_from = id, values_from = Y) %>%
  column_to_rownames("sind_inx")

# observation up to 195

for(j in 1: length(obs_t)){
  
  t <- obs_t[j]
  pred_t <- which(mid > t) # unobserved bin index (binned grid)
  pred_t_full <- (t+1):max(mid) # unobserved time index (original grid)
  
  # fast GFPCA prediction container
  eta_hat_mat <- matrix(NA, nrow=length(mid), ncol = N)
  eta_hat_mat_full <- matrix(NA, nrow=J, ncol = N) # latent Gaussian function on original grid by interpolation
 
  # prediction
  tstart <- Sys.time()
  for(i in 1:N){
    eta_hat_mat[ , i] <- out_samp_dyn_pred(df_new = df %>% filter(id==i & sind_inx<=t) %>%
                                             select(-eta_i),
                                           fpca_fit = fpca_fit)$eta_pred
    eta_hat_mat_full[, i] <- approx(x=mid, y=eta_hat_mat[ , i], xout = 1:J)$y
  }
  tend <- Sys.time()
  
  ## ISE
  ise_fgfpca[j] <- mean(apply((eta_hat_mat[pred_t, ]-eta_true_mat[pred_t, ])^2, 2, sum))
  ## AUC after interpolation
  prob_mat <- exp(eta_hat_mat_full[pred_t_full, ])/(1+exp(eta_hat_mat_full[pred_t_full, ]))
  auc_perf <- performance(prediction(prob_mat, y_mat[pred_t_full,]), "auc")
  auc_fgfpca[j] <-mean(unlist(auc_perf@y.values))
  ## computation time
  comp_t_fgfpca[j] <- tend-tstart
  
  # reference method: adaptive GLMM
  # https://drizopoulos.github.io/GLMMadaptive/
  # fit model with full data (time as fixed effect, random intercept and slopt for time)
  # only evaluate SE at bin midpoint
  tstart2 <- Sys.time()
  preds_adglmm <- predict(adglmm, 
                       newdata = df %>% filter(sind_inx <= t),
                       newdata2 =  df %>% filter(sind_inx > t), 
                       type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  tend2 <- Sys.time()
  
  eta_hat_mat2 <- preds_adglmm$newdata2 %>% select(id, sind_inx, pred) %>% 
    filter(sind_inx %in% mid) %>%
    pivot_wider(names_from = id, values_from = pred) %>% 
    column_to_rownames("sind_inx") %>% as.matrix()
  eta_hat_mat_full2 <- preds_adglmm$newdata2 %>% select(id, sind_inx, pred) %>% 
    filter(sind_inx <= max(mid)) %>% 
    pivot_wider(names_from = id, values_from = pred) %>% 
    column_to_rownames("sind_inx") %>% as.matrix()

 ## ISE
 ise_adglmm[j] <- mean(apply((eta_hat_mat2-eta_true_mat[pred_t, ])^2, 2, sum)) 
 ## AUC after interpolation
 prob_mat2 <- exp(eta_hat_mat_full2)/(1+exp(eta_hat_mat_full2))
 auc_perf2 <- performance(prediction(prob_mat2, y_mat[pred_t_full, ]), "auc")
 auc_adglmm[j] <-mean(unlist(auc_perf2@y.values))
 # computation time
 comp_t_adglmm[j] <- tend2-tstart2

}




