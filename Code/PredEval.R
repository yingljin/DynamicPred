# This script codes the computaiton of evaluation of predictive performance

library(ROCR)

##### Interpolated ISE and AUC #####

# true eta and true y
eta_true_mat <- df %>%
  select(id, sind_inx, eta_i) %>%
  pivot_wider(names_from = id, values_from = eta_i) %>%
  column_to_rownames("sind_inx")
y_mat <- df %>%
  select(id, sind_inx, Y) %>%
  pivot_wider(names_from = id, values_from = Y) %>%
  column_to_rownames("sind_inx")

##### Observations up to 195 #####
pred_t <- 196:max(mid)
eta_hat_mat <- matrix(NA, nrow=length(pred_t), ncol = N) # prediction container
dim(eta_hat_mat)
for(i in 1:N){
  pred <- out_samp_dyn_pred(df_new = df %>% 
                              filter(id==i & sind_inx<=195) %>%
                              select(-eta_i),
                            fpca_fit = fpca_fit)$eta_pred
  interp_pred <- approx(x=mid, y=pred, xout = pred_t)$y
  eta_hat_mat[, i] <- interp_pred
}
# ISE
ise_t195 <- mean(apply((eta_hat_mat-eta_true_mat[pred_t, ])^2, 2, sum))
# AUC
prob_mat <- exp(eta_hat_mat)/(1+exp(eta_hat_mat))
auc_perf <- performance(prediction(prob_mat, y_mat[pred_t, ]), "auc")
mean(unlist(auc_perf@y.values))

for(i in 1:5){
  auc_perf <- performance(prediction(prob_mat[, i], y_mat[pred_t, i]), "auc")
  print(auc_perf@y.values[[1]])
}
summary(auc_perf)



##### Reference methods #####

# binary FPCA for dynamic prediction
# non-convergence
library(registr)

bfpca_fit <- bfpca(Y= df %>% select(id, Y, sind_inx) %>% rename(index=sind_inx, value = Y), 
      npc_varExplained = 0.95, t_min=1, t_max = J)

# A standard mxied model
# time as the only fixed effect
# random intercept for ID and random slope for time
# slow, also with numeric issues
glmm_fit <- glmer(Y ~ sind_inx + (id | sind_inx),
                  data = df %>% select(id, Y, sind_inx), family = binomial)

# functional mixed models
# generalized function-on-scalar regression
# marginal constrants? 
fit_gam <- gam(Y ~ s(sind_inx, k=20, bs = "cr") + ti(id, sind_inx, bs=c("re","cr"), mc=c(TRUE,FALSE), k=c(5,5)), 
              family = "binomial",
              data =  df %>% select(id, Y, sind_inx))



##### Midpoint ISE ##### 

# not include AUC because it is difficult to determin "midpoint"

# only true latent functions at bin midpoints are used
eta_true_mat <- df_true_eta %>% pivot_wider(names_from = id, values_from = true) %>%
  column_to_rownames("mid")
eta_true_mat <- as.matrix(eta_true_mat)

# observation up to 195
eta_hat_mat <- matrix(NA, nrow=length(mid), ncol = N) # prediction container
for(i in 1:N){
  eta_hat_mat[ , i] <- out_samp_dyn_pred(df_new = df %>% filter(id==i & sind_inx<=195) %>%
                      select(-eta_i),
                    fpca_fit = fpca_fit)$eta_pred
}
## do not include observed bins
pred_rows <- which(mid > 195) # unobserved bin index
eta_hat_mat[pred_rows,]
## ISE
ise_t195 <- mean(apply((eta_hat_mat[pred_rows, ]-eta_true_mat[pred_rows, ])^2, 2, sum))

# observation up to 395
eta_hat_mat <- matrix(NA, nrow=length(mid), ncol = N) # prediction container
for(i in 1:N){
  eta_hat_mat[ , i] <- out_samp_dyn_pred(df_new = df %>% filter(id==i & sind_inx<=395) %>%
                                           select(-eta_i),
                                         fpca_fit = fpca_fit)$eta_pred
}
## do not include observed bins
pred_rows <- which(mid > 395) # unobserved bin index
eta_hat_mat[pred_rows,]
ise_t395 <- mean(apply((eta_hat_mat[pred_rows, ]-eta_true_mat[pred_rows, ])^2, 2, sum))

# observation up to 595
eta_hat_mat <- matrix(NA, nrow=length(mid), ncol = N) # prediction container
for(i in 1:N){
  eta_hat_mat[ , i] <- out_samp_dyn_pred(df_new = df %>% filter(id==i & sind_inx<=595) %>%
                                           select(-eta_i),
                                         fpca_fit = fpca_fit)$eta_pred
}
## do not include observed bins
pred_rows <- which(mid > 595) # unobserved bin index
eta_hat_mat[pred_rows,]
ise_t395 <- mean(apply((eta_hat_mat[pred_rows, ]-eta_true_mat[pred_rows, ])^2, 2, sum))







