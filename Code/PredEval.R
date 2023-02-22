# This script codes the computaiton of evaluation of predictive performance

library(ROCR)

##### Interpolated ISE and AUC #####

# true eta and true y
eta_true_mat <- df %>%
  select(id, sind_inx, eta_i) %>%
  pivot_wider(names_from = id, values_from = eta_i) %>%
  column_to_rownames("sind_inx")

# View(eta_true_mat)

y_mat <- df %>%
  select(id, sind_inx, Y) %>%
  pivot_wider(names_from = id, values_from = Y) %>%
  column_to_rownames("sind_inx")

# View(y_mat)

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


# Adaptive GLMM

library(GLMMadaptive)

View(df)
fm1 <- mixed_model(Y ~ sind, random = ~ sind | id, 
                   data =  df,
                   family = binomial())

preds_fm1 <- predict(fm1, 
                     newdata = df %>% filter(sind_inx <= 195),
                     newdata2 =  df %>% filter(sind_inx > 195), 
                     type = "subject_specific", type_pred = "link",
                     se.fit = TRUE, return_newdata = TRUE)
head(preds_fm1$newdata2)

ggplot(preds_fm1$newdata2%>% filter(id==1))+
  geom_line(aes(x=sind_inx, y=pred))
#large coefficient issue

# functional mixed models
# generalized function-on-scalar regression
# marginal constrants? 
fit_gam <- gam(Y ~ s(sind_inx, k=20, bs = "cr") + ti(id, sind_inx, bs=c("re","cr"), mc=c(TRUE,FALSE), k=c(5,5)), 
              family = "binomial",
              data =  df %>% select(id, Y, sind_inx) %>% filter(id %in% 1:100))



##### Midpoint ISE ##### 

# true eta and true y
head(df_true_eta)
table(df_true_eta$id)
unique(df_true_eta$mid)

eta_true_mat <- df_true_eta %>%
  pivot_wider(names_from = id, values_from = true) %>%
  column_to_rownames("mid")

y_mat <- df %>%
  select(id, sind_inx, Y) %>% 
  filter(sind_inx %in% mid) %>% 
  pivot_wider(names_from = id, values_from = Y) %>%
  column_to_rownames("sind_inx")

# observation up to 195
eta_hat_mat <- matrix(NA, nrow=length(mid), ncol = N) # prediction container
for(i in 1:N){
  eta_hat_mat[ , i] <- out_samp_dyn_pred(df_new = df %>% filter(id==i & sind_inx<=195) %>%
                      select(-eta_i),
                    fpca_fit = fpca_fit)$eta_pred
}

## do not include observed bins
pred_t <- which(mid > 195) # unobserved bin index

## ISE
ise_t195 <- mean(apply((eta_hat_mat[pred_t, ]-eta_true_mat[pred_t, ])^2, 2, sum))

## AUC
prob_mat <- exp(eta_hat_mat)/(1+exp(eta_hat_mat))
auc_perf <- performance(prediction(prob_mat[pred_t, ], y_mat[pred_t, ]), "auc")
mean(unlist(auc_perf@y.values))

## Reference model: 
### functional mixed models/generalized function-on-scalar regression
# marginal constrants? 
fit_gam <- gam(Y ~ s(sind_inx, k=20, bs = "cr") + ti(id, sind_inx, bs=c("re","cr"), mc=c(TRUE,FALSE), k=c(5,5)), 
               family = "binomial",
               data =  df %>%
                 select(id, sind_inx, Y) %>% 
                 filter(sind_inx %in% mid))



