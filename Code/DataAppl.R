# This script runs the data application section

set.seed(330)

library(here)
library(here)
library(tidyverse)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)

# code
source(here("Code/GLMM-FPCA.R")) # use pred_latent function to estimate latent function 
source(here("Code/OutSampMLE.R"))

# data
df <- read.csv("Data/Daily_Weigh_ins.csv")

#### Data clean ####

J <- ncol(df %>% select(starts_with("day")))  # number of observations points
N <- nrow(df) # sample size
rand_id <- sample(N, 4) # sample used to visualization


# move from wide to long
df_long <- df %>% select(participant_id, starts_with("day")) %>% 
  pivot_longer(2:401, values_to = "Y", names_to = "sind_inx") %>%
  mutate(sind_inx = gsub("day_", "", sind_inx)) %>%
  mutate(sind_inx = as.numeric(sind_inx)) %>%
  rename("id" = "participant_id")

df_long %>%
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind_inx, y=Y))+
  facet_wrap(~id)


##### fGFPCA #####

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

df_long$bin <- cut(df_long$sind_inx, breaks = brks, include.lowest = T, labels = mid)
df_long$bin <- as.numeric(as.character(df_long$bin))

df_bin_lst <- split(df_long, f = df_long$bin)

lapply(df_bin_lst, function(x)table(x$Y))

# fit local GLMM and estimate latent function
# near-unidentifiability issues !!!
df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
df_est_latent <- bind_rows(df_est_latent) %>% select(id, bin, eta_hat) %>% distinct() # on the binned grid, one unique value each bin
df_est_latent %>%
  left_join(df_long, by = c("id", "bin")) %>%
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_line(aes(x=bin, y= exp(eta_hat)/(1+exp(eta_hat))))+
  geom_point(aes(x=bin, y = Y))+
  facet_wrap(~id)

# fPCA
mat_est_unique <- matrix(df_est_latent$eta_hat, nrow=N, ncol=n_bin, byrow = F) # row index subject, column binned time
fpca_mod <- fpca.face(mat_est_unique, pve = 0.95, argvals = mid, knots=20, var=T)
K <- ncol(fpca_mod$efunctions) # number of eigenfunctions

# Out-of-sample prediction
df_est_latent[, 'pred_t100'] <- df_est_latent[, 'pred_t200'] <- df_est_latent[, 'pred_t300'] <- NA 
score_out_mat <- array(NA, dim = c(N, K, 3)) # dim indexes subject, eigenfunction and max obs time respectively

# prediction for a single subject
for(i in 1:N){
  df_i <- df_long %>% filter(id==i) 
  # prediction 
  pred_t100 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=100), fpca_fit = fpca_mod, K = K)
  pred_t200 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=200), fpca_fit = fpca_mod, K = K)
  pred_t300 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=300), fpca_fit = fpca_mod, K = K)
  # prediction in container
  df_est_latent[df_est_latent$id==i, 'pred_t100'] <- pred_t100$eta_pred
  df_est_latent[df_est_latent$id==i, 'pred_t200'] <- pred_t200$eta_pred
  df_est_latent[df_est_latent$id==i, 'pred_t300'] <- pred_t300$eta_pred
  # score in container
  score_out_mat[i, ,1] <- pred_t100$score_out
  score_out_mat[i, ,2] <- pred_t200$score_out
  score_out_mat[i, ,3] <- pred_t300$score_out
}

save(df_est_latent, score_out_mat, file = here("Data/ApplOutput_fGFPCA.RData"))


##### GLMMadaptive #####

adglmm <- mixed_model(Y ~ sind_inx, random = ~ sind_inx | id, 
                      data =  df_long, family = binomial())
  
# container
df_est_latent_ref <- df_long
df_est_latent_ref[, 'pred_t100'] <- df_est_latent_ref[, 'pred_t200'] <- df_est_latent_ref[, 'pred_t300'] <- NA 
  
# time up to t
pred_t100 <- predict(adglmm, newdata = df_long %>% filter(sind_inx <= 100),
                     newdata2 =  df_long %>% filter(sind_inx > 100), type = "subject_specific", 
                     type_pred = "link", se.fit = TRUE, return_newdata = TRUE)
df_est_latent_ref$pred_t100[df_est_latent_ref$sind_inx>100] <- pred_t100$newdata2$pred

pred_t200 <- predict(adglmm, newdata = df_long %>% filter(sind_inx <= 200),
                     newdata2 =  df_long %>% filter(sind_inx > 200), type = "subject_specific", 
                     type_pred = "link", se.fit = TRUE, return_newdata = TRUE)
df_est_latent_ref$pred_t200[df_est_latent_ref$sind_inx>200] <- pred_t200$newdata2$pred

pred_t300 <- predict(adglmm, newdata = df_long %>% filter(sind_inx <= 300),
                     newdata2 =  df_long %>% filter(sind_inx > 300), type = "subject_specific", 
                     type_pred = "link", se.fit = TRUE, return_newdata = TRUE)
df_est_latent_ref$pred_t300[df_est_latent_ref$sind_inx>300] <- pred_t300$newdata2$pred



save(df_est_latent_ref, file = here("Data/ApplOutput_GLMMadaptive.RData"))
