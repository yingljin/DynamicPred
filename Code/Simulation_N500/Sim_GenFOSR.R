# This script implements the full simulation
# with a second reference method
# a generalized function-on-funciton regression
# on the simulated dataset generated from Code/GenData.R

#### set up ####

set.seed(1114)

library(here)
library(tidyverse)
library(ggpubr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(kableExtra)
library(knitr)
library(mvtnorm)
library(mgcv)
library(splines)
library(LaplacesDemon)
theme_set(theme_minimal())


#### load simulated data ####

load(here("Data/sim_data.RData"))

N <- length(unique(sim_data[[1]]$id))
J <- length(unique(sim_data[[1]]$t))
t <- unique(sim_data[[1]]$t)
M <- length(sim_data)
K <- 4 # number of eigenfunctions to use

# overview of a few samples
rand_id <- sample(N, size = 4)
sim_data[[357]] %>% filter(id %in% rand_id) %>% 
  ggplot()+
  geom_point(aes(x=t, y=Y))+
  geom_line(aes(x=t, y=plogis(eta_i)), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)


N_train <- 500 # training sample size
N_test <- 100 # testing sample szie 
# L <- 5 # number of last observations used as time-fixed predictor
L <- 1
windows <- seq(0, 1, by = 0.2) # prediction window

#### Simulation ####
t_vec_gfofr <- rep(NA, M)
pred_list_gfofr <- list()

pb <- txtProgressBar(min=0, max = M, initial = 0, style = 3)

for(m in 1:M){
  df_m <- sim_data[[m]]
  
  # prediction window 
  df_m$window <- cut(df_m$t, breaks=windows, labels = 1:5, include.lowest = T)
  
  # split data
  train_df <- df_m %>% filter(id %in% 1:N_train)
  test_df <- df_m %>% filter(!id %in% 1:N_train)
  train_df$window <- as.numeric(as.character(train_df$window))
  test_df$window <- as.numeric(as.character(test_df$window))
  
  t1 <- Sys.time()
  for(w in 1:4){
    # model fit
    # predictor (training set)
    y_obs_max <- train_df %>% filter(window==w) %>% group_by(id) %>%
      slice_max(t, n=L) %>% select(id, Y, sind) %>%
      mutate(name = sind-min(sind)+1) %>%
      pivot_wider(id_cols = id, values_from = Y, names_from = "name", names_prefix = "yl")
    # outcome (training set)
    df_pred_tr <- train_df %>% filter(window > w) %>%
      left_join(y_obs_max, by = "id") %>% 
      mutate_at(vars(Y, starts_with("yl")), as.factor) 
    # fit_gen_fosr <- bam(Y ~ s(t, bs="cr", k=20) + 
    #                       s(t, bs="cr", k=20, by = yl5)+
    #                       s(t, bs="cr", k=20, by = yl4)+
    #                       s(t, bs="cr", k=20, by = yl3)+
    #                       s(t, bs="cr", k=20, by = yl2)+
    #                       s(t, bs="cr", k=20, by = yl1),
    #                     family = binomial, data=df_pred_tr, 
    #                     method = "fREML",
    #                     discrete = TRUE)
    fit_gen_fosr <- bam(Y ~ s(t, bs="cr", k=20) + 
                          s(t, bs="cr", k=20, by = yl1),
                        family = binomial, data=df_pred_tr, 
                        method = "fREML",
                        discrete = TRUE)
    
    # prediction
    # predictor (testing set)
    y_obs_max_te <- test_df %>% filter(window==w) %>% group_by(id) %>%
      slice_max(t, n=L) %>% select(id, Y, sind) %>% 
      mutate(name = sind-min(sind)+1) %>%
      pivot_wider(id_cols = id, values_from = Y, names_from = name, names_prefix = "yl")
    # outcome (testing set)
    df_pred_te <- test_df %>% filter(window > w) %>%
      left_join(y_obs_max_te, by = "id") %>% 
      mutate_at(vars(Y, starts_with("yl")), as.factor)
    
    # predict using the FOSR model fit above
    pred_name <- paste0("pred_w", w)
    test_df[, pred_name] <- NA
    test_df[test_df$window>w, pred_name] <- predict(fit_gen_fosr, newdata = df_pred_te, type = "link")
  }
  t2 <- Sys.time()
  
  t_vec_gfofr[m] <- t2-t1 # seconds
  pred_list_gfofr[[m]] <- test_df
  
  setTxtProgressBar(pb, m)
}

close(pb)

save(pred_list_gfofr, t_vec_gfofr, file = here("Data/SimOutput_GFOSR_L1.RData"))

pred_list_gfofr[[399]] %>%
  filter(t>0.55 & t<0.65)
  View()
  
pred_list_gfofr[[399]] %>% filter(id %in% c(501, 511, 521, 531)) %>% 
    ggplot()+
    geom_point(aes(x=t, y=Y))+
    geom_line(aes(x=t, y=plogis(eta_i), col = "red"))+
    geom_line(aes(x=t, y=plogis(pred_w1), col = "w1"))+
    geom_line(aes(x=t, y=plogis(pred_w2), col = "w2"))+
    geom_line(aes(x=t, y=plogis(pred_w3), col = "w3"))+
    geom_line(aes(x=t, y=plogis(pred_w4), col = "w4"))+
    facet_wrap(~id)

# test_df %>% filter(sind > 400)
