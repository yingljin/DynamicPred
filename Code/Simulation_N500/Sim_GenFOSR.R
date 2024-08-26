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
  geom_line(aes(x=t, y=exp(eta_i)/(1+exp(eta_i))), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)


N_train <- 500 # training sample size
N_test <- 100 # testing sample szie 
L <- 1 # number of last observations used as time-fixed predictor
windows <- seq(0, 1, by = 0.2) # prediction window

#### Simulation ####
# M <-5
time_gfofr <- rep(NA, M)
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
    ## and prediction interval
    pred_mw <- predict(fit_gen_fosr, newdata = df_pred_te, 
                       type = "link",
                       se.fit = T)
    pred_int_lb <- pred_mw$fit-qnorm(0.975)*pred_mw$se.fit
    pred_int_ub <- pred_mw$fit+qnorm(0.975)*pred_mw$se.fit
    
    pred_name <- paste0(c("pred", "pred_lb", "pred_ub"), w)
    test_df[test_df$window>w, pred_name[1]] <- pred_mw$fit
    test_df[test_df$window>w, pred_name[2]] <- pred_int_lb
    test_df[test_df$window>w, pred_name[3]] <- pred_int_ub
  }
  t2 <- Sys.time()
  
  time_gfofr[m] <- difftime(t2, t1, units = "mins")
  pred_list_gfofr[[m]] <- test_df
  
  setTxtProgressBar(pb, m)
}

close(pb)

#### Check output ####
pred_list_gfofr[[1]] %>% filter(t>=0.35) %>% View()

rand_id <- sample(test_df$id, 4)

df_exp <- pred_list_gfofr[[3]] %>% 
  filter(id %in% rand_id) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  
df_exp[df_exp$t<=0.2, c("pred1", "pred_lb1", "pred_ub1")] <- NA
df_exp[df_exp$t<=0.4, c("pred2", "pred_lb2", "pred_ub2")] <- NA
df_exp[df_exp$t<=0.6, c("pred3", "pred_lb3", "pred_ub3")] <- NA
df_exp[df_exp$t<=0.8, c("pred1", "pred_lb4", "pred_ub4")] <- NA

df_exp %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred1, col = "0.2"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb1, ymax = pred_ub1,
                  col = "0.2", fill = "0.2", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred2, col = "0.4"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb2, ymax = pred_ub2,
                  col = "0.4", fill = "0.4", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred3, col = "0.6"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb3, ymax = pred_ub3,
                  col = "0.6", fill = "0.6", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred4, col = "0.8"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb4, ymax = pred_ub4,
                  col = "0.8", fill = "0.8", alpha = 0.1),
              linetype="dashed")+
  facet_grid(rows = vars(id))+
  guides(alpha = "none", col="none")+
  labs(title = "GFOSR (L=5)")
ggsave(here("Images/IntervalExp1_3.jpeg"), height=12, width = 5)  

#### Save results ####

pred_list_gfofr_l1 <- pred_list_gfofr
time_gfofr_l1 <- time_gfofr


save(pred_list_gfofr, time_gfofr, 
     file = here("Data//TrialRun/SimOutput_GFOSR_L1.RData"))
save(pred_list_gfofr_l1, time_gfofr_l1, 
     file = here("Data/SimN500/SimOutput_GFOSR_L1.RData"))

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
