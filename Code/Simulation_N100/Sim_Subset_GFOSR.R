set.seed(1114)

library(here)
library(tidyverse)
library(ggplot2)
theme_set(theme_minimal())
library(fcr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(mgcv)
library(ggpubr)
library(gridExtra)
library(kableExtra)
library(splines)




#### load and reduce simulated data ####

load(here("Data/sim_data.RData"))

# reduce data
## training data
M <- 500 # number of simulations

sim_data <- lapply(sim_data[1:M], function(x){x %>% filter(id %in% c(1:200))})

lapply(sim_data, dim)
head(sim_data[[1]])
unique(sim_data[[1]]$id)
unique(sim_data[[1]]$t)

Ntr <- Nte <- 100 # trainig and testing sample size
N <- Ntr+Nte
J <- 1000

t <- unique(sim_data[[1]]$t)

K <- 4 # number of PCs


# overview of a few samples
rand_id <- sample(N, size = 4)
sim_data[[357]] %>% filter(id %in% rand_id) %>% 
  ggplot()+
  geom_point(aes(x=t, y=Y))+
  geom_line(aes(x=t, y=plogis(eta_i)), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)

#### model set up ####

L <- 1
# L <- 5 # number of last observations used as time-fixed predictor
windows <- seq(0, 1, by = 0.2) # prediction window


#### Simulation ####
t_vec_gfofr_subset <- rep(NA, M)
pred_list_gfofr_subset <- list()

pb <- txtProgressBar(min=0, max = M, initial = 0, style = 3)

for(m in 1:M){
  df_m <- sim_data[[m]]
  df_m$id <- droplevels(df_m$id)
  
  # prediction window 
  df_m$window <- cut(df_m$t, breaks=windows, labels = 1:5, include.lowest = T)
  
  # split data
  train_df <- df_m %>% filter(id %in% 1:Ntr)
  test_df <- df_m %>% filter(!id %in% 1:Ntr)
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
  
  t_vec_gfofr_subset[m] <- t2-t1 # seconds
  pred_list_gfofr_subset[[m]] <- test_df
  
  setTxtProgressBar(pb, m)
}

close(pb)

save(pred_list_gfofr_subset, t_vec_gfofr_subset, 
     file = here("Data/SubSimOutput_GFOSR_L1.RData"))

pred_list_gfofr_subset[[375]] %>% 
  filter(t>0.55 & t<0.65) %>% View()
  
pred_list_gfofr_subset[[396]] %>% filter(id %in% c(101, 111, 121, 131)) %>% 
    ggplot()+
    geom_point(aes(x=t, y=Y))+
    geom_line(aes(x=t, y=plogis(eta_i), col = "red"))+
    geom_line(aes(x=t, y=plogis(pred_w1), col = "w1"))+
    geom_line(aes(x=t, y=plogis(pred_w2), col = "w2"))+
    geom_line(aes(x=t, y=plogis(pred_w3), col = "w3"))+
    geom_line(aes(x=t, y=plogis(pred_w4), col = "w4"))+
    facet_wrap(~id)
