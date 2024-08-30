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
L <- 5
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
    fit_gen_fosr <- bam(Y ~ s(t, bs="cr", k=20) +
                          s(t, bs="cr", k=20, by = yl5)+
                          s(t, bs="cr", k=20, by = yl4)+
                          s(t, bs="cr", k=20, by = yl3)+
                          s(t, bs="cr", k=20, by = yl2)+
                          s(t, bs="cr", k=20, by = yl1),
                        family = binomial, data=df_pred_tr,
                        method = "fREML",
                        discrete = TRUE)
    # fit_gen_fosr <- bam(Y ~ s(t, bs="cr", k=20) +
    #                       s(t, bs="cr", k=20, by = yl1),
    #                     family = binomial, data=df_pred_tr,
    #                     method = "fREML",
    #                     discrete = TRUE)
    
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
  
  t_vec_gfofr_subset[m] <- difftime(t2, t1, units = "mins")
  pred_list_gfofr_subset[[m]] <- test_df
  
  setTxtProgressBar(pb, m)
}

close(pb)

#### check result ####

pred_list_gfofr_subset[[1]] %>% filter(t>=0.35) %>% View()
pred_list_gfofr_subset[[1]] %>% filter(id %in% 101:104) %>%
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  %>%
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
  labs(title = "GFOSR (L=1)")

mean(t_vec_gfofr_subset)

#### save results ####

pred_subset_gfofr_l5 <- pred_list_gfofr_subset
time_subset_gfofr_l5 <- t_vec_gfofr_subset

save(pred_subset_gfofr_l5, time_subset_gfofr_l5, 
     file = here("Data/SimN100/SubSimOutput_GFOSR_L5.RData"))
