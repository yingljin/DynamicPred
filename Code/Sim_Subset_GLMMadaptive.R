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




#### load simulated data ####

load(here("Data/sim_data.RData"))

# reduce data
## training data
M <- 100 # number of simulations

sim_data <- lapply(sim_data[1:M], function(x){x %>% filter(id %in% c(1:200))})
lapply(sim_data, dim)
head(sim_data[[1]])
unique(sim_data[[1]]$id)
unique(sim_data[[1]]$t)

Ntr <- Nte <- 100 # trainig and testing sample size
J <- 1000


#### Data reduction set up ####

t_bin <- seq(0, J, by = 10)

## model fit
## in training, take only 1 observation for every 10 observation
## for example, only the bin midpoints
## It seems that bin width does not affect computational speed too much
## at least much less than model flexibility

#### Containers ####


pred_list_all <- list()
fit_time <- pred_time <- rep(NA, M)

#### GLMMadaptive ####

# The simulation process run into a numeric problem at the 95th iteration
# (singularity)

pb = txtProgressBar(min = 95, max = M, initial = 0, style = 3) 
for(m in 95:M){
  
  # data split and training data reduction
  try_fit <- NA
  df_test_m <- sim_data[[m]] %>% filter(!id %in% 1:100)
  
  # model fit
  while(is.na(try_fit)){
    df_train_m <- sim_data[[m]] %>% filter(id %in% 1:100 & sind %in% t_bin) 
    try_fit <- tryCatch(expr={t1 <- Sys.time()
                              fit_adglmm <- mixed_model(fixed = Y ~ 0+ns(t, df = 4), 
                                                        random = ~ 0+ns(t, df = 4) | id, 
                                                        data =  df_train_m, family = binomial())
                              t2 <- Sys.time()},
                        error = function(e){NA}
                        )
    t_bin <- t_bin+1
    
  }
  
  fit_time[m] <- t2-t1 ## minutes

  
  # prediction
  t1 <- Sys.time()
  pred_m <- list(
    predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.2),
            newdata2 =  df_test_m %>% filter(t>0.2), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.4),
            newdata2 =  df_test_m %>% filter(t>0.4), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.6),
            newdata2 =  df_test_m %>% filter(t>0.6), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.8),
            newdata2 =  df_test_m %>% filter(t>0.8), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2)
  t2 <- Sys.time()
  pred_time[m] <- t2-t1 # seconds
  
  pred_m <- lapply(pred_m, function(x){x %>% select(id, t, pred)})
  df_pred_m <- df_test_m %>% 
    left_join(pred_m[[1]]) %>% 
    rename(pred0.2 = pred) %>% 
    left_join(pred_m[[2]]) %>% 
    rename(pred0.4 = pred) %>% 
    left_join(pred_m[[3]]) %>% 
    rename(pred0.6 = pred) %>% 
    left_join(pred_m[[4]]) %>% 
    rename(pred0.8 = pred)
  
  pred_list_all[[m]] <- df_pred_m
  
  setTxtProgressBar(pb, m)
}

close(pb)

mean(fit_time)
mean(pred_time)

pred_list_all[[96]] %>% filter(t>0.2) %>% View()
pred_list_all[[99]] %>% filter(t>0.2) %>% View()

#### Save results ####

fit_time_subset_adglmm <- fit_time
pred_time_subset_adglmm <- pred_time
pred_subset_adglmm <- pred_list_all

save(fit_time_subset_adglmm, pred_time_subset_adglmm, pred_subset_adglmm, 
     file = here("Data/SubSimOutput_GLMMadaptive.RData"))






## separate by window
window <- seq(0, 1, by = 0.2)
M <- 1

#### ISE ####
# ise_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))
ise_mat <- df_pred_m %>%
  mutate(err1 = (pred0.2-eta_i)^2,
         err2 = (pred0.4-eta_i)^2,
         err3 = (pred0.6-eta_i)^2,
         err4 = (pred0.8-eta_i)^2) %>% 
  select(id, t, starts_with("err")) %>% 
  mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
  group_by(window, id) %>% 
  summarise_at(vars(err1, err2, err3, err4), sum) %>% 
  group_by(window) %>% 
  summarize_at(vars(err1, err2, err3, err4), mean) %>%
  filter(window != "[0,0.2]") %>% 
  select(starts_with("err"))
# ise_mat[, ,m] <- as.matrix(ise_tb)


# dims: prediction window, max obs time, simulation iter
mean_ise <- apply(ise_mat, c(1, 2), mean)
mean_ise <- data.frame(mean_ise) %>% 
  mutate(Window = c("(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "(0.8, 1.0]"),
         .before = 1)
colnames(mean_ise) <- c("Window", "0.2", "0.4", "0.6", "0.8")
mean_ise


#### Calculate AUC ####

## auc container 
# auc_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))

## a function to calculate AUC
get_auc <- function(y, pred){
  if(sum(is.na(y))>0 | sum(is.na(pred))>0){
    auc <- NA
  }
  else{
    this_perf <- performance(prediction(pred, y), measure = "auc")
    auc <- this_perf@y.values[[1]]
  }
  return(auc)
}


## calcualte AUC
auc_tb <- df_pred_m %>% 
  mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>%
  group_by(window) %>%
  summarise(auc1 = get_auc(Y, pred0.2),
            auc2 = get_auc(Y, pred0.4),
            auc3 = get_auc(Y, pred0.6),
            auc4 = get_auc(Y, pred0.8)) %>%
  filter(window != "[0,0.2]") %>% 
  select(starts_with("auc"))

#### save data ####
df_subset_adglmm <- df_pred_m
save(df_subset_adglmm, 
     file = here("Data/SubSimOutput_GLMMadaptive.RData"))


