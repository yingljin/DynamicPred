# This script implements the full simulation 
# for dynamic prediction using GLMMadaptive
# as a completing method 
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

# load simulated data set. From Code/GenData.R
load(here("Data/sim_data.RData"))

#### Data overview ####

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


#### GLMMadaptive model estimation ####

# split data into train and test set
N_train <- 500
N_test <- 100

# containers
pred_list_ref <- list()
fit_time_ref <- pred_time_ref <- rep(NA, M)


pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

for(m in 1:M){
  
  # data plit
  df <- sim_data[[m]]
  train_df <- df %>% filter(id %in% 1:N_train)
  test_df <- df %>% filter(!id %in% 1:N_train)
  
  t1 <- Sys.time()
  # fit GLMM adaptive model: linear time fixed and random effect
  adglmm <- mixed_model(Y ~ t, random = ~ t | id,
                        data =  train_df, family = binomial())
  t2 <- Sys.time()
  fit_time_ref[m] <- t2-t1 # minutes
  
  # prediction on test set
  t1 <- Sys.time()
  pred_m <- list(
    predict(adglmm, newdata = test_df %>% filter(t <= 0.2),
            newdata2 =  test_df %>% filter(t>0.2), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(adglmm, newdata = test_df %>% filter(t <= 0.4),
            newdata2 =  test_df %>% filter(t>0.4), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(adglmm, newdata = test_df %>% filter(t <= 0.6),
            newdata2 =  test_df %>% filter(t>0.6), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2,
    
    predict(adglmm, newdata = test_df %>% filter(t <= 0.8),
            newdata2 =  test_df %>% filter(t>0.8), 
            type = "subject_specific", type_pred = "link",
            se.fit = FALSE, return_newdata = TRUE)$newdata2
  )
  t2 <- Sys.time()
  pred_time_ref[m] <- t2-t1 # seconds
  
  # format
  pred_m <- lapply(pred_m, function(x){x %>% select(id, t, pred)})
  pred_list_ref[[m]] <- test_df %>% 
    left_join(pred_m[[1]]) %>% 
    rename(pred0.2 = pred) %>% 
    left_join(pred_m[[2]]) %>% 
    rename(pred0.4 = pred) %>% 
    left_join(pred_m[[3]]) %>% 
    rename(pred0.6 = pred) %>% 
    left_join(pred_m[[4]]) %>% 
    rename(pred0.8 = pred)
  
  
  setTxtProgressBar(pb, m)
}



save(pred_list_ref, fit_time_ref, pred_time_ref, 
     file = here("Data/SimOutput_GLMMadaptive.RData"))


#### Calculate ISE ####

# load simulation results
load(here("Data/SimOutput_GLMMadaptive.RData"))
mean(fit_time_ref)
mean(pred_time_ref)

head(pred_list_ref[[1]])

## prediction window 
window <- seq(0, 1, by = 0.2)
M <- length(pred_list_ref)

## ISE container 
ise_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))
# dims: prediction window, max obs time, simulation iter

## calculation
for(m in 1:M){
  this_df <- pred_list_ref[[m]]
  ise_tb <- pred_list_ref[[m]] %>%
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
  ise_mat[, ,m] <- as.matrix(ise_tb)
  
  
}

mean_ise <- apply(ise_mat, c(1, 2), mean)
mean_ise <- data.frame(mean_ise) %>% 
  mutate(Window = c("(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "(0.8, 1.0]"),
         .before = 1)
colnames(mean_ise) <- c("Window", "0.2", "0.4", "0.6", "0.8")
mean_ise


#### Calculate AUC ####

## auc container 
auc_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))

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

get_auc(pred_list_ref[[1]]$Y[pred_list_ref[[1]]$t>0.2], 
        pred_list_ref[[1]]$pred0.2[pred_list_ref[[1]]$t>0.2])

## calcualte AUC

for(m in 1:M){
  this_df <- pred_list_ref[[m]]
  auc_tb <- this_df %>% 
    mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
    select(Y, starts_with("pred"), window) %>%
    group_by(window) %>%
    summarise(auc1 = get_auc(Y, pred0.2),
              auc2 = get_auc(Y, pred0.4),
              auc3 = get_auc(Y, pred0.6),
              auc4 = get_auc(Y, pred0.8)) %>%
    filter(window != "[0,0.2]") %>% 
    select(starts_with("auc"))
  auc_mat[, ,m] <- as.matrix(auc_tb)
  
}


mean_auc <- apply(auc_mat, c(1, 2), mean)
mean_auc <- data.frame(mean_auc) %>% 
  mutate(Window = c("(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "(0.8, 1.0]"),
         .before = 1)
colnames(mean_auc) <- c("Window", "0.2", "0.4", "0.6", "0.8")
mean_auc
