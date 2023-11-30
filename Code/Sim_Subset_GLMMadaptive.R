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




#### load simulated data ####

load(here("Data/sim_data.RData"))

# reduce data
## training data
sub_sim_data <- sim_data[1:100] 
M <- 100 # number of simulations

sub_sim_data <- lapply(sub_sim_data, function(x){x %>% filter(id %in% c(1:10,501:510))})
Ntr <- Nte <- 10 # trainig and testing sample size

sub_sim_data <- lapply(sub_sim_data, function(x){x %>% filter(t > 0.25 & t <= 0.75)})
J <- length(unique(sub_sim_data[[1]]$t))

lapply(sub_sim_data, dim)
head(sub_sim_data[[1]])


#### One iteration example ####

df_train_m <- sub_sim_data[[1]] %>% filter(id %in% 1:10)
df_test_m <- sub_sim_data[[1]] %>% filter(!id %in% 1:10)

#### GLMM adaptive ####
## model fit
t1 <- Sys.time()
fit_adglmm <- mixed_model(fixed = Y ~ 0+ns(t, df = 5), 
                          random = ~ 0+ns(t, df = 5) | id, 
                          data =  df_train_m, family = binomial())
t2 <- Sys.time()
t_sub_fit <- t2-t1 ## 22.5 minutes

## prediction
t1 <- Sys.time()
pred_m <- list(
  predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.35),
          newdata2 =  df_test_m %>% filter(t>0.35), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.45),
          newdata2 =  df_test_m %>% filter(t>0.45), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.55),
          newdata2 =  df_test_m %>% filter(t>0.55), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.65),
          newdata2 =  df_test_m %>% filter(t>0.65), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2)
t2 <- Sys.time()
t2-t1 # 0.16 seconds

pred_m <- lapply(pred_m, function(x){x %>% select(id, t, pred)})
head(pred_m[[1]])
head(df_test_m)
df_pred_m <- df_test_m %>% 
  left_join(pred_m[[1]]) %>% 
  rename(pred0.35 = pred) %>% 
  left_join(pred_m[[2]]) %>% 
  rename(pred0.45 = pred) %>% 
  left_join(pred_m[[3]]) %>% 
  rename(pred0.55 = pred) %>% 
  left_join(pred_m[[4]]) %>% 
  rename(pred0.65 = pred)

## separate by window
window <- seq(0.25, 0.75, by = 0.1)
M <- 1

#### ISE ####
# ise_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))
ise_mat <- df_pred_m %>%
  mutate(err1 = (pred0.35-eta_i)^2,
         err2 = (pred0.45-eta_i)^2,
         err3 = (pred0.55-eta_i)^2,
         err4 = (pred0.65-eta_i)^2) %>% 
  select(id, t, starts_with("err")) %>% 
  mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
  group_by(window, id) %>% 
  summarise_at(vars(err1, err2, err3, err4), sum) %>% 
  group_by(window) %>% 
  summarize_at(vars(err1, err2, err3, err4), mean) %>%
  filter(window != "[0.25,0.35]") %>% 
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
  summarise(auc1 = get_auc(Y, pred0.35),
            auc2 = get_auc(Y, pred0.45),
            auc3 = get_auc(Y, pred0.55),
            auc4 = get_auc(Y, pred0.65)) %>%
  filter(window != "[0.25,0.35]") %>% 
  select(starts_with("auc"))

#### save data ####
df_pred_m_adglmm <- df_pred_m
save(df_pred_m_adglmm, 
     file = here("Data/SubSimOutput_GLMMadaptive.RData"))


