# data application of NHANES

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
library(LaplacesDemon)
library(splines)

#### load data ####
df <- read_rds(here("Data/nhanes_bi.rds"))
df <- df %>% rename(id=SEQN, Y=Z)
head(df)

N <- length(unique(df$id)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject
t <- unique(df$sind)

#### Data split ####
# 60% (5257) subjects for training, 40% (3506) for out-of-sample prediction
train_id <- sample(unique(df$id), size = N*0.6)
test_id <- setdiff(unique(df$id), train_id)

train_df <- df %>% filter(id %in% train_id)
test_df <- df %>% filter(id %in% test_id)

Ntr <- length(train_id)
Nte <- length(test_id)

#### GFOSR model ####

window <- seq(0, J, by = 360)
train_df$window <- cut(train_df$sind, breaks=window, 
                       labels = 1:4, include.lowest = T)
train_df$window <- as.numeric(as.character(train_df$window))
test_df$window <- cut(test_df$sind, breaks=window, 
                       labels = 1:4, include.lowest = T)
test_df$window <- as.numeric(as.character(test_df$window))

# L <- 1
L <- 5

# model fit
t1 <- Sys.time()
for(w in 1:3){
  # model fit
  # predictor (training set)
  y_obs_max <- train_df %>% filter(window==w) %>% group_by(id) %>%
    slice_max(sind, n=L) %>% select(id, Y, sind) %>%
    mutate(name = sind-min(sind)+1) %>%
    pivot_wider(id_cols = id, values_from = Y, names_from = "name", names_prefix = "yl")
  # outcome (training set)
  df_pred_tr <- train_df %>% filter(window > w) %>%
    left_join(y_obs_max, by = "id") %>% 
    mutate_at(vars(Y, starts_with("yl")), as.factor) 
  fit_gen_fosr <- bam(Y ~ s(sind, bs="cr", k=20) +
                        s(sind, bs="cr", k=20, by = yl5)+
                        s(sind, bs="cr", k=20, by = yl4)+
                        s(sind, bs="cr", k=20, by = yl3)+
                        s(sind, bs="cr", k=20, by = yl2)+
                        s(sind, bs="cr", k=20, by = yl1),
                      family = binomial, data=df_pred_tr,
                      method = "fREML",
                      discrete = TRUE)
  # fit_gen_fosr <- bam(Y ~ s(sind, bs="cr", k=20) +
  #                       s(sind, bs="cr", k=20, by = yl1),
  #                     family = binomial, data=df_pred_tr,
  #                     method = "fREML",
  #                     discrete = TRUE)
  
  # prediction
  # predictor (testing set)
  y_obs_max_te <- test_df %>% filter(window==w) %>% group_by(id) %>%
    slice_max(sind, n=L) %>% select(id, Y, sind) %>% 
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
  

# pred_appl_gfosr_l1 <- test_df
# time_appl_gfosr_l1 <- t2-t1
pred_appl_gfosr_l5 <- test_df
time_appl_gfosr_l5 <- t2-t1

# save(pred_appl_gfosr_l1, time_appl_gfosr_l1, file = here("Data/Appl_GFOSR_L1.RData"))
save(pred_appl_gfosr_l5, time_appl_gfosr_l5, file = here("Data/Appl_GFOSR_L5.RData"))
  

  
test_df %>% filter(sind>350)


test_df %>% 
  filter(id %in% sample(test_id, 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred_w1, col = "6am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_w2, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_w3, col = "6pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

nhanes_pred_adglmm <- df_pred



save(nhanes_pred_fgfpca, nhanes_pred_adglmm,
     file = here("Data/ApplOutput.RData"))
