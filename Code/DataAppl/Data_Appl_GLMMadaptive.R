# data application of NHANES
# using GLMMadaptive
# corresponding to manuscript Section 5

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
df <- read_rds(here("DataRaw/nhanes_bi_sub.rds"))
df <- df %>% rename(id=SEQN, Y=Z)

N <- length(unique(df$id)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject
t <- unique(df$sind)
K <- 4 # number of eigenfunctions to use


#### Data split ####

# 60% subjects for training, 40% for out-of-sample prediction
train_id <- sample(unique(df$id), size = N*0.6)
test_id <- setdiff(unique(df$id), train_id)

train_df <- df %>% filter(id %in% train_id)
test_df <- df %>% filter(id %in% test_id)


#### fit GLMMadaptvie ####
# on training set
# only able to fit a random intecept
head(train_df)

t1 <- Sys.time()
adglmm_mod <- mixed_model(Y ~ t, random = ~ 1 | id, 
                          data = train_df %>% mutate(t=sind/J),
                          family = binomial())
t2 <- Sys.time()
t_est_adglmm <- t2-t1 # model fitting took 20.82 mins
summary(adglmm_mod)

# prediction
head(test_df)
test_df <- test_df %>% mutate(t=sind/J)


df_pred1 <- predict(adglmm_mod, 
                    newdata = test_df %>% filter(sind <= 360),
                    newdata2 =  test_df %>% filter(sind > 360), 
                    type = "subject_specific", type_pred = "link",
                    se.fit = FALSE, return_newdata = TRUE)$newdata2
head(df_pred1)
df_pred1 <- df_pred1 %>% rename(pred360=pred)

df_pred2 <- predict(adglmm_mod, 
                    newdata = test_df %>% filter(sind <= 720),
                    newdata2 =  test_df %>% filter(sind > 720), 
                    type = "subject_specific", type_pred = "link",
                    se.fit = FALSE, return_newdata = TRUE)$newdata2
df_pred2 <- df_pred2 %>% rename(pred720=pred)
head(df_pred2)

df_pred3 <- predict(adglmm_mod, 
                    newdata = test_df %>% filter(sind <= 1080),
                    newdata2 =  test_df %>% filter(sind > 1080), 
                    type = "subject_specific", type_pred = "link",
                    se.fit = FALSE, return_newdata = TRUE)$newdata2
df_pred3 <- df_pred3 %>% rename(pred1080=pred)
head(df_pred3)


#### check results ###
pred_nhanes_adglmm <- test_df %>% 
  left_join(df_pred1 %>% select(id, sind, pred360)) %>% 
  left_join(df_pred2 %>% select(id, sind, pred720)) %>% 
  left_join(df_pred3 %>% select(id, sind, pred1080)) 
  

pred_nhanes_adglmm %>% 
  filter(id %in% sample(test_id, 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred360, col = "6am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred720, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred1080, col = "6pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

save(pred_nhanes_adglmm,
     file = here("Data/ApplOutput_GLMMadaptive.RData"))
