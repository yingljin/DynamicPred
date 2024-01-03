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


#### Model ####
# g(E[Y_i(s)|Y_i(t), t<s]) = \beta_0(s) + \beta_1 (s)* Y_i(t) for t < s
# This is essential a function-on-function regression
# given history as a function covariate
# to predict future as a functional outcome

# Let's start with a one-simulation example
df_m <- sim_data[[1]]
N_train <- 500
N_test <- 100

# split data
train_df <- df_m %>% filter(id %in% 1:N_train)
test_df <- df_m %>% filter(!id %in% 1:N_train)

# fit model on the training set
# let's say we wanna fit the concurrent functional model 
# the outcome must be evaluated at the same time as the covariate
# which means prediction window and sobservation window must be equal-length
# therefore, for each dataset, we would need to fit 10 models:
# observed 0-0.2, predict 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0
# remaining scenarios not clear. Let's wait for now. 


# let me try a FOF regression first
## time window
windows <- seq(0, 1, by = 0.2)
train_df$window <- cut(train_df$t, breaks=windows, labels = 1:5,
                       include.lowest = T)
class(train_df$window)
table(train_df$window, useNA = "always") # each window has equal number of observations

# I would like to try a gam model for generalized FOF
# first some data formatting
## observed 0-0.2, predict 0.2-0.4
dfm1 <- train_df %>% filter(window==2) %>% 
  rename(Yf=Y) %>%
  left_join(train_df %>% filter(window==1) %>% select(id, sind, Y) %>% 
              mutate(sind=sind+200) %>%
              rename(Yh=Y))

## gam model
fitm1 <- bam(Yf ~ s(t, bs="cc", k=10) + s(t, bs="cc", k=10, by = Yh) , 
                   family = binomial, data=dfm1, 
                   method = "fREML",
                   discrete = TRUE)

## prediction?
test_df$window <- cut(test_df$t, breaks=windows, labels = 1:5,
                       include.lowest = T)
table(test_df$window, useNA = "always")
dfm1_te <- test_df %>% filter(window==2) %>% 
  rename(Yf=Y) %>%
  left_join(test_df %>% filter(window==1) %>% select(id, sind, Y) %>% 
              mutate(sind=sind+200) %>%
              rename(Yh=Y))
head(dfm1_te)

pred_m1 <- predict(fitm1, newdata = dfm1_te)
length(pred_m1)
dfm1_te$pred_eta <- pred_m1

dfm1_te %>% 
  filter(id %in% 501:504) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta_i, col = "true"))+
  geom_line(aes(x=t, y=pred_eta, col = "pred"))+
  facet_wrap(~id)
