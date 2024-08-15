
# In this script I will try to investigate the prediction interval of GFOSR model

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

load(here("SimOneIter.RData"))


N <- length(unique(df_m$id))
J <- length(unique(df_m$t))
t <- unique(df_m$t)
K <- 4 # number of eigenfunctions to use

# overview of a few samples
df_m %>% filter(id %in% sample(N, size = 4)) %>% 
  ggplot()+
  geom_point(aes(x=t, y=Y))+
  geom_line(aes(x=t, y=exp(eta_i)/(1+exp(eta_i))), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)


N_train <- 500 # training sample size
N_test <- 100 # testing sample szie 
# L <- 5 # number of last observations used as time-fixed predictor
L <- 1
windows <- seq(0, 1, by = 0.2) # prediction window


# prediction window 
df_m$window <- cut(df_m$t, breaks=windows, labels = 1:5, include.lowest = T)

# split data
train_df <- df_m %>% filter(id %in% 1:N_train)
test_df <- df_m %>% filter(!id %in% 1:N_train)
train_df$window <- as.numeric(as.character(train_df$window))
test_df$window <- as.numeric(as.character(test_df$window))

# let's use the 'observation up to 0.2' case as an example

# first, format training data
## Predictor
y_obs_max <- train_df %>% filter(window==1) %>% 
  group_by(id) %>%
  slice_max(t, n=L) %>% # only keep the last L observations of each subject
  select(id, Y, sind) %>%
  mutate(name = sind-min(sind)+1) %>%
  pivot_wider(id_cols = id, values_from = Y, names_from = "name", names_prefix = "yl")
  
## outcome (training set)
df_pred_tr <- train_df %>% filter(window >1) %>%
  left_join(y_obs_max, by = "id")

## fit model
fit_gen_fosr <- bam(Y ~ s(t, bs="cr", k=20) + 
                      s(t, bs="cr", k=20, by = yl1),
                    family = binomial, data=df_pred_tr, 
                    method = "fREML",
                    discrete = TRUE)

summary(fit_gen_fosr)
  
# now we can do out-of-sample prediction
# for simplicity, let's use only one subject
## predictor (testing set)
last_obs <- test_df %>% filter(id == 501 & window==1) %>%
  slice_max(t, n=L) %>% select(Y)
# y_obs_max_te <- test_df %>% filter(id == 501 & window==1) %>%
#   slice_max(t, n=L) %>% select(id, Y, sind) %>% 
#   mutate(name = sind-min(sind)+1) %>%
#   pivot_wider(id_cols = id, values_from = Y, names_from = name, names_prefix = "yl")

## outcome (testing set)
df_pred_te <- test_df %>% filter(id == 501 & (window > 1)) %>%
  # left_join(y_obs_max_te, by = "id") %>% 
  # mutate_at(vars(Y, starts_with("yl")), as.factor)
  mutate(yl1 = unlist(last_obs))

## predict using the FOSR model fit above
## include point prediction and standard deviation
pred_mw <- predict(fit_gen_fosr, newdata = df_pred_te, 
                   type = "link", # on the linear predictor scale
                   se.fit = T)


# plot the standard deviation of latent track
# a very strange U shape
plot(pred_mw$se.fit)
  
#### Manual calculation ####
# here I try to calculate standard deviation of latent track 
# I can extract the variance of coefficients from the GFOSR model
coef_var <- diag(vcov(fit_gen_fosr))

# since the last observation of subject 501 is 0, we only need the 
# functional intercept part
coef_var[1:20]

# and  I need to know the spline basis functions
## functional intercept
Xmat <- model.matrix.gam(fit_gen_fosr)
Xmat <- Xmat[!duplicated(Xmat[, 2]), ]

# calculate standard deviation
sd_eta <- sqrt(Xmat[,1:20]^2 %*% coef_var[1:20])
plot(sd_eta)
# this hand-calculated SD still looks strange

# compare the two SD estiamtes
plot(sd_eta)
points(pred_mw$se.fit, col = "red")
