# This script writes the fGFPCA algorithm into a function

#### set up ####

set.seed(1114)

# options(mc.cores=parallel::detectCores())

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
# library(LaplacesDemon)
library(rstan)
library(fastGFPCA)
theme_set(theme_minimal())



##### Function for prediction #####

# prediction for one subject


# data_i: include column id, index, mu, value
# efuncs
# evals
# mu
# tmax: end of observation track
# J: number of total observation

DynPred <- function(data_i, efuncs, evals, mu, tmax,
                    nchains = 2, nwarmup = 1000, niter = 2000, ncores = 1,
                    interval = T, int_qt = c(0.025, 0.075)){
  
  df_it <- data_i %>% filter(index <= tmax)
  J <- nrow(efuncs)
  K <- ncol(efuncs)

  # stan data    
  stanData <- list(
      J = J, # full measurement grid
      Ju = nrow(df_it), # number of observations
      Y = df_it$value, # outcome
      K = K, # number of PC functions
      efuncs = efuncs, # estimated PC functions
      b0 = mu, # estimated mean
      lambda = evals # estimated eigenvalues
  )
      
  fit <- stan(
    file = here("Code/prediction.stan"),  # Stan program
    data = stanData,    # named list of data
    chains = nchains,             # number of Markov chains
    warmup = nwarmup,          # number of warmup iterations per chain
    iter = niter,            # total number of iterations per chain
    cores = ncores,              # number of cores (could use one per chain)
    refresh = 0             # no progress shown
  )
  
  # point prediction
  
  scores_tmax <- summary(fit)$summary[1:K, "mean"]
  pred_eta <- mu+efuncs%*%scores_tmax
  
  if(interval == TRUE){
  # credible prediction interval using sampling quantiles
    score_draws <- as.matrix(fit)[, 1:K]
    eta_draws <- mu+efuncs %*% t(score_draws)
    eta_draws <- apply(eta_draws, 1, quantile, probs = c(0.025, 0.975))
    pred = data.frame(pred=pred_eta, 
                      pred_lb = eta_draws[1, ],
                      pred_ub = eta_draws[2, ])   
  }
  else{
    pred = data.frame(pred=pred_eta)
  }
    
  
  return(list(pred = pred_eta, scores = scores_tmax))

}

##### Test with toy data #####

# binomial data, with bins that do not overlap
df_gfpca <- sim_gfpca(N = 200, J = 200, case = 2)$df_gfpca
df_gfpca2 <- sim_gfpca(N = 200, J = 200, case = 2)$df_gfpca
df_test <- df_gfpca2 %>% filter(id==sample(1:200, 1) & index <= 0.75)


gfpca_mod <- fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 10, family = binomial)

gfpca_mod$efunctions
gfpca_mod$evalues

range(df_gfpca2$index)

test_out <- DynPred(data_i = df_test, efuncs = gfpca_mod$efunctions, 
        evals = gfpca_mod$evalues, mu = gfpca_mod$mu,
        tmax = max(df_test$index))

test_out$scores

df_gfpca2 %>% filter(id==unique(df_test$id)) %>% 
  add_column(test_out$pred) %>% 
  ggplot()+
  geom_line(aes(x=index, y=eta, col = "true"))+
  geom_line(aes(x=index, y = pred, col = "pred"))+
  geom_line(aes(x=index, y = pred_lb, col = "pred"), linetype = "dashed")+
  geom_line(aes(x=index, y = pred_ub, col = "pred"), linetype = "dashed")


##### test with simulated data #####

load(here("Data/sim_data.RData"))
df <- sim_data[[1]]
unique(df$id)

df_train <- df %>% filter(id %in% 1:500)
df_test <- df %>% filter(id %in% 501:600)

# fGFPCA
head(df_train)
head(df_gfpca)
df_train <- df_train %>% rename(index=t, value = Y) 
df_test <- df_test %>% rename(index=t, value = Y) 
gfpca_fit <- fast_gfpca(df_train, npc = 4, family = "binomial")

# estimates
gfpca_fit$efunctions
gfpca_fit$evalues
gfpca_fit$mu

mid <- sample(unique(df_test$id), 1)

test_out <- DynPred(data_i = df_test %>% 
                      filter(id==mid & index <= 0.35),
                    efuncs = gfpca_fit$efunctions, 
                    evals = gfpca_fit$evalues, mu = gfpca_fit$mu,
                    tmax = 0.55)

df_test %>% 
  filter(id==mid) %>% 
  add_column(test_out$pred) %>% 
  ggplot()+
  geom_line(aes(x=index, y=eta_i, col = "true"))+
  geom_line(aes(x=index, y = pred, col = "pred"))+
  geom_line(aes(x=index, y = pred_lb, col = "pred"), linetype = "dashed")+
  geom_line(aes(x=index, y = pred_ub, col = "pred"), linetype = "dashed")

#### Test on NHANES data ####
df_nhanes <- read_rds(here("Data/nhanes_bi.rds"))
head(df_nhanes)

df_nhanes <- df_nhanes %>% 
  rename(id = SEQN, value = Z, index = sind)

N <- length(unique(df_nhanes$id))
train_id <- sample(unique(df_nhanes$id), size = 500)

df_train <- df_nhanes %>% filter(id %in% train_id)
df_test <- df_nhanes %>% filter(!id %in% train_id)

# fGFPCA
head(df_train)
gfpca_fit <- fast_gfpca(df_train, pve = 0.9, family = "binomial")

# estimates
dim(gfpca_fit$efunctions)
gfpca_fit$evalues
gfpca_fit$mu

mid <- sample(unique(df_test$id), 1)

test_out <- DynPred(data_i = df_test %>% 
                      filter(id==mid & index <= 1000),
                    efuncs = gfpca_fit$efunctions, 
                    evals = gfpca_fit$evalues, mu = gfpca_fit$mu,
                    tmax = 1000)

df_test %>% 
  filter(id==mid) %>% 
  add_column(test_out$pred) %>% 
  mutate_at(vars(starts_with("pred")), function(x){exp(x)/(1+exp(x))}) %>%
  ggplot()+
  geom_point(aes(x=index, y=value), alpha = 0.2)+
  geom_line(aes(x=index, y = pred, col = "pred"))+
  geom_line(aes(x=index, y = pred_lb, col = "pred"), linetype = "dashed")+
  geom_line(aes(x=index, y = pred_ub, col = "pred"), linetype = "dashed")


