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



#### fGFPCA ####

# functions
K <- 4
source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))


# binning 
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(min(df_train_m$sind)-1, max(df_train_m$sind), by = bin_w) # cutoff points on time index
mid <- (brks+bin_w/2)[1:n_bin] # interval mid point time index
t <- unique(df_train_m$t)
mid_t <- unique(df_train_m$t[df_train_m$sind %in% mid]) # actually time corresponding to bin mid point

# values used for projection
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)


# result container
# M <- 10
M <- 1
# pred_list_all <- list()
# converge_state_list <- list()
# fit_time <- pred_time <- rep(NA, M)

#### Model fit ####

## Step 1: bin observations
df_train_m$bin <- cut(df_train_m$sind, breaks = brks, 
                      include.lowest = T, labels = mid)
df_train_m$bin <- as.numeric(as.character(df_train_m$bin))


## fit local GLMM model on the training set
t1 <- Sys.time()
## Step 2: local GLMM and estimate latent function
## on training set
bin_lst <- split(df_train_m, f = df_train_m$bin)
df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) # singularity issue
df_est_latent <- bind_rows(df_est_latent) 
range(df_est_latent$eta_hat)

## Step 3: FPCA
uni_eta_hat <- df_est_latent %>% filter(bin==sind)
mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                         nrow=Ntr, 
                         ncol=n_bin, byrow = F) 
fpca_mod <- fpca.face(mat_est_unique, pve=0.999, argvals = mid_t, var=T)
K <- 4
fpca_mod$efunctions
fpca_mod$evalues

## Step 4: project and debias
## Projection
B <- spline.des(knots = knots_values, x = mid_t, ord = p + 1,
                outer.ok = TRUE)$design  # evaluate B-splines on binned grid
Bnew <- spline.des(knots = knots_values, x = t, ord = p + 1,
                   outer.ok = TRUE)$design  # evaluate B-splines on original grid
df_phi <- matrix(NA, J, K) 
for(k in 1:K){
  lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
  df_phi[,k] <- Bnew %*% coef(lm_mod)
}# project binned eigenfunctions onto the original grid
## debias
df_phi <- data.frame(t = t, df_phi)
colnames(df_phi) <- c("t", paste0("phi", 1:4))
df_train_m <- df_train_m %>% left_join(df_phi, by = "t")
df_train_m$id <- droplevels(df_train_m$id)
debias_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                     s(id, by=phi1, bs="re")+
                     s(id, by=phi2, bs="re")+
                     s(id, by=phi3, bs="re")+
                     s(id, by=phi4, bs="re"), 
                   family = binomial, data=df_train_m, 
                   method = "fREML",
                   discrete = TRUE)
new_mu <- predict(debias_glmm, type = "terms")[1:J, 1] # extract re-evaluated mean
new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
# rescale
new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
range(new_phi)
new_phi <- as.matrix(new_phi)

new_lambda <- new_lambda/n_bin
t2 <- Sys.time()
t2-t1 # less than 1 second


#### Prediction ####
# Prediction pf test sample using Laplace approximation
pred_list_m <- list()
converge_state_m <- matrix(NA, Nte, 4)
t1 <- Sys.time()
# per subject
df_test_m$id <- droplevels(df_test_m$id)
test_id <- unique(df_test_m$id)

for(i in seq_along(test_id)){
  df_i <- df_test_m %>% filter(id==test_id[i])
  
  # per max obs time
  for(tmax in c(0.35, 0.45, 0.55, 0.65)){
    df_it <- df_i %>% filter(t <= tmax)
    max_rid <- nrow(df_it)
    
    # into a list
    MyData <- list(K=K, 
                   X=new_phi[1:max_rid, ], 
                   mon.names=mon.names,
                   parm.names=parm.names, 
                   pos.xi=pos.xi, 
                   y=df_it$Y, 
                   tao=diag(new_lambda), f0=new_mu[1:max_rid])
    
    
    # fit laplace approximation
    Fit <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData, Method = "NM", Iterations = 1000,
                                CovEst = "Identity")
    converge_state_m[i, which(c(0.35, 0.45, 0.55, 0.65)==tmax)] <- Fit$Converged
    score <- Fit$Summary1[, "Mode"]
    
    # prediction
    eta_pred_out <- new_mu+new_phi%*%score
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
  }
  
  pred_list_m[[i]] <- df_i
}
t2 <- Sys.time()
t2-t1 # About 8 seconds
# pred_time[m] <- t2-t1

# save predictions
df_pred <- bind_rows(pred_list_m)
df_pred$pred0.35[df_pred$t<=0.35] <- NA
df_pred$pred0.45[df_pred$t<=0.45] <- NA
df_pred$pred0.55[df_pred$t<=0.55] <- NA
df_pred$pred0.65[df_pred$t<=0.65] <- NA

mean(converge_state_m) # all converged

# pred_list_all[[m]] <- df_pred
df_pred_m_fGFPCA <- df_pred
save(df_pred_m_fGFPCA, 
     file = here("Data/SubSimOutput_fGFPCA.RData"))




# save convergence status
converge_state_list[[m]] <- converge_state_m


