# This script implements the full simulation 
# for dynamic prediction of fGFPCA
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

#### fGFPCA Model estimation ####

# functions
source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))

# split data into train and test set
N_train <- 500
N_test <- 100

# binning 
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points on time index
mid <- (brks+bin_w/2)[1:n_bin] # interval mid point time index
mid_t <- t[mid] # actually time corresponding to bin mid point

# values used for projection
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)


# result container
# M <- 10
M
pred_list_all <- list()
converge_state_list <- list()
fit_time <- pred_time <- rep(NA, M)



pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

# The whole process
for(m in 1:M){
  
  # for every simulated dataset
  df <- sim_data[[m]]
  
  ## Step 1: bin observations
  df$bin <- cut(df$sind, breaks = brks, include.lowest = T, labels = mid)
  df$bin <- as.numeric(as.character(df$bin))
  
  # split into training and test set
  train_df <- df %>% filter(id %in% 1:N_train)
  test_df <- df %>% filter(!id %in% 1:N_train)
  
  # fit fGFPCA model on the training set
  t1 <- Sys.time()
  ## Step 2: local GLMM and estimate latent function
  ## on training set
  bin_lst <- split(train_df, f = train_df$bin)
  df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) 
  df_est_latent <- bind_rows(df_est_latent) 
  ## Step 3: FPCA
  uni_eta_hat <- df_est_latent %>% filter(bin==sind)
  mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                           nrow=N_train, 
                           ncol=n_bin, byrow = F) 
  fpca_mod <- fpca.face(mat_est_unique, argvals = mid_t, var=T)
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
  train_df <- train_df %>% left_join(df_phi, by = "t")
  train_df$id <- as.factor(train_df$id)
  debias_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                       s(id, by=phi1, bs="re")+
                       s(id, by=phi2, bs="re")+
                       s(id, by=phi3, bs="re")+
                       s(id, by=phi4, bs="re"), 
                     family = binomial, data=train_df, 
                     method = "fREML",
                     discrete = TRUE)
  new_mu <- predict(debias_glmm, type = "terms")[1:J, 1] # extract re-evaluated mean
  new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
  # rescale
  new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
  new_phi <- as.matrix(new_phi)
  new_lambda <- new_lambda/n_bin
  t2 <- Sys.time()
  fit_time[m] <- t2-t1
  
  # Prediction pf test sample using Laplace approximation
  # container
  pred_list_m <- list()
  converge_state_m <- matrix(NA, N_test, 4)
  t1 <- Sys.time()
  # per subject
  test_id <- unique(test_df$id)
  for(i in seq_along(test_id)){
    df_i <- test_df %>% filter(id==test_id[i])
    
    # per max obs time
    for(tmax in c(0.2, 0.4, 0.6, 0.8)){
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
      converge_state_m[i, which(c(0.2, 0.4, 0.6, 0.8)==tmax)] <- Fit$Converged
      score <- Fit$Summary1[, "Mode"]
      
      # prediction
      eta_pred_out <- new_mu+new_phi%*%score
      df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
    }
    
    pred_list_m[[i]] <- df_i
  }
  t2 <- Sys.time()
  t_pred <- t2-t1 # About 3.5 minutes
  pred_time[m] <- t2-t1
  
  # save predictions
  df_pred <- bind_rows(pred_list_m)
  df_pred$pred0.2[df_pred$t<=0.2] <- NA
  df_pred$pred0.4[df_pred$t<=0.4] <- NA
  df_pred$pred0.6[df_pred$t<=0.6] <- NA
  df_pred$pred0.8[df_pred$t<=0.8] <- NA
 
  pred_list_all[[m]] <- df_pred
  
  # save convergence status
  converge_state_list[[m]] <- converge_state_m
  
  
  setTxtProgressBar(pb, m)
}

close(pb)




#### Check output ####

mean(fit_time) # average time for model fitting: 26 secs
mean(pred_time) # average time for prediction: 133 secs 

pred_list_all %>% lapply(dim)
length(pred_list_all)

# convergence status
mean(sapply(converge_state_list, mean)) # all dataset converged

# visualize test samples
rand_id <- sample(test_id, 4)

# prediction results
pred_list_all[[476]] %>% 
  filter(id %in% rand_id) %>% 
  mutate_at(vars(eta_i, pred0.2, pred0.4, pred0.6, pred0.8), function(x){exp(x)/(1+exp(x))}) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.2)+
  geom_line(aes(x=sind, y=eta_i, col = "True"))+
  geom_line(aes(x=sind, y=pred0.2, col = "0.2"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred0.4, col = "0.4"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred0.6, col = "0.6"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred0.8, col = "0.8"), linetype="dashed", na.rm=T)+
  facet_wrap(~id)




#### Save results ####

save(pred_list_all, fit_time, pred_time, converge_state_list,
     file = here("Data/SimOutput_fGFPCA.RData"))

#### Calculate ISE ####

# load simulation results
load(here("Data/SimOutput_fGFPCA.RData"))

dim(pred_list_all[[1]])

## prediction window 
window <- seq(0, 1, by = 0.2)
M <- length(pred_list_all)

## ISE container 
ise_mat <- array(NA, dim = c(length(window)-2, length(window)-2, M))
# dims: prediction window, max obs time, simulation iter

## calculation
for(m in 1:M){
  this_df <- pred_list_all[[m]]
  ise_tb <- pred_list_all[[m]] %>%
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
mean_ise <- data.frame(mean_ies) %>% 
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

get_auc(pred_list_all[[1]]$Y[pred_list_all[[1]]$t>0.2], 
        pred_list_all[[1]]$pred0.2[pred_list_all[[1]]$t>0.2])

## calcualte AUC

for(m in 1:M){
  this_df <- pred_list_all[[m]]
  auc_tb <- df_pred %>% 
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

#### Convergence ####
mean(sapply(converge_state_list, mean)) # 100% precent convergence 
