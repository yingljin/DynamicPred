# This script implements the full simulation for dynamic prediction of fGFPCA
# correponding to Section 4.4.1 in the manuscript

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
library(rstan)
theme_set(theme_minimal())


#### load simulated data ####

load(here("Data/sim_data.RData")) # data generated from Code/GenData.R

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
  geom_line(aes(x=t, y=exp(eta_i)/(1+exp(eta_i))), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)

#### fGFPCA Model estimation ####

# functions
source(here("Code/Functions/GLMM-FPCA.R")) 
# source(here("Code/OutSampBayes.R"))

# split data into train and test set
N_train <- 500
N_test <- 100

# binning 
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points on time index
mid <- (brks+bin_w/2)[1:n_bin] # interval mid point time index
mid_t <- t[mid] # actually time corresponding to bin mid point

# values used for projection in step 4
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)


# result container
M <- 3 # try on a few datasets first
# M
pred_list_fGFPCA <- list()
# score_list_all <- list()
# sd_list_all <- list()
# converge_state_list <- list()
time_fGFPCA <- rep(NA, M)


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
  new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1] # extract re-evaluated mean
  ## extract variacne from the debiased model
  new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
  
  # rescale
  new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
  new_phi <- as.matrix(new_phi)
  new_lambda <- new_lambda/n_bin
  # t2 <- Sys.time()
  # fit_time[m] <- difftime(t2, t1, units = "mins")
  
  # Prediction pf test sample using stan
  # container
  pred_list_m <- list()
  # converge_state_m <- matrix(NA, N_test, 4)
  # score_list_m <- array(NA, dim = c(N_test, K, 4))
  # sd_list_m <- array(NA, dim = c(N_test, K, 4)) # subject by PC by case
  test_id <- unique(test_df$id)
  tmax_vec <- c(0.2, 0.4, 0.6, 0.8)
  
  # t1 <- Sys.time()
  for(i in seq_along(test_id)){
    
    df_i <- test_df %>% filter(id==test_id[i])
    
    # container for one subject
    # scores <- matrix(NA, K, length(tmax_vec))
    # sd_scores <- matrix(NA, K, length(tmax_vec))
    
    for(tid in seq_along(tmax_vec)){
      
      tmax <- tmax_vec[tid]
      df_it <- df_i %>% filter(t <= tmax)
      
      # into a list
      stanData <- list(
        J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
        efuncs = new_phi, 
        b0 = new_mu,
        lambda = new_lambda
      )
      
      fit <- stan(
        file = here("Code/Functions/prediction.stan"),  # Stan program
        data = stanData,    # named list of data
        chains = 2,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0             # no progress shown
      )
      
      # point prediction
      scores_tmax <- summary(fit)$summary[1:K, "mean"]
      eta_pred_out <- new_mu+new_phi%*%scores_tmax
      df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
      
      # prediction interval using sampling quantiles
      score_draws <- as.matrix(fit)[, 1:4]
      eta_draws <- new_phi %*% t(score_draws)
      eta_draws <- apply(eta_draws, 1, quantile, probs = c(0.025, 0.975))
      # quantile interval
      df_i[ , paste0("pred", tmax, "_lb")] <- as.vector(new_mu+eta_draws[1, ])
      df_i[ , paste0("pred", tmax, "_ub")] <- as.vector(new_mu+eta_draws[2, ])
      
      
    }
    
   pred_list_m[[i]] <- df_i
    
  }
  t2 <- Sys.time()
  # t_pred <- t2-t1 # About 3.5 minutes
  time_fGFPCA[m] <- difftime(t2, t1, units = "mins")
  
  # save results
  pred_list_fGFPCA[[m]]<- bind_rows(pred_list_m)
  # converge_state_list[[m]] <- converge_state_m
  
  setTxtProgressBar(pb, m)
}

close(pb)


#### Check output ####

# time 
mean(time_fGFPCA)

# prediction
# figure
rand_id <- sample(test_id, 4)

df_exp <- pred_list_fGFPCA[[1]] %>%
  filter(id %in% rand_id) %>%
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  
df_exp[df_exp$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub",
                        "pred0.2_qlb", "pred0.2_qub")] <- NA
df_exp[df_exp$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub",
                        "pred0.4_qlb", "pred0.4_qub")] <- NA
df_exp[df_exp$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub",
                        "pred0.6_qlb", "pred0.6_qub")] <- NA
df_exp[df_exp$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub",
                        "pred0.8_qlb", "pred0.8_qub")] <- NA

df_exp %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.2_lb, ymax = pred0.2_ub,
                  col = "0.2", fill = "0.2", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.4_lb, ymax = pred0.4_ub,
                  col = "0.4", fill = "0.4", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.6_lb, ymax = pred0.6_ub,
                  col = "0.6", fill = "0.6", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.8_lb, ymax = pred0.8_ub,
                  col = "0.8", fill = "0.8", alpha = 0.1),
              linetype="dashed")+
  facet_grid(rows = vars(id))+
  guides(alpha = "none", col="none")
  

#### Save results ####

save(pred_list_fGFPCA, time_fGFPCA,
     file = here("Data/SimOutput_fGFPCA.RData"))

# save(pred_list_fGFPCA, time_fGFPCA,
#      file = here("Data/SimN500/SimOutput_fGFPCA.RData"))
