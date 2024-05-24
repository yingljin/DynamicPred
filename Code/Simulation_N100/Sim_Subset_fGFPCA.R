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
library(mvtnorm)
library(splines)




#### load simulated data ####

load(here("Data/sim_data.RData"))

M <- 500 # number of simulations

# reduce data
sim_data <- lapply(sim_data[1:M], function(x){x %>% filter(id %in% c(1:200))})
lapply(sim_data, dim)
head(sim_data[[1]])
unique(sim_data[[1]]$id)
unique(sim_data[[1]]$t)

Ntr <- Nte <- 100 # trainig and testing sample size
J <- 1000



#### set up ####
K <- 4 # number of eigenfunctions
source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))

# values for binning 
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points on time index
mid <- (brks+bin_w/2)[1:n_bin] # interval mid point time index
t <- unique(sim_data[[1]]$t)
mid_t <- t[mid]

# values used for projection
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)
## Projection
B <- spline.des(knots = knots_values, x = mid_t, ord = p + 1,
                outer.ok = TRUE)$design  # evaluate B-splines on binned grid
Bnew <- spline.des(knots = knots_values, x = t, ord = p + 1,
                   outer.ok = TRUE)$design  # evaluate B-splines on original grid

window <- seq(0.2, 0.8, by = 0.2)

#### containers ####

pred_list_all <- list()
converge_state_list <- list()
fit_time <- pred_time <- rep(NA, M)

# M <- 10
M

#### fGFPCA

pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 
for(m in 1:M){
  
  # data split
  df_train_m <- sim_data[[m]] %>% filter(id %in% 1:100)
  df_test_m <- sim_data[[m]] %>% filter(!id %in% 1:100)
  
  t1 <- Sys.time()
  # step 1: bin 
  df_train_m$bin <- cut(df_train_m$sind, breaks = brks, 
                        include.lowest = T, labels = mid)
  df_train_m$bin <- as.numeric(as.character(df_train_m$bin))
  
  # step 2: local GLMMs
  ## on training set
  bin_lst <- split(df_train_m, f = df_train_m$bin)
  df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) # singularity issue
  df_est_latent <- bind_rows(df_est_latent) 
  
  ## Step 3: FPCA
  uni_eta_hat <- df_est_latent %>% filter(bin==sind)
  mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                           nrow=Ntr, 
                           ncol=n_bin, byrow = F) 
  fpca_mod <- fpca.face(mat_est_unique, argvals = mid_t, var=T)
  
  ## Step 4: project and debias
  df_phi <- matrix(NA, J, K) 
  for(k in 1:K){
    lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
    df_phi[,k] <- Bnew %*% coef(lm_mod)
  } # project binned eigenfunctions onto the original grid
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
  new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1] # extract re-evaluated mean
  new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
  # rescale
  new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
  new_phi <- as.matrix(new_phi)
  new_lambda <- new_lambda/n_bin
  t2 <- Sys.time()
  fit_time[m] <- t2-t1 # less than 1 second
  
  
  ## Prediction
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
    for(tmax in window){
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
      converge_state_m[i, which(window==tmax)] <- Fit$Converged
      score <- Fit$Summary1[, "Mode"]
      
      # prediction
      eta_pred_out <- new_mu+new_phi%*%score
      df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
    }
    
    pred_list_m[[i]] <- df_i
  }
  t2 <- Sys.time()
  pred_time[m] <- t2-t1
  
  # clean 
  # save predictions
  df_pred_m <- bind_rows(pred_list_m)
  df_pred_m$pred0.2[df_pred_m$t<=0.2] <- NA
  df_pred_m$pred0.4[df_pred_m$t<=0.4] <- NA
  df_pred_m$pred0.6[df_pred_m$t<=0.6] <- NA
  df_pred_m$pred0.8[df_pred_m$t<=0.8] <- NA
  
  pred_list_all[[m]] <- df_pred_m
  converge_state_list[[m]] <- converge_state_m
  
  setTxtProgressBar(pb, m)
}

close(pb)


# results
mean(fit_time, na.rm = T)/60 # in minutes
class(fit_time)
mean(pred_time, na.rm=T)
sum(lapply(converge_state_list, mean)!=1) # all datasets converged
pred_list_all[[1]] %>%
  filter(t > 0.4) %>%
  View()

fit_time_subset_fGFPCA <- fit_time
pred_time_subset_fGFPCA <- pred_time
pred_subset_fGFPCA <- pred_list_all
save(fit_time_subset_fGFPCA, pred_time_subset_fGFPCA, pred_subset_fGFPCA, 
     file = here("Data/SubSimOutput_fGFPCA.RData"))


#### ISE ####
window <- seq(0, 1, by = 0.2)
ise_fgfpca2 <- array(NA, dim = c(length(window)-2, length(window)-2, 10))


for(m in 1:length(pred_subset_fGFPCA)){
  this_df <- pred_subset_fGFPCA[[m]]
  ise_tb_m <- this_df %>%
    mutate(err1 = (pred0.2-eta_i)^2,
           err2 = (pred0.4-eta_i)^2,
           err3 = (pred0.6-eta_i)^2,
           err4 = (pred0.8-eta_i)^2) %>% 
    select(id, t, starts_with("err")) %>% 
    mutate(window = cut(t, breaks = seq(0, 1, by=0.2), include.lowest = T)) %>% 
    group_by(window, id) %>% 
    summarise_at(vars(err1, err2, err3, err4), sum) %>% 
    group_by(window) %>% 
    summarize_at(vars(err1, err2, err3, err4), mean) %>%
    filter(window != "[0,0.2]") %>% 
    select(starts_with("err")) %>% as.matrix()
  ise_fgfpca2[,,m] <- ise_tb_m
}

mean_ise_fgfpca2 <- apply(ise_fgfpca2, c(1, 2), mean)

colnames(mean_ise_fgfpca2) <- c("0.2", "0.4", "0.6", "0.8")


