# This script implements the full simulation 
# for dynamic prediction of fGFPCA

#### set up ####

set.seed(516)

library(here)
library(tidyverse)
theme_set(theme_minimal())

library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(mgcv)
library(splines)



#### load simulated data ####

load(here("Data/sim_data.RData"))

# data structure
M <- length(sim_data) # number of simulation
N <- length(unique(sim_data[[1]]$id)) # sample size 
J <- length(unique(sim_data[[1]]$sind)) # number of measurements for each subject

# data preview
rand_id <- sample(as.numeric(unique(sim_data[[1]]$id)), size=4)

sim_data[[1]] %>% 
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.5)+
  geom_line(aes(x=sind, y=exp(eta_i)/(1+exp(eta_i))))+
  facet_wrap(~id)


##### fGFPCA Model estimation #####

source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))
# use pred_latent function to estimate latent function 

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

# result container
pred_list_all <- list()
# score <- array()
# score_lst <- list()
fit_time <- pred_time <- rep(NA, M)

K <- 4 # number of eigenfunctions to use

pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

# The whole process
for(m in 1:M){
  
  # for every simulated dataset
  df <- sim_data[[m]]
  
  ## Step 1: bin observations
  df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
  df$bin <- as.numeric(as.character(df$bin))
  
  # split into training and test set
  # 80% subjects for training, 20% for testing
  train_id <- sample(unique(df$id), size = N*0.8)
  test_id <- setdiff(unique(df$id), train_id)
  train_df <- df %>% filter(id %in% train_id)
  test_df <- df %>% filter(id %in% test_id)
  
  # fit fGFPCA model
  t1 <- Sys.time()
  ## Step 2: local GLMM and estimate latent function
  ## on training set
  bin_lst <- split(train_df, f = train_df$bin)
  df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) 
  df_est_latent <- bind_rows(df_est_latent) 
  ## Step 3: FPCA
  mat_est_unique <- matrix(df_est_latent$eta_hat[df_est_latent$sind_inx==df_est_latent$bin],
                           nrow=N*0.8, ncol=n_bin, byrow = F) # row-subject, column-time
  fpca_mod <- fpca.face(mat_est_unique, argvals = mid, knots=20, var=T)
  
  ## Step 4: project and debias
  ## We have not yet decided on project method, so let's do debias on the binned grid for now
  df_phi <- data.frame(bin = mid, fpca_mod$efunctions[, 1:K])
  colnames(df_phi) <- c("bin", paste0("phi", 1:4))
  train_df <- train_df %>% left_join(df_phi, by = "bin")
  train_df$id <- as.factor(train_df$id)
  debias_glmm <- bam(Y ~ s(bin, bs="cr", k=10)+
                       s(id, by=phi1, bs="re", k=10)+
                       s(id, by=phi2, bs="re", k=10)+
                       s(id, by=phi3, bs="re", k=10)+
                       s(id, by=phi4, bs="re", k=10), 
                     family = binomial, 
                     data=train_df %>% filter(sind_inx==bin), 
                     method = "fREML",
                     discrete = T)
  t2 <- Sys.time()
  fit_time[m] <- difftime(t2, t1, units = "secs")
  ## get re-evaluated eigenvaleus
  new_evals <- 1/debias_glmm$sp[2:5]

  # Out-of-sample prediction with Laplace approximation
  ## Observation track: 0-200, 0-400, 0-600, 0-800
  ## prediction window: 200-400, 400-600, 600-800, 800-1000
  ### containers
  pred_list_m <- list()
  # score_out_mat <- array(NA, dim = c(N, K, 4)) # dim: subject, eigenfunction, obs track
  
  t1 <- Sys.time()
  for(i in seq_along(test_id)){
    df_i <- test_df %>% filter(id==test_id[i])
    # prediction 
    pred1 <- out_pred_laplace(mu = fpca_mod$mu, 
                              evalues = new_evals, 
                              phi_mat = fpca_mod$efunctions[, 1:K], 
                              df_new = df_i %>% filter(sind<=200), kpc = K)
    pred2 <- out_pred_laplace(mu = fpca_mod$mu, 
                              evalues = fpca_mod$evalues[1:K], 
                              phi_mat = fpca_mod$efunctions[, 1:K],
                              df_i %>% filter(sind<=400), kpc = K)
    pred3 <- out_pred_laplace(mu = fpca_mod$mu, 
                              evalues = fpca_mod$evalues[1:K], 
                              phi_mat = fpca_mod$efunctions[, 1:K],
                              df_i %>% filter(sind<=600), kpc = K)
    pred4 <- out_pred_laplace(mu = fpca_mod$mu, 
                              evalues = fpca_mod$evalues[1:K], 
                              phi_mat = fpca_mod$efunctions[, 1:K],
                              df_i %>% filter(sind<=800), kpc = K)
    
    
    pred_i <- data.frame(id = test_id[i], bin=mid, 
                         pred1=pred1$eta_pred, pred2=pred2$eta_pred,
                         pred3=pred3$eta_pred, pred4=pred4$eta_pred)
    
    pred_i$pred1[pred_i$bin<=200] <- NA
    pred_i$pred2[pred_i$bin<=400] <- NA
    pred_i$pred3[pred_i$bin<=600] <- NA
    pred_i$pred4[pred_i$bin<=800] <- NA
    
    pred_list_m[[i]] <- pred_i
  }
  t2 <- Sys.time()
  pred_time[m] <- difftime(t2, t1, units = "secs")
 
  pred_list_all[[m]] <- test_df %>% left_join(bind_rows(pred_list_m), by = c("id", "bin")) 
  
  
  setTxtProgressBar(pb, m)
}

close(pb)




#### Check output ####

mean(fit_time) # average time for model fitting: 26 secs
mean(pred_time) # average time for prediction: 133 secs 

pred_list_all %>% lapply(dim)

rand_id <- sample(pred_list_all[[158]]$id, 4)

# prediction results
pred_list_all[[158]] %>% 
  filter(id %in% sample(pred_list_all[[158]]$id, 4)) %>% 
  mutate_at(vars(eta_i, pred1, pred2, pred3, pred4), function(x){exp(x)/(1+exp(x))}) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.2)+
  geom_line(aes(x=sind, y=eta_i, col = "True"))+
  geom_line(aes(x=sind, y=pred1, col = "200"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred2, col = "400"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred3, col = "600"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=sind, y=pred4, col = "800"), linetype="dashed", na.rm=T)+
  facet_wrap(~id)




#### Save results ####

save(pred_list_all, fit_time, pred_time, file = here("Data/SimOutput_fGFPCA.RData"))







