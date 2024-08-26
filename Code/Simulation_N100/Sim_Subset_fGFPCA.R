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
# library(LaplacesDemon)
library(mvtnorm)
library(splines)
library(rstan)
library(RColorBrewer)
cols <- c(brewer.pal(4, "Set2"), "#000000") # define a color palette 
names(cols) <- c("0.2", "0.4", "0.6", "0.8", "True")


#### load simulated data ####

load(here("Data/sim_data.RData"))

M <- 500 # number of simulations

# reduce data
sim_data <- lapply(sim_data[1:M], function(x){x %>% filter(id %in% c(1:200))})
# lapply(sim_data, dim)
# head(sim_data[[1]])
# unique(sim_data[[1]]$id)
# unique(sim_data[[1]]$t)

Ntr <- Nte <- 100 # trainig and testing sample size
J <- 1000



#### set up ####
K <- 4 # number of eigenfunctions
source(here("Code/GLMM-FPCA.R")) 
# source(here("Code/OutSampBayes.R"))

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

# M <- 5
pred_list_all <- list()
# converge_state_list <- list()
time_vec <- rep(NA, M)


M

#### fGFPCA ####

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
  # t1 <- Sys.time()
  debias_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                       s(id, by=phi1, bs="re")+
                       s(id, by=phi2, bs="re")+
                       s(id, by=phi3, bs="re")+
                       s(id, by=phi4, bs="re"), 
                     family = binomial, data=df_train_m, 
                     method = "fREML",
                     discrete = TRUE)
  # t2 <- Sys.time()
  new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1] # extract re-evaluated mean
  new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
  # rescale
  new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
  new_phi <- as.matrix(new_phi)
  new_lambda <- new_lambda/n_bin
  # t2 <- Sys.time()
  # fit_time[m] <- t2-t1 # less than 1 second
  
  
  ## Prediction
  # Prediction pf test sample using Laplace approximation
  pred_list_m <- list()
  # converge_state_m <- matrix(NA, Nte, 4)
  # t1 <- Sys.time()
  # per subject
  df_test_m$id <- droplevels(df_test_m$id)
  test_id <- unique(df_test_m$id)
  
  for(i in seq_along(test_id)){
    df_i <- df_test_m %>% filter(id==test_id[i])
    
    # per max obs time
    for(tid in seq_along(window)){
      tmax <- window[tid]
      df_it <- df_i %>% filter(t <= tmax)
      # max_rid <- nrow(df_it)
      
      # into a list
      stanData <- list(
        J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
        efuncs = new_phi, 
        b0 = new_mu,
        lambda = new_lambda
      )
      
      
      # fit stan
      fit <- stan(
        file = here("Code/prediction.stan"),  # Stan program
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
  
  time_vec[m] <- difftime(t2, t1, units = "mins")
  
  # clean 
  pred_list_all[[m]]<- bind_rows(pred_list_m)
  
  setTxtProgressBar(pb, m)
}

close(pb)


#### check result ####
# simulation is finished
# results
mean(time_vec)


pred_list_all[[1]] %>% 
  filter(t > 0.4) %>%
  View()

rand_id <- sample(test_id, 4)


df_plot <- pred_list_all[[1]] %>% filter(id %in% rand_id) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            .funs = function(x){exp(x)/(1+exp(x))}) 

remove_obs <- function(x){
   x[x$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
   x[x$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
   x[x$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
   x[x$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA
   return(x)}



remove_obs(df_plot) %>%
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
  guides(alpha = "none", col="none")+
  facet_wrap(~id)+
  labs(col = "Maximum observation time", x = "Time", y="")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)

#### save results ####

length(time_vec)
time_subset_fGFPCA <- time_vec
# pred_time_subset_fGFPCA <- pred_time
pred_subset_fGFPCA <- pred_list_all

save(time_subset_fGFPCA, pred_subset_fGFPCA, 
     file = here("Data/SimN100/SubSimOutput_fGFPCA.RData"))


#### ISE ####
window <- seq(0, 1, by = 0.2)
ise_fgfpca2 <- array(NA, dim = c(length(window)-2, length(window)-2, M))


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


