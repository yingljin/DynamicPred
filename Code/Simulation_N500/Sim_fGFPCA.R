# This script implements the full simulation 
# for dynamic prediction of fGFPCA
# on the simulated dataset generated from Code/GenData.R

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
theme_set(theme_minimal())

my_logis <- function(x){exp(x)/(1+exp(x))}

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
  geom_line(aes(x=t, y=my_logis(eta_i)), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)

#### fGFPCA Model estimation ####

# functions
source(here("Code/GLMM-FPCA.R")) 
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

# values used for projection
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)


# result container
M <- 5
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
        file = here("Code/prediction.stan"),  # Stan program
        data = stanData,    # named list of data
        chains = 2,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0             # no progress shown
      )
      
      scores_tmax <- summary(fit)$summary[1:K, "mean"]
      sd_scores_tmax <- summary(fit)$summary[1:K, "sd"]
      
      # latent function predictions
      eta_pred_out <- new_mu+new_phi%*%scores_tmax
      df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
      
      # prediction interval
      sd_eta <- sqrt((new_phi^2) %*% sd_scores_tmax^2)
      df_i[ , paste0("pred", tmax, "_lb")] <- as.vector(eta_pred_out[, 1]-qnorm(0.975)*sd_eta)
      df_i[ , paste0("pred", tmax, "_ub")] <- as.vector(eta_pred_out[, 1]+qnorm(0.975)*sd_eta)
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

# score_list_m


#### Check output ####

# time 
time_fGFPCA


# prediction
pred_list_all %>% lapply(dim)
head(pred_list_all[[3]])

# figure
rand_id <- sample(test_id, 4)

df_exp <- pred_list_fGFPCA[[3]] %>% 
  filter(id %in% rand_id) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  
df_exp[df_exp$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
df_exp[df_exp$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
df_exp[df_exp$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
df_exp[df_exp$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA

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
ggsave(here("Images/IntervalExp1.jpeg"), height=12, width = 5)  
  
 
  
# another figure
ggarrange(
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.2), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.2_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.2_ub), linetype="dashed",na.rm=T, col = "red"),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.4), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.4_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.4_ub), linetype="dashed",na.rm=T, col = "red"),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.6), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.6_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.6_ub), linetype="dashed",na.rm=T, col = "red"),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.8), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.8_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.8_ub), linetype="dashed",na.rm=T, col = "red"),
  nrow = 1
)
ggsave(here("Images/IntervalExp.jpeg"), height=12, width = 12)






df_exp %>% 
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"), linetype = "dashed")+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), na.rm=T, linewidth = 1.0)+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), na.rm=T, linewidth = 1.0)+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), na.rm=T, linewidth = 1.0)+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), na.rm=T, linewidth = 1.0)+


#### Save results ####

save(pred_list_fGFPCA, time_fGFPCA,
     file = here("Data/TrialRun/SimOutput_fGFPCA.RData"))

#### Calculate ISE ####

# load simulation results
# load(here("Data/SimOutput_fGFPCA.RData"))

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

ise_mat

mean_ise <- apply(ise_mat, c(1, 2), mean)
mean_ise <- data.frame(mean_ise) %>% 
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
  auc_tb <- this_df %>% 
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

## prediction interval
## Maybe start with if the prediction interval covers the truth?

list_cover <- list()

for(m in 1:M){

  list_cover[[m]] <- pred_list_all[[m]] %>%
    mutate(
      cover0.2 = pred0.2_lb<=eta_i & pred0.2_ub>=eta_i,
      cover0.4 = pred0.4_lb<=eta_i & pred0.4_ub>=eta_i,
      cover0.6 = pred0.6_lb<=eta_i & pred0.6_ub>=eta_i,
      cover0.8 = pred0.8_lb<=eta_i & pred0.8_ub>=eta_i
    ) %>% 
    group_by(t) %>% 
    summarize_at(vars(starts_with("cover")), mean)
}
  
### one simulation
list_cover[[137]] %>% 
  mutate(cover0.2=ifelse(t<=0.2, NA, cover0.2),
         cover0.4=ifelse(t<=0.4, NA, cover0.4),
         cover0.6=ifelse(t<=0.6, NA, cover0.6),
         cover0.8=ifelse(t<=0.8, NA, cover0.8)) %>%
  ggplot()+
  geom_line(aes(x=t, y=cover0.2, col="0.2"))+
  geom_line(aes(x=t, y=cover0.4, col="0.4"))+
  geom_line(aes(x=t, y=cover0.6, col="0.6"))+
  geom_line(aes(x=t, y=cover0.8, col="0.8"))+
  labs(x="t", y="coverage rate", title = "One simulation")
ggsave(here("Images/CoverRate.jpeg"), height=5, width = 5)

### All simulation
list_cover %>% 
  bind_rows(.id="iter") %>%
  group_by(t) %>% 
  summarize_at(vars(starts_with("cover")), mean) %>% 
  mutate(cover0.2=ifelse(t<=0.2, NA, cover0.2),
         cover0.4=ifelse(t<=0.4, NA, cover0.4),
         cover0.6=ifelse(t<=0.6, NA, cover0.6),
         cover0.8=ifelse(t<=0.8, NA, cover0.8)) %>%
  ggplot()+
  geom_line(aes(x=t, y=cover0.2, col="0.2"))+
  geom_line(aes(x=t, y=cover0.4, col="0.4"))+
  geom_line(aes(x=t, y=cover0.6, col="0.6"))+
  geom_line(aes(x=t, y=cover0.8, col="0.8"))+
  labs(x="t", y="coverage rate", title = "All simulation")
  
ggsave(here("Images/CoverRate_All.jpeg"), height=5, width = 5)
