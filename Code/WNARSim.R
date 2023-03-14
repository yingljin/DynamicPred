# This script implements simulation for WNAR manuscript

set.seed(1025)

library(here)
library(tidyverse)
# library(ggplot2)
# theme_set(theme_minimal())
# library(fcr)
library(refund)
library(lme4)
# library(ggpubr)
# library(gridExtra)
# library(kableExtra)




#### Generate data ####
source(here("Code/ToySimulation.R"))


sim_data <- list()

M <- 500

# generate 500 subject with four eigenfunctions, observation over 1000 time points
for(m in 1:M){
  sim_data[[m]] <- gen_data(N = 500, J = 1000, sind = seq(0,1,len=1000) , lambda = 0.5^(0:(4-1)))
}

save(sim_data, file = here("Data/sim_data.RData"))

load(here("Data/sim_data.RData"))

#### fGFPCA ####
source(here("Code/GLMM-FPCA.R"))
source(here("Code/OutSampMLE.R"))

J <- 1000
N <- 500
K <- 4


# result container
pred_lst <- list()
score_lst <- list()



t1 <- Sys.time()
for(m in 1:50){
  
  df <- sim_data[[m]]
  
  # bin observations
  bin_w <- 10 # bin width
  n_bin <- J/bin_w # number of bins
  brks <- seq(0, J, by = bin_w) # cutoff points
  mid <- (brks+bin_w/2)[1:n_bin] # mid points
  df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
  df$bin <- as.numeric(as.character(df$bin))
  
  df_bin_lst <- split(df, f = df$bin)
  
  # fit local GLMM and estimate latent function
  df_pred_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
  df_pred_latent <- bind_rows(df_pred_latent)

  # put predicted latent function to wide format
  df_pred_latent <- df_pred_latent %>% arrange(id, sind_inx)
  df_pred_unique <- df_pred_latent %>% select(id, bin, eta_hat) %>% distinct()
  mat_pred_unique <- matrix(df_pred_unique$eta_hat, nrow=N, ncol=n_bin, byrow = T)
  
  # fPCA
  fpca_fit <- fpca.face(mat_pred_unique, pve = 0.95, argvals = unique(df_pred_unique$bin), 
                        knots=20, var=T)
  
  # Out-of-sample prediction
  df_pred_unique[, 'pred_t200'] <- df_pred_unique[, 'pred_t400'] <- df_pred_unique[, 'pred_t600'] <- df_pred_unique[, 'pred_t800'] <- NA 
  score_out_mat <- array(NA, dim = c(N, K, 4))
  
  # prediction for a single subject
  for(i in 1:N){
    df_i <- df %>% filter(id==i) %>% select(-eta_i)
    # prediction 
    pred_t200 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=200), fpca_fit = fpca_fit)
    pred_t400 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=400), fpca_fit = fpca_fit)
    pred_t600 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=600), fpca_fit = fpca_fit)
    pred_t800 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=800), fpca_fit = fpca_fit)
    # prediction in container
    df_pred_unique[df_pred_unique$id==i, 'pred_t200'] <- pred_t200$eta_pred
    df_pred_unique[df_pred_unique$id==i, 'pred_t400'] <- pred_t400$eta_pred
    df_pred_unique[df_pred_unique$id==i, 'pred_t600'] <- pred_t600$eta_pred
    df_pred_unique[df_pred_unique$id==i, 'pred_t800'] <- pred_t800$eta_pred
    # score in container
    score_out_mat[i, ,1] <- pred_t200$score_out
    score_out_mat[i, ,2] <- pred_t400$score_out
    score_out_mat[i, ,3] <- pred_t600$score_out
    score_out_mat[i, ,4] <- pred_t800$score_out
  }
  
  pred_lst[[m]] <- df_pred_unique
  score_lst[[m]] <- score_out_mat
}

t2 <- Sys.time()

RunTime <- (t2-t1)

save(pred_lst, score_lst, RunTime, file = here("Data/SimOutput_fGFPCA.RData"))

#### GLMMadaptive ####

#### Plots #####

pred_lst[[1]] %>% filter(id==5) %>% 
  left_join(df %>% filter(id==5) %>% filter(sind_inx %in% mid) %>% select(bin, eta_i)) %>%
  ggplot()+
  geom_line(aes(x=bin, y=eta_i))+
  geom_line(aes(x=bin, y = pred_t200, col = "200"), linetype = "dashed")+
  geom_line(aes(x=bin, y = pred_t400, col = "400"), linetype = "dashed")+
  geom_line(aes(x=bin, y = pred_t600, col = "600"), linetype = "dashed")+
  geom_line(aes(x=bin, y = pred_t800, col = "800"), linetype = "dashed")
