# This script implements simulation for WNAR manuscript

set.seed(516)

library(here)
library(tidyverse)
theme_set(theme_minimal())
# library(fcr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(mgcv)
library(splines)
# library(ggpubr)
# library(gridExtra)
# library(kableExtra)




#### load simulated data ####

load(here("Data/sim_data.RData"))

# data structure
M <- length(sim_data) # number of simulation
N <- length(unique(sim_data[[1]]$id))
J <- length(unique(sim_data[[1]]$sind))

# data preview
rand_id <- sample(as.numeric(unique(sim_data[[1]]$id)), size=4)

sim_data[[1]] %>% 
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.5)+
  geom_line(aes(x=sind, y=exp(eta_i)/(1+exp(eta_i))))+
  facet_wrap(~id)


##### fGFPCA Model estimation #####

source(here("Code/GLMM-FPCA.R")) # use pred_latent function to estimate latent function 
# source(here("Code/OutSampMLE.R"))


# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points


# result container
pred_lst <- list()
score_lst <- list()
runtime <- rep(NA, M)

m <- 1

pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 


for(m in 1:M){
  
  df <- sim_data[[m]]
  
  t1 <- Sys.time()
  # bin observations
  df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
  df$bin <- as.numeric(as.character(df$bin))
  
  df_bin_lst <- split(df, f = df$bin)
  
  # fit local GLMM and estimate latent function
  df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x, n_node = 0)}) 
  df_est_latent <- bind_rows(df_est_latent) %>% select(id, bin, eta_hat) %>% distinct() # on the binned grid, one unique value each bin
  mat_est_unique <- matrix(df_est_latent$eta_hat, nrow=N, ncol=n_bin, byrow = F) # row index subject, column binned time
  
  # fPCA
  fpca_mod <- fpca.face(mat_est_unique, argvals = mid, knots=20, var=T)
  
  ###### debias? #######
  # first, try on the binned grid
  df_phi <- data.frame(bin=mid, fpca_mod$efunctions[, 1:4])
  colnames(df_phi) <- c("bin", paste0("phi", 1:4))
  
  df <- df %>% left_join(df_phi, by = "bin")
  df$id <- as.factor(df$id)
  # length(unique(df$id))
  
  # usethis::edit_r_environ()
  
  t1 <- Sys.time()
  debias_glmm <- bam(Y ~ s(bin, bs="cr", k=10)+
                       s(id, by=phi1, bs="re", k=10)+
                       s(id, by=phi2, bs="re", k=10)+
                       s(id, by=phi3, bs="re", k=10)+
                       s(id, by=phi4, bs="re", k=10), 
                     family = binomial, data=df %>% filter(id %in% sample(N, size = 100)), method = "fREML")
  t2 <- Sys.time()
  t_debias <- t2=t1
  ## still takes a really long time (more than 15 minutes), probably won't finish fitting. 
  ## 10 subjects are doable
  
  
  # then, try on the original grid? 
  #####################
  
  
  # Out-of-sample prediction
  df_est_latent[, 'pred_t200'] <- df_est_latent[, 'pred_t400'] <- df_est_latent[, 'pred_t600'] <- df_est_latent[, 'pred_t800'] <- NA 
  score_out_mat <- array(NA, dim = c(N, K, 4)) # dim indexes subject, eigenfunction and max obs time respectively
  
  # prediction for a single subject
  for(i in 1:N){
    df_i <- df %>% filter(id==i) %>% select(-eta_i)
    # prediction 
    pred_t200 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=200), fpca_fit = fpca_mod)
    pred_t400 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=400), fpca_fit = fpca_mod)
    pred_t600 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=600), fpca_fit = fpca_mod)
    pred_t800 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind_inx<=800), fpca_fit = fpca_mod)
    # prediction in container
    df_est_latent[df_est_latent$id==i, 'pred_t200'] <- pred_t200$eta_pred
    df_est_latent[df_est_latent$id==i, 'pred_t400'] <- pred_t400$eta_pred
    df_est_latent[df_est_latent$id==i, 'pred_t600'] <- pred_t600$eta_pred
    df_est_latent[df_est_latent$id==i, 'pred_t800'] <- pred_t800$eta_pred
    # score in container
    score_out_mat[i, ,1] <- pred_t200$score_out
    score_out_mat[i, ,2] <- pred_t400$score_out
    score_out_mat[i, ,3] <- pred_t600$score_out
    score_out_mat[i, ,4] <- pred_t800$score_out
  }
  t2 <- Sys.time()
  
  pred_lst[[m]] <- df_est_latent
  score_lst[[m]] <- score_out_mat
  runtime[m] <- t2-t1
  
  setTxtProgressBar(pb, m)
}


save(pred_lst, score_lst, runtime, file = here("Data/SimOutput_fGFPCA.RData"))

#### GLMMadaptive ####

pred_lst_ref <- list()
runtime_ref <- rep(NA, M)



pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

for(m in 1:M){
  df <- sim_data[[m]]
  
  t1 <- Sys.time()
  # GLMM adaptive model
  adglmm <- mixed_model(Y ~ sind, random = ~ sind | id, data =  df, family = binomial())
  
  # container
  df_est_latent <- df
  df_est_latent[, 'pred_t200'] <- df_est_latent[, 'pred_t400'] <- df_est_latent[, 'pred_t600'] <- df_est_latent[, 'pred_t800'] <- NA 
  
  # time up to t
  pred_t200 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 200),
                        newdata2 =  df %>% filter(sind_inx > 200), type = "subject_specific", type_pred = "link",
                          se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t200[df_est_latent$sind_inx>200] <- pred_t200$newdata2$pred
  
  pred_t400 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 400),
                       newdata2 =  df %>% filter(sind_inx > 400), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t400[df_est_latent$sind_inx>400] <- pred_t400$newdata2$pred
  
  pred_t600 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 600),
                       newdata2 =  df %>% filter(sind_inx > 600), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t600[df_est_latent$sind_inx>600] <- pred_t600$newdata2$pred
  
  pred_t800 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 800),
                       newdata2 =  df %>% filter(sind_inx > 800), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t800[df_est_latent$sind_inx>800] <- pred_t800$newdata2$pred
  
  t2 <- Sys.time()
  
  
  pred_lst_ref[[m]] <- df_est_latent
  runtime_ref[m] <- t2-t1
  
  setTxtProgressBar(pb, m)
  
}

save(pred_lst_ref, runtime_ref, file = here("Data/SimOutput_GLMMadaptive.RData"))

#### Plots #####

par(mfrow=c(2, 2))
plot(mid, fpca_mod$efunctions[,1])
plot(mid, fpca_mod$efunctions[,2])
plot(mid, fpca_mod$efunctions[,3])
plot(mid, fpca_mod$efunctions[,4])



pred_lst_ref[[1]] %>% filter(id==5) %>%
  left_join(df %>% filter(id==5) %>%  select(sind, eta_i)) %>%
  ggplot()+
  geom_line(aes(x=sind, y=eta_i))+
  geom_line(aes(x=sind, y = pred_t200, col = "200"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t400, col = "400"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t600, col = "600"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t800, col = "800"), linetype = "dashed")


mean(runtime_ref)
