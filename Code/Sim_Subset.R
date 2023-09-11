set.seed(314)

library(here)
library(tidyverse)
# library(ggplot2)
# theme_set(theme_minimal())
# library(fcr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(mgcv)
# library(ggpubr)
# library(gridExtra)
# library(kableExtra)




#### load simulated data ####

load(here("Data/sim_data.RData"))


#### simulation set up ####

J <- 1000 # number of observations points
N <- 500 # sample size
M <- 50 


#### GLMM adaptive ####

# cubic spline basis functions
df_basis <- smoothCon(s(sind, bs = "cr", k = 20), data=sim_data[[1]])[[1]]$X
colnames(df_basis) <- paste("B", 1:20, sep="")
this_df <- data.frame(sim_data[[1]] %>% select(id, Y), df_basis)
fix_f <- paste("Y ~", 
                paste("B", 1:20, sep = "", collapse = " + "))
fix_f <- as.formula(fix_f)
rand_f <- paste

adglmm <- mixed_model(this_f, random = ~ 1 | id, data =  this_df, family = binomial())


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
