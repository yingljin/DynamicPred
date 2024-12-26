
# This script includes code for simulation on the reduced sample
# using the competing method GLMMadaptive
# corresponding to manuscript Section 4.4.2

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
library(splines)




#### load and reduce simulated data ####

load(here("Data/sim_data.RData"))

# reduce data
## training data
M <- 500 # number of simulations

sim_data <- lapply(sim_data[1:M], function(x){x %>% filter(id %in% c(1:200))})
# lapply(sim_data, dim)
# head(sim_data[[1]])
# unique(sim_data[[1]]$id)
# unique(sim_data[[1]]$t)

Ntr <- Nte <- 100 # trainig and testing sample size
J <- 1000


#### Data reduction set up ####

t_bin <- seq(0, J, by = 10)

## model fit
## in training, take only 1 observation for every 10 observation
## for example, only the bin midpoints
## It seems that bin width does not affect computational speed too much
## at least much less than model flexibility

#### Containers ####

M <- 5

pred_list_all <- list()
time_vec <- rep(NA, M)

# incorporate a convergence check since some iteration would run into numeric issues
# such as singularity
# we skip such datasets, saving results only of the numeric-free datsets
num_probs <- rep(NA, M)


#### GLMMadaptive ####


pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 
for(m in 1:M){
  
  # data split and training data reduction
  df_train_m <- sim_data[[m]] %>% filter(id %in% 1:100 & sind %in% t_bin) 
  df_test_m <- sim_data[[m]] %>% filter(!id %in% 1:100)
  
  # model fit with a numeric check
  t1 <- Sys.time()
  fit_adglmm <- tryCatch(expr={mixed_model(fixed = Y ~ 0+ns(t, df = 4), 
                                           random = ~ 0+ns(t, df = 4) | id, 
                                           data =  df_train_m, family = binomial())},
                      error = function(e){list(converged = FALSE)},
                      warning = function(w){list(converged = FALSE)}
                      )
  
  # skip datasets with numeric problems
  num_probs[m] <- ifelse(fit_adglmm$converged, 0, 1)  
  if(fit_adglmm$converged){
    
  # if the model was successfully fit, then
      # fit_time[m] <- t2-t1 ## minutes
      
      # prediction
      # t1 <- Sys.time()
      pred_m <- list(
        predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.2),
                newdata2 =  df_test_m %>% filter(t>0.2), 
                type = "subject_specific", type_pred = "link",
                se.fit = T, return_newdata = TRUE)$newdata2,
        
        predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.4),
                newdata2 =  df_test_m %>% filter(t>0.4), 
                type = "subject_specific", type_pred = "link",
                se.fit = T, return_newdata = TRUE)$newdata2,
        
        predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.6),
                newdata2 =  df_test_m %>% filter(t>0.6), 
                type = "subject_specific", type_pred = "link",
                se.fit = T, return_newdata = TRUE)$newdata2,
        
        predict(fit_adglmm, newdata = df_test_m %>% filter(t <= 0.8),
                newdata2 =  df_test_m %>% filter(t>0.8), 
                type = "subject_specific", type_pred = "link",
                se.fit = T, return_newdata = TRUE)$newdata2)
      t2 <- Sys.time()
      time_vec[m] <- difftime(t2, t1, units = "mins") # seconds
      
      pred_m <- lapply(pred_m, function(x){x %>% select(id, t, pred, low, upp)})
      pred_list_all[[m]] <- df_test_m %>% 
        left_join(pred_m[[1]], by = c("id", "t")) %>% 
        rename(pred0.2 = pred,
               pred_lb0.2 = low, 
               pred_ub0.2 = upp) %>% 
        left_join(pred_m[[2]], by = c("id", "t")) %>% 
        rename(pred0.4 = pred,
               pred_lb0.4 = low, 
               pred_ub0.4 = upp) %>% 
        left_join(pred_m[[3]], by = c("id", "t")) %>% 
        rename(pred0.6 = pred,
               pred_lb0.6 = low, 
               pred_ub0.6 = upp) %>% 
        left_join(pred_m[[4]], by = c("id", "t")) %>% 
        rename(pred0.8 = pred,
               pred_lb0.8 = low, 
               pred_ub0.8 = upp)
      
      setTxtProgressBar(pb, m)
  }
  
}

close(pb)

##### Check results #####

mean(time_vec)
sum(num_probs)

# pred_list_all[[1]] %>% filter(t>0.39) %>% View()

pred_list_all[[1]] %>% 
  filter(id %in% 101:104) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb0.2, ymax = pred_ub0.2,
                  col = "0.2", fill = "0.2", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb0.4, ymax = pred_ub0.4,
                  col = "0.4", fill = "0.4", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb0.6, ymax = pred_ub0.6,
                  col = "0.6", fill = "0.6", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred_lb0.8, ymax = pred_ub0.8,
                  col = "0.8", fill = "0.8", alpha = 0.1),
              linetype="dashed")+
  facet_grid(rows = vars(id))+
  guides(alpha = "none", col="none")+
  labs(title = "GLMMadaptive")

#### Save results ####

fit_time_subset_adglmm <- time_vec
pred_subset_adglmm <- pred_list_all

save(fit_time_subset_adglmm, pred_subset_adglmm, 
     num_probs,
     file = here("Data/SubSimOutput_GLMMadaptive.RData"))






