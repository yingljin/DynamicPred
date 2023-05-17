# data application of NHANES

set.seed(516)

library(here)
library(tidyverse)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
theme_set(theme_minimal())

# code
source(here("Code/GLMM-FPCA.R")) # use pred_latent function to estimate latent function 
# source(here("Code/OutSampMLE.R"))
source(here("Code/OutsampBayes.R"))

# data 
df <- read_rds(here("Data/nhanes_bi.rds"))

N <- length(unique(df$SEQN)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject



#### visulization ####
rand_id <- sample(unique(df$SEQN), size = 4)

df %>% 
  filter(SEQN %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Z), size = 0.5)+
  facet_wrap(~SEQN)+
  labs(x="Time", y = "Activity")+
  theme(plot.margin = margin(0, 1, 0, 0))
ggsave(here("Images/NHANES/data.png"), width=7, height=5, bg="white")



##### fGFPCA #####

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points


df$bin <- cut(df$sind, breaks = brks, include.lowest = T, labels = mid)
df$bin <- as.numeric(as.character(df$bin))
head(df, 20)

df <- df %>% rename(id = SEQN, Y=Z)

df_bin_lst <- split(df, f = df$bin)

lapply(df_bin_lst, function(x)table(x$Y))

# fit local GLMM and estimate latent function
# near-unidentifiability issues 
df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
lapply(df_est_latent, function(x)range(x$eta_hat))

# on the binned grid, one unique value each bin
df_est_latent <- bind_rows(df_est_latent) %>% select(id, bin, eta_hat) %>% distinct() 

df_est_latent %>%
  left_join(df %>% filter(sind==bin), by = c("id", "bin")) %>%
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_line(aes(x=bin, y= exp(eta_hat)/(1+exp(eta_hat))), col = "Red")+
  geom_point(aes(x=bin, y = Y), size = 0.5)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")
ggsave(here("Images/NHANES/est_eta.png"), width=7, height=5, bg="white")



# fPCA
mat_est_unique <- matrix(df_est_latent$eta_hat, nrow=N, ncol=n_bin, byrow = F) # row index subject, column binned time
dim(mat_est_unique)

fpca_mod <- fpca.face(mat_est_unique, pve = 0.95, argvals = mid, var=T)

# estimator parameters
dim(fpca_mod$efunctions)

K <- ncol(fpca_mod$efunctions) # 13 eigenfunctions
# Out-of-sample prediction
# observations time: 360, 720, 1080
df_est_latent[, 'pred_t360'] <- df_est_latent[, 'pred_t720'] <- df_est_latent[, 'pred_t1080'] <- NA 
score_out_mat <- array(NA, dim = c(N, K, 3)) # dim indexes subject, eigenfunction and max obs time respectively
dim(score_out_mat)

# unique id

id_vec <- unique(df_est_latent$id)

skip_id <- c()

pb = txtProgressBar(min = 0, max = length(id_vec), initial = 0, style = 3) 

# prediction for a single subject
for(i in seq_along(id_vec)){
  df_i <- df %>% filter(id==id_vec[i]) 
  
  # prediction 
  ## up to 360
  pred_t1 <- tryCatch({out_pred_laplace(fpca_mod, df %>% filter(id==id_vec[i] & sind<=360))},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t1))==0){
    eta_pred_t1 <- fpca_mod$mu+fpca_mod$efunctions%*%pred_t1$score_out
    df_est_latent[df_est_latent$id == id_vec[i], 'pred_t360'] <- eta_pred_t1
    score_out_mat[i, ,1] <- pred_t1$score_out
  } else{
    skip_id <- append(skip_id, id_vec[i])
  }
  
  ## up to 720
  pred_t2 <- tryCatch({out_pred_laplace(fpca_mod, df %>% filter(id==id_vec[i] & sind<=720))},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t2))==0){
    eta_pred_t2 <- fpca_mod$mu+fpca_mod$efunctions%*%pred_t2$score_out
    df_est_latent[df_est_latent$id == id_vec[i], 'pred_t720'] <- eta_pred_t2
    score_out_mat[i, ,2] <- pred_t2$score_out
  } else{
    skip <- c(skip, id_vec[i])
  }
  
  ## up to 1080
  pred_t3 <- tryCatch({out_pred_laplace(fpca_mod, df %>% filter(id==id_vec[i] & sind<=1080))},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t3))==0){
    eta_pred_t3 <- fpca_mod$mu+fpca_mod$efunctions%*%pred_t3$score_out
    df_est_latent[df_est_latent$id == id_vec[i], 'pred_t1080'] <- eta_pred_t3
    score_out_mat[i, ,3] <- pred_t3$score_out
  } else{
    skip <- c(skip, id_vec[i])
  }
  
  setTxtProgressBar(pb, i)
}


#### check results ####
# score
score_out_mat[1,,]

# prediction on probability scale

df_est_latent %>%
  filter(id %in% rand_id) %>%
  mutate(pred_t1080 = ifelse(bin<=1080, NA, pred_t1080),
         pred_t720 = ifelse(bin<=720, NA, pred_t720),
         pred_t360 = ifelse(bin<=360, NA, pred_t360)) %>%
  mutate_at(vars(pred_t360, pred_t720, pred_t1080), function(x)exp(x)/(1+exp(x))) %>%
  left_join(df %>% filter(sind==bin), by = c("id", "bin")) %>%
  ggplot()+
    geom_line(aes(x=bin, y = pred_t360, col = "360"), linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t720, col = "720"),linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t1080, col = "1080"),linetype = "dashed")+
    geom_point(aes(x=bin, y = Y, col = "Outcome"), size = 0.5)+
    facet_wrap(~id)+
    labs(x = "Time", y = "Estimated latent function (probablity scale)")
ggsave(here("Images/NHANES/pred_eta_prob.png"), width=7, height=5, bg="white")

# prediction on latent function scale
df_est_latent %>%
  filter(id %in% rand_id) %>%
  mutate(pred_t1080 = ifelse(bin<=1080, NA, pred_t1080),
         pred_t720 = ifelse(bin<=720, NA, pred_t720),
         pred_t360 = ifelse(bin<=360, NA, pred_t360)) %>%
  ggplot()+
    geom_line(aes(x=bin, y = pred_t360, col = "360"), linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t720, col = "720"),linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t1080, col = "1080"),linetype = "dashed")+
    geom_line(aes(x=bin, y = eta_hat, col = "eta_hat"), linetype = "solid" )+
    facet_wrap(~id)+
    labs(x = "Time", y = "Estimated latent function")
ggsave(here("Images/NHANES/pred_eta_latent.png"), width=7, height=5, bg="white")


# save results
save(df_est_latent, score_out_mat, skip_id, file = here("Data/ApplOutput_fGFPCA.RData"))


#### Numeric issues with laplace approximation ####

load(here("Data/ApplOutput_fGFPCA.RData"))

# seems like the numeric issues happy when we don't have a long enough observed track

df %>% filter(SEQN %in% skip_id[1:4]) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Z))+
  facet_wrap(~SEQN)
head(df)
# compare between subjects with or without numeric issues

com_id <- c(rand_id, skip_id[1:4])



df_est_latent %>% filter(id %in% com_id) %>%
  left_join(df %>% select(SEQN, Z, sind) %>%
              rename(id = SEQN, Y=Z, bin = sind), by = c("id", "bin")) %>%
  mutate(type = ifelse(id %in% rand_id, "Prediction", "Failed prediction")) %>%
  ggplot(aes(x=bin, y=eta_hat, col = type))+
  geom_line()+
  facet_wrap(~id, ncol = 2, nrow = 4)
ggsave(here("Images/NHANES/num_issue.png"), width=7, height=5, bg="white")
  
mis_col_id <- rep(NA, length(skip_id))  
for(i in seq_along(skip_id)){
  
  num_NA <- apply(df_est_latent%>%filter(id==skip_id[i]), 2, function(x){sum(is.na(x))})
  print(num_NA)
  
}

sum(mis_col_id!="pred_t360")
