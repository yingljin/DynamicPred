# data application of NHANES

set.seed(404)

library(here)
library(tidyverse)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
theme_set(theme_minimal())

# code
source(here("Code/GLMM-FPCA.R")) # use pred_latent function to estimate latent function 
source(here("Code/OutSampMLE.R"))

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
  labs(x="Time", y = "Activity")


##### fGFPCA #####

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points


df$bin <- cut(df$sind, breaks = brks, include.lowest = T, labels = mid)
df$bin <- as.numeric(as.character(df$bin))

df <- df %>% rename(id = SEQN, Y=Z)

df_bin_lst <- split(df, f = df$bin)

lapply(df_bin_lst, function(x)table(x$Y))

# fit local GLMM and estimate latent function
# near-unidentifiability issues !!!
df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
df_est_latent <- bind_rows(df_est_latent) %>% select(id, bin, eta_hat) %>% distinct() # on the binned grid, one unique value each bin
df_est_latent %>%
  left_join(df, by = c("id", "bin")) %>%
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_line(aes(x=bin, y= exp(eta_hat)/(1+exp(eta_hat))), col = "Red")+
  geom_point(aes(x=bin, y = Y), size = 0.5)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

# fPCA
mat_est_unique <- matrix(df_est_latent$eta_hat, nrow=N, ncol=n_bin, byrow = F) # row index subject, column binned time
dim(mat_est_unique)

fpca_mod <- fpca.face(mat_est_unique, pve = 0.95, argvals = mid, knots=20, var=T)
K <- ncol(fpca_mod$efunctions) # number of eigenfunctions

# Out-of-sample prediction
# observations time: 360, 720, 1080
df_est_latent[, 'pred_t360'] <- df_est_latent[, 'pred_t720'] <- df_est_latent[, 'pred_t1080'] <- NA 
score_out_mat <- array(NA, dim = c(N, K, 3)) # dim indexes subject, eigenfunction and max obs time respectively


# unique id

id_vec <- unique(df_est_latent$id)

# prediction for a single subject
for(i in seq_along(rand_id)){
  df_i <- df %>% filter(id==rand_id[i]) 
  # prediction 
  pred_t1 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind<=360), fpca_fit = fpca_mod, K = K)
  pred_t2 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind<=720), fpca_fit = fpca_mod, K = K)
  pred_t3 <- out_samp_dyn_pred(df_new = df_i %>% filter(sind<=1080), fpca_fit = fpca_mod, K = K)
  # prediction in container
  df_est_latent[df_est_latent$id == rand_id[i], 'pred_t360'] <- pred_t1$eta_pred
  df_est_latent[df_est_latent$id== rand_id[i], 'pred_t720'] <- pred_t2$eta_pred
  df_est_latent[df_est_latent$id== rand_id[i], 'pred_t1080'] <- pred_t3$eta_pred
  # score in container
  score_out_mat[i, ,1] <- pred_t1$score_out
  score_out_mat[i, ,2] <- pred_t2$score_out
  score_out_mat[i, ,3] <- pred_t3$score_out
}



save(df_est_latent, score_out_mat, file = here("Data/ApplOutput_fGFPCA.RData"))

df_est_latent %>%
  filter(id %in% rand_id) %>%
  mutate_at(vars(pred_t360, pred_t720, pred_t1080), function(x)exp(x)/(1+exp(x))) %>%
  # filter(id %in% rand_id) %>%
  ggplot()+
  geom_line(aes(x=bin, y = pred_t360, col = "360"))+
  geom_line(aes(x=bin, y = pred_t720, col = "720"))+
  geom_line(aes(x=bin, y = pred_t1080, col = "1080"))+
  # geom_point(aes(x=bin, y = ), size = 0.5)+
   facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")
