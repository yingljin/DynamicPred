# This script implements simulation 
# but only for one iterarion
# as an explorative example

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

df <- sim_data[[1]]
source(here("Code/GLMM-FPCA.R")) # use pred_latent function to estimate latent function 
# source(here("Code/OutSampMLE.R"

# Step 1: bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
table(df$bin)
df$bin <- as.numeric(as.character(df$bin))

## preview of binned data
df %>% filter(id %in% rand_id) %>%
  group_by(id, bin) %>%
  summarise(num = sum(Y)) %>%
  ggplot(aes(x=bin, y=num))+
  geom_point(size = 0.5)+
  facet_wrap(~id)

# Step 2: Local GLMMs
df_bin_lst <- split(df, f = df$bin)
lapply(df_bin_lst, dim)

t1 <- Sys.time()
df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x, n_node = 0)}) 
t2 <- Sys.time()
t_local_glmm <- t2-t1 # 4.72 secs

df_est_latent <- bind_rows(df_est_latent) 

## Compare estimated latent function to true latent function
df_est_latent %>% filter(id %in% rand_id) %>%
  ggplot()+
  geom_line(aes(x=sind, y=eta_i, col="True"))+
  geom_line(aes(x=sind, y=eta_hat, col="Estimated"))+
  facet_wrap(~id)

# Step 3: fPCA
est_latent_unique <- df_est_latent %>% select(id, bin, eta_hat) %>% distinct(.) 
mat_est_unique <- matrix(est_latent_unique$eta_hat, 
                         nrow=N, ncol=n_bin, byrow = F) 
  # row index subject, column binned time
fpca_mod <- fpca.face(mat_est_unique, argvals = mid, knots=20, var=T)

## summary of FPCA results
plot(mid, fpca_mod$mu)
fpca_mod$evalues
K <- 4 # number of eigenfunctions to use
plot(mid, fpca_mod$efunctions[, 1])
plot(mid, fpca_mod$efunctions[, 2])
plot(mid, fpca_mod$efunctions[, 3])
plot(mid, fpca_mod$efunctions[, 4])
fpca_mod$VarMats[[1]][1:5, 1:5]
fpca_mod$VarMats[[5]][1:5, 1:5]
length(fpca_mod$VarMats)

heatmap(cov2cor(fpca_mod$VarMats[[1]]), Rowv = NA, Colv = NA) 
  # individual correlation


###### debias? #######

# first, try on the binned grid with mgcv::bam
df_phi <- data.frame(bin=mid, fpca_mod$efunctions[, 1:4])
colnames(df_phi) <- c("bin", paste0("phi", 1:4))

df <- df %>% left_join(df_phi, by = "bin")
df$id <- as.factor(df$id)
# length(unique(df$id))

rm(df_bin_lst, est_latent_unique, mat_est_unique, sim_data)

# usethis::edit_r_environ()

t1 <- Sys.time()
debias_glmm <- bam(Y ~ s(bin, bs="cr", k=10)+
                     s(id, by=phi1, bs="re", k=10)+
                     s(id, by=phi2, bs="re", k=10)+
                     s(id, by=phi3, bs="re", k=10)+
                     s(id, by=phi4, bs="re", k=10), 
                   family = binomial, data=df, method = "fREML", 
                   discrete = TRUE)
t2 <- Sys.time()
t_debias <- t2-t1
## took 3.04 hours to fit, and this is just one iteration!
## This is not realistic. Need to find another way.
## Set discrete=TRUE shortens the time to 36.29 secs! 


## next, extract the re-evaluated eigenvalues
## from the smoothing factors
# https://sites.stat.washington.edu/courses/stat527/s14/readings/Wand_MixedModels.pdf
# https://www.degruyter.com/document/doi/10.1515/ijb-2015-0026/html?lang=en
summary(debias_glmm)
debias_glmm$sp # smoothing parameters
sqrt(1/debias_glmm$sp[2:5]) # not close to the simulation eigenvaluesa at all! 
xi_var_true <- 0.5^(0:(K-1))

## maybe because sind and sind_inx are on different scale? 





# interpolate eigenfunctions back to the original grid? 





######## Out-of-sample Prediction #############


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
# preview of binned data
sim_data[[1]]

