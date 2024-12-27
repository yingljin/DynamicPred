# data application of NHANES
# using fGFPCA
# corresponding to manuscript Section 5


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
library(LaplacesDemon)
library(splines)
library(mvtnorm)
library(rstan)

#### load data ####
df <- read_rds(here("DataRaw/nhanes_bi_sub.rds"))
df <- df %>% rename(id=SEQN, Y=Z)

N <- length(unique(df$id)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject
t <- unique(df$sind)
K <- 4 # number of eigenfunctions to use

# functions
source(here("Code/Functions/GLMM-FPCA.R")) 


#### Data split ####

# 60% (5257) subjects for training, 40% (3506) for out-of-sample prediction
train_id <- sample(unique(df$id), size = N*0.6)
test_id <- setdiff(unique(df$id), train_id)

train_df <- df %>% filter(id %in% train_id)
test_df <- df %>% filter(id %in% test_id)


#### fGFPCA ####

# Step 1: Bin every 10 observations
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

train_df$bin <- cut(train_df$sind, breaks = brks, include.lowest = T, labels = mid)
train_df$bin <- as.numeric(as.character(train_df$bin))

## checks
head(train_df)
table(train_df$bin)


# Step 2:  Local GLMM
train_bin_lst <- split(train_df, f = train_df$bin)

t1=Sys.time()
df_est_latent <- lapply(train_bin_lst, function(x){pred_latent(x, n_node = 0)}) 
t2= Sys.time()
t2-t1 

## example estimated latent function
df_est_latent <- bind_rows(df_est_latent) 
rand_id <- sample(train_id, 4)
df_est_latent %>% 
  filter(id %in% rand_id) %>%
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_line(aes(x=sind, y=eta_hat, group = id, col = "estimated"))+
  geom_point(aes(x=sind, y = Y, group = id, col = "outcome"), size = 0.5)+
  facet_wrap(~id, scales = "free")+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")



# Step 3: FPCA
uni_eta_hat <- df_est_latent %>% filter(bin==sind)
mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                         nrow=length(train_id), 
                         ncol=n_bin, byrow = F) 
# row index subject, column binned time

t1 <- Sys.time()
fpca_mod <- fpca.face(mat_est_unique, argvals = mid, var=T, npc=K) # keep 4 PCs
t2 <- Sys.time()
t2-t1 

## check eigenvaluse and eigenfunctions
fpca_mod$evalues
data.frame(sind = mid, 
           fpca_mod$efunctions) %>%
  pivot_longer(2:5) %>%
  ggplot()+
  geom_line(aes(x=sind, y=value))+
  facet_wrap(~name)

## check mean
plot(mid, fpca_mod$mu, type = "l")

# save(fpca_mod, file = here("Data/Appl_fpca_model.RData"))

# Step 4: Re-evaluation
## grid extension
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid) - min(mid)) + min(mid)

B <- spline.des(knots = knots_values, x = mid, ord = p + 1,
                outer.ok = TRUE)$design  # evaluate B-splines on binned grid
Bnew <- spline.des(knots = knots_values, x = 1:J, ord = p + 1,
                   outer.ok = TRUE)$design  # evaluate B-splines on original grid

df_phi <- matrix(NA, J, K) 
for(k in 1:K){
  lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
  df_phi[,k] <- Bnew %*% coef(lm_mod)
}# project binned eigenfunctions onto the original grid

# plot to compare eigenfunctions before and after re-evaluation
df_pc1 <- data.frame(t=1:J,  df_phi) %>%
  rename(PC1=X1, PC2=X2, PC3=X3, PC4=X4) %>%
  pivot_longer(2:5, names_to = "PC") 

df_pc2 <- data.frame(t = mid, fpca_mod$efunctions) %>%
  rename(PC1=X1, PC2=X2, PC3=X3, PC4=X4) %>%
  pivot_longer(2:5, names_to = "PC") 

left_join(df_pc1, df_pc2, by = c("t", "PC")) %>%
  ggplot()+
  geom_line(aes(x=t, y=value.x, col = PC), linewidth = 0.2)+
  geom_point(aes(x=t, y=value.y, col = PC), na.rm = T, alpha = 0.5, size = 0.5)+
  facet_wrap(~PC)+
  labs(x="Time", y="", title = "Re-evaluate eigenfunctions on the original grid")

## debias
df_phi <- data.frame(sind = 1:J, df_phi)
colnames(df_phi) <- c("sind", paste0("phi", 1:K))
train_df <- train_df %>% 
  select(!starts_with("phi")) %>%
  left_join(df_phi, by = "sind")
train_df$id <- as.factor(train_df$id)

head(train_df)
head(df_phi)

# don't run this part unless it's absolutely necessary!
# It takes 20 hours!
t1 <- Sys.time()
debias_glmm <- bam(Y ~ s(sind, bs = "bs")+
                     s(id, by=phi1, bs="re")+
                     s(id, by=phi2, bs="re")+
                     s(id, by=phi3, bs="re")+
                     s(id, by=phi4, bs="re"),
                   family = binomial,
                   data=train_df,
                   method = "fREML",
                   discrete = TRUE)
t2 <- Sys.time()
t2-t1 # 21 hours
# save model
# save(debias_glmm, file = here("Data/DataAppl/global_bam.RData"))

# if we don't wish to retrain the model, we can simply reload the debias GLMM model
# load(here("Data/DataAppl/global_bam.RData"))

# check mean function
new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1]# extract re-evaluated mean
plot(t, new_mu)
lines(mid, fpca_mod$mu, col = "Red")
summary(new_mu)

# check eigenvalues
new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
new_lambda/sum(new_lambda)
# the 2nd and 3rd eigenvalues flipped! 

## rescale eigenfunctions using the number of bins
head(df_phi) 
new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
new_phi <- as.matrix(new_phi)
head(new_phi)

new_lambda <- new_lambda/n_bin


#### prediction ####
# test_id <- test_id[1:500]

window <- seq(0, J, by = 360) 

pred_list_m <- list()

pb = txtProgressBar(min = 0, max = length(test_id), initial = 0, style = 3) 

t1 <- Sys.time()
for(i in seq_along(test_id)){
  df_i <- test_df %>% filter(id==test_id[i])
  
  # per max obs time
  # t1 <- Sys.time()
  for(tmax in window[2:4]){
    df_it <- df_i %>% filter(t <= tmax)
    
    # into a list
    stanData <- list(
      J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
      efuncs = new_phi, 
      b0 = new_mu,
      lambda = new_lambda
    )
    
    # fit stan model
    fit <- stan(
      file = here("Code/Functions/prediction.stan"),  # Stan program
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
  # t2 <- Sys.time()
  
  pred_list_m[[i]] <- df_i
  
  setTxtProgressBar(pb, i)
}
t2 <- Sys.time()
close(pb)



#### Results ####

head(pred_list_m[[1]])
length(pred_list_m)
pred_nhanes_fgfpca <- bind_rows(pred_list_m)

#### save results ####
save(pred_nhanes_fgfpca, 
     file = here("Data/ApplOutput_fGFPCA.RData"))

#### Check results ####

# Prediction tracks
head(pred_nhanes_fgfpca)
pred_nhanes_fgfpca$pred360[pred_nhanes_fgfpca$sind<=360] <- NA
pred_nhanes_fgfpca$pred720[pred_nhanes_fgfpca$sind<=720] <- NA
pred_nhanes_fgfpca$pred1080[pred_nhanes_fgfpca$sind<=1080] <- NA


pred_nhanes_fgfpca %>% 
  select(!ends_with("b")) %>% 
  filter(id %in% sample(unique(pred_nhanes_fgfpca$id), 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred360, col = "6am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred720, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred1080, col = "6pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")


