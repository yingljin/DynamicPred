# data application of NHANES

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

#### load data ####
df <- read_rds(here("Data/nhanes_bi.rds"))
df <- df %>% rename(id=SEQN, Y=Z)

N <- length(unique(df$id)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject
t <- unique(df$sind)
K <- 4 # number of eigenfunctions to use

# functions
source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))

#### Data split ####

# 60% (5257) subjects for training, 40% (3506) for out-of-sample prediction
train_id <- sample(unique(df$id), size = N*0.6)
test_id <- setdiff(unique(df$id), train_id)

train_df <- df %>% filter(id %in% train_id)
test_df <- df %>% filter(id %in% test_id)


##### fGFPCA #####

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



# FPCA
uni_eta_hat <- df_est_latent %>% filter(bin==sind)
mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                         nrow=length(train_id), 
                         ncol=n_bin, byrow = F) 
# row index subject, column binned time

t1 <- Sys.time()
fpca_mod <- fpca.face(mat_est_unique, argvals = mid, var=T, npc=K) # keep 4 PCs
t2 <- Sys.time()
t2-t1 

## checks
fpca_mod$evalues
data.frame(sind = mid, 
           fpca_mod$efunctions) %>%
  pivot_longer(2:5) %>%
  ggplot()+
  geom_line(aes(x=sind, y=value))+
  facet_wrap(~name)

plot(mid, fpca_mod$mu, type = "l")

save(fpca_mod, file = here("Data/Appl_fpca_model.RData"))

# Re-evaluation
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
# save(debias_glmm, file = here("Data/Appl_debias_model.RData"))


# reload the debias GLMM model
load(here("Data/Appl_debias_model.RData"))

new_mu <- predict(debias_glmm, type = "terms")[1:J, 1]+coef(debias_glmm)[1]# extract re-evaluated mean
plot(t, new_mu)
lines(mid, fpca_mod$mu, col = "Red")
summary(new_mu)

# fpca_mod$evalues
new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda
new_lambda/sum(new_lambda)
# the 2nd and 3rd eigenvalues flipped! 

## rescale using the number of bins
head(df_phi) 
new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
new_phi <- as.matrix(new_phi)
head(new_phi)

new_lambda <- new_lambda/n_bin


# prediction
# test_id <- test_id[1:500]

window <- seq(0, J, by = 360) 

converge_state_m <- matrix(NA, nrow = length(test_id), 4)
pred_list_m <- list()

t1 <- Sys.time()
for(i in seq_along(test_id)){
  df_i <- test_df %>% filter(id==test_id[i])
  
  # per max obs time
  for(tmax in window[2:5]){
    df_it <- df_i %>% filter(t <= tmax)
    max_rid <- nrow(df_it)
    
    # into a list
    MyData <- list(K=K, 
                   X=new_phi[1:max_rid, ], 
                   mon.names=mon.names,
                   parm.names=parm.names, 
                   pos.xi=pos.xi, 
                   y=df_it$Y, 
                   tao=diag(new_lambda), f0=new_mu[1:max_rid])
    
    
    # fit laplace approximation
    Fit <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData, Method = "NM", 
                                Iterations = 1000,
                                CovEst = "Identity")
    converge_state_m[i, which(window[2:6]==tmax)] <- Fit$Converged
    score <- Fit$Summary1[, "Mode"]
    
    # prediction
    eta_pred_out <- new_mu+new_phi%*%score
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
  }
  
  pred_list_m[[i]] <- df_i
}
t2 <- Sys.time()
t2-t1 



# numeric problems
mean(converge_state_m)

# check results
head(pred_list_m[[1]])
length(pred_list_m)
pred_nhanes_fgfpca <- bind_rows(pred_list_m)

# some clean up
pred_nhanes_fgfpca$pred360[pred_nhanes_fgfpca$sind<=360] <- NA
pred_nhanes_fgfpca$pred720[pred_nhanes_fgfpca$sind<=720] <- NA
pred_nhanes_fgfpca$pred1080[pred_nhanes_fgfpca$sind<=1080] <- NA
pred_nhanes_fgfpca <- pred_nhanes_fgfpca %>% select(-pred1440)
head(pred_nhanes_fgfpca)

pred_nhanes_fgfpca %>% 
  filter(id %in% sample(test_id, 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred360, col = "6am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred720, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred1080, col = "6pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

# check AUC

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

auc_appl_fgfpca <- pred_nhanes_fgfpca %>% 
  mutate(window = cut(sind, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>%
  group_by(window) %>% 
  summarise(auc1 = get_auc(Y, pred360),
            auc2 = get_auc(Y, pred720),
            auc3 = get_auc(Y, pred1080)) %>%
  filter(window != "[0,360]") %>% 
  select(starts_with("auc"))

##### save prediction #####
save(pred_nhanes_fgfpca, 
     file = here("Data/ApplOutput_fGFPCA.RData"))

