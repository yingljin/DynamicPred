# data application of NHANES

set.seed(121)

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

#### load data ####
df <- read_rds(here("Data/nhanes_bi.rds"))
df <- df %>% rename(id=SEQN, Y=Z)

N <- length(unique(df$id)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject
t <- unique(df$sind)
K <- 4 # number of eigenfunctions to use


#### visulization ####
rand_id <- sample(unique(df$id), size = 4)

df %>% 
  filter(id %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.5)+
  facet_wrap(~id)+
  labs(x="Time", y = "Activity")+
  theme(plot.margin = margin(0, 1, 0, 0))


#### Data split ####
# 80% (7010) subjects for training, 20% (1753) for out-of-sample prediction
train_id <- sample(unique(df$id), size = N*0.8)
test_id <- setdiff(unique(df$id), train_id)

train_df <- df %>% filter(id %in% train_id)
test_df <- df %>% filter(id %in% test_id)


##### fGFPCA #####

# functions
source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutSampBayes.R"))

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points
mid_t <- t[mid] # time value corresponding to mid points

train_df$bin <- cut(train_df$sind, breaks = brks, include.lowest = T, labels = mid)
train_df$bin <- as.numeric(as.character(train_df$bin))

# local GLMM
t1 <- Sys.time()
bin_lst <- split(train_df, f = train_df$bin)
df_est_latent <- lapply(bin_lst, function(x){pred_latent(x, n_node = 0)}) 
df_est_latent <- bind_rows(df_est_latent) 
head(df_est_latent)

# FPCA
uni_eta_hat <- df_est_latent %>% filter(bin==sind)
head(uni_eta_hat)
mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                         nrow=length(train_id), 
                         ncol=n_bin, byrow = F) 
fpca_mod <- fpca.face(mat_est_unique, argvals = mid_t, var=T)
ncol(fpca_mod$efunctions)
fpca_mod$evalues[1:4]

# Re-evaluation
## grid extension
p <- 3 # order of b splines 
knots <- 35 # number of knots (same from FPCA model)
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)

B <- spline.des(knots = knots_values, x = mid_t, ord = p + 1,
                outer.ok = TRUE)$design  # evaluate B-splines on binned grid
Bnew <- spline.des(knots = knots_values, x = t, ord = p + 1,
                   outer.ok = TRUE)$design  # evaluate B-splines on original grid

df_phi <- matrix(NA, J, K) 
for(k in 1:K){
  lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
  df_phi[,k] <- Bnew %*% coef(lm_mod)
}# project binned eigenfunctions onto the original grid


plot(mid_t, fpca_mod$efunctions[,2], pch = ".")
lines(t, df_phi[,2])

## debias
df_phi <- data.frame(t = t, df_phi)
colnames(df_phi) <- c("t", paste0("phi", 1:4))
train_df <- train_df %>% rename(t=sind) %>%
  left_join(df_phi, by = "t")
train_df$id <- as.factor(train_df$id)
debias_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                     s(id, by=phi1, bs="re")+
                     s(id, by=phi2, bs="re")+
                     s(id, by=phi3, bs="re")+
                     s(id, by=phi4, bs="re"), 
                   family = binomial, data=train_df, 
                   method = "fREML",
                   discrete = TRUE)
new_mu <- predict(debias_glmm, type = "terms")[1:J, 1] # extract re-evaluated mean
plot(t, new_mu)

new_lambda <- 1/debias_glmm$sp[2:5] # extract re-evaluated lambda

## rescale
new_phi <- df_phi %>% select(starts_with("phi"))*sqrt(n_bin)
new_phi <- as.matrix(new_phi)

new_lambda <- new_lambda/n_bin
t2 <- Sys.time()
t2-t1


# prediction
t1 <- Sys.time()
# per subject
test_id
window <- seq(0, J, by = 240)

converge_state_m <- matrix(NA, nrow = length(test_id), 5)
pred_list_m <- list()

for(i in seq_along(test_id)){
  df_i <- test_df %>% filter(id==test_id[i])
  
  # per max obs time
  for(tmax in window[2:6]){
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
t_pred <- t2-t1 # About 3.5 minutes
pred_time[m] <- t2-t1


# numeric problems
mean(converge_state_m)

# check results
df_pred <- bind_rows(pred_list_m)
df_pred$pred240[df_pred$sind<=240] <- NA
df_pred$pred480[df_pred$sind<=480] <- NA
df_pred$pred720[df_pred$sind<=720] <- NA
df_pred$pred960[df_pred$sind<=960] <- NA
df_pred$pred1200[df_pred$sind<=1200] <- NA


df_pred %>% 
  filter(id %in% sample(test_id, 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred240, col = "4am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred480, col = "8am"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred720, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred960, col = "4pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred1200, col = "8pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

nhanes_pred_fgfpca <- df_pred

#### Reference method: GLMMadaptive ####

# fit GLMMadaptvie model on the training set
head(train_df)

t1 <- Sys.time()
adglmm_mod <- mixed_model(Y ~ t, random = ~ 1 | id, 
                          data = train_df %>% mutate(t=t/J),
                      family = binomial())
t2 <- Sys.time()
t_est_adglmm <- t2-t1 # model fitting took 20.82 mins
summary(adglmm_mod)

# prediction
test_df <- test_df %>% rename(t=sind) %>%
  mutate(t=t/J)

window/J

pred_list_m <- list(
  predict(adglmm_mod, newdata = test_df %>% filter(t <= 240/J),
          newdata2 =  test_df %>% filter(t>240/J), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(adglmm_mod, newdata = test_df %>% filter(t <= 480/J),
          newdata2 =  test_df %>% filter(t>480/J), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(adglmm_mod, newdata = test_df %>% filter(t <= 720/J),
          newdata2 =  test_df %>% filter(t>720/J), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(adglmm_mod, newdata = test_df %>% filter(t <= 960/J),
          newdata2 =  test_df %>% filter(t>960/J), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2,
  
  predict(adglmm_mod, newdata = test_df %>% filter(t <= 1200/J),
          newdata2 =  test_df %>% filter(t>1200/J), 
          type = "subject_specific", type_pred = "link",
          se.fit = FALSE, return_newdata = TRUE)$newdata2
)

# check results
pred_list_m <- lapply(pred_list_m, function(x){x %>% select(id, t, pred)})
lapply(pred_list_m, head)
df_pred <- test_df %>% 
  left_join(pred_list_m[[1]]) %>% 
  rename(pred240 = pred) %>% 
  left_join(pred_list_m[[2]]) %>% 
  rename(pred480 = pred) %>% 
  left_join(pred_list_m[[3]]) %>% 
  rename(pred720 = pred) %>% 
  left_join(pred_list_m[[4]]) %>% 
  rename(pred960 = pred) %>% 
  left_join(pred_list_m[[5]]) %>% 
  rename(pred1200 = pred)

head(df_pred)

df_pred %>% filter(t>(480/J)) %>% head()

df_pred %>% 
  rename(sind=t) %>%
  filter(id %in% sample(test_id, 4)) %>%
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred240, col = "4am"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred480, col = "8am"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred720, col = "12pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred960, col = "4pm"),linetype = "dashed")+
  geom_line(aes(x=sind, y = pred1200, col = "8pm"),linetype = "dashed")+
  geom_point(aes(x=sind, y = Y, col = "Outcome"), size = 0.2)+
  facet_wrap(~id)+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

nhanes_pred_adglmm <- df_pred



save(nhanes_pred_fgfpca, nhanes_pred_adglmm,
     file = here("Data/ApplOutput.RData"))
