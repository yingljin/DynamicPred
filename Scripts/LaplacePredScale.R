
set.seed(516)

#### Set up ####

library(here)
library(tidyverse)
library(ggpubr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(kableExtra)
library(knitr)
library(mvtnorm)
library(mgcv)
library(splines)

source(here("Code/GLMM-FPCA.R")) 
source(here("Code/OutsampBayes.R"))

#### generated training data ####

N <- 500 # sample size
J <- 1000 # number of observation points

t = seq(0,1,len=J) # observations points

# mean function
f_0 <- function(s) 0 

#eigenfunctions 
K <- 4 # number of eigenfunctions
phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
                     sin(4*pi*t),cos(4*pi*t))

# eigenvalues
lambda = rep(10, 4)

# score
xi <- matrix(rnorm(N*K),N,K)
xi <- xi %*% diag(sqrt(lambda))

# subject-specific random effect
b_i <- xi %*% t(phi) # of size N by J

# latent gaussian function
eta_i <- t(vapply(1:N, function(x){
  f_0(t) + b_i[x,]
}, numeric(J)))

# outcome binary function
Y_i <- matrix(rbinom(N*J, size=1, prob=plogis(eta_i)), 
              N, J, byrow=FALSE)

train_df <- data.frame(id = factor(rep(1:N, each=J)),
                       t = rep(t, N), 
                       Y = as.vector(t(Y_i)),
                       eta_i = as.vector(t(eta_i)),
                       sind = rep(1:J, N))

# visualization
rand_id <- sample(N, size = 4)
train_df %>% filter(id %in% rand_id) %>% 
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.5)+
  geom_line(aes(x=sind, y=plogis(eta_i)), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)

#### Bin data ####

bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

train_df$bin <- cut(train_df$sind, breaks = brks, include.lowest = T, labels = mid)
train_df$bin <- as.numeric(as.character(train_df$bin))
# unique(df$bin)

train_df %>% 
  filter(id %in% rand_id) %>%
  group_by(id, bin) %>%
  summarise(num = sum(Y)) %>%
  ggplot()+
  geom_point(aes(x=bin, y=num), size = 0.5)+
  facet_wrap(~id)+
  labs(x="Time", y = "Activity", title = "Number of active nimutes within each bin")


#### Local GLMM ####

train_bin_lst <- split(train_df, f = train_df$bin)

t1=Sys.time()
df_est_latent <- lapply(train_bin_lst, function(x){pred_latent(x, n_node = 0)}) 
t2= Sys.time()
t2-t1 

df_est_latent <- bind_rows(df_est_latent) 
# head(df_est_latent)

# example estimated latent function
train_id <- unique(train_df$id)
rand_id <- sample(train_id, 4)

df_est_latent %>% 
  filter(id %in% rand_id) %>%
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  mutate(eta_i = exp(eta_i)/(1+exp(eta_i))) %>%
  ggplot()+
  geom_line(aes(x=t, y=eta_hat, group = id, col = "estimated"))+
  geom_line(aes(x=t, y=eta_i, group = id, col = "true"))+
  geom_point(aes(x=t, y = Y, group = id), size = 0.5)+
  facet_wrap(~id, scales = "free")+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

#### FPCA ####
uni_eta_hat <- df_est_latent %>% filter(bin==sind)
mid_t <- df_est_latent$t[df_est_latent$bin==df_est_latent$sind]
mid_t <- unique(mid_t) # the time points correspoinding to bin midpoints

mat_est_unique <- matrix(uni_eta_hat$eta_hat,
                         nrow=length(train_id), 
                         ncol=n_bin, byrow = F) 

fpca_mod <- fpca.face(mat_est_unique, argvals = mid_t, var=T)

# mean
plot(mid_t, fpca_mod$mu, type = "l", xlab = "bin", ylab = "Mean")

# eigenfunctions
par(mfrow=c(2,2))
plot(mid_t, fpca_mod$efunctions[, 1], type="l", xlab="bin", ylab="PC1")
plot(mid_t, fpca_mod$efunctions[, 2], type="l", xlab="bin", ylab="PC2")
plot(mid_t, fpca_mod$efunctions[, 3], type="l", xlab="bin", ylab="PC3")
plot(mid_t, fpca_mod$efunctions[, 4], type="l", xlab="bin", ylab="PC4")

# eigenvalues
fpca_mod$evalues[1:K]

#### Projection ####

# order of b splines
p <- 3 

# knots
knots <- 35
knots_values <- seq(-p, knots + p, length = knots + 1 + 2 *p)/knots
knots_values <- knots_values * (max(mid_t) - min(mid_t)) + min(mid_t)

# evaluate B-splines on binned grid
B <- spline.des(knots = knots_values, x = mid_t, ord = p + 1,
                outer.ok = TRUE)$design
# evaluate B-splines on original grid
Bnew <- spline.des(knots = knots_values, x = unique(train_df$t), ord = p + 1,
                   outer.ok = TRUE)$design

# project binned eigenfunctions onto the original grid
efunctions_new <- matrix(NA, J, K)
for(k in 1:K){
  lm_mod <- lm(fpca_mod$efunctions[,k] ~ B-1)
  efunctions_new[,k] <- Bnew %*% coef(lm_mod)
}

# visualization
par(mfrow=c(2,2))
plot(mid_t, fpca_mod$efunctions[, 1], xlab="bin", ylab="PC1", pch=20, cex = 0.5)
lines(t, efunctions_new[, 1], col="blue")

plot(mid_t, fpca_mod$efunctions[, 2], xlab="bin", ylab="PC2", pch=20, cex = 0.5)
lines(t, efunctions_new[, 2], col="blue")

plot(mid_t, fpca_mod$efunctions[, 3], xlab="bin", ylab="PC3", pch=20, cex = 0.5)
lines(t, efunctions_new[, 3], col="blue")

plot(mid_t, fpca_mod$efunctions[, 4], xlab="bin", ylab="PC4", pch=20, cex = 0.5)
lines(t, efunctions_new[, 4], col="blue")

#### Debias ####

df_phi <- data.frame(t = t, efunctions_new)
colnames(df_phi) <- c("t", paste0("phi", 1:4))

train_df <- train_df %>% left_join(df_phi, by = "t")
train_df$id <- as.factor(train_df$id)

t1 <- Sys.time()
debias_glmm <- bam(Y ~ s(t, bs="cc", k=10)+
                     s(id, by=phi1, bs="re")+
                     s(id, by=phi2, bs="re")+
                     s(id, by=phi3, bs="re")+
                     s(id, by=phi4, bs="re"), 
                   family = binomial, data=train_df, 
                   method = "fREML",
                   discrete = TRUE)
t2 <- Sys.time()
t2-t1 # less than 45 seconds

# mean function
new_mu <- predict(debias_glmm, type = "terms")

plot(t, new_mu[1:1000, 1], type = "l")

# deviased eigenvalues
new_lambda <- 1/debias_glmm$sp[2:5]
new_lambda

# rescale
scaled_lambda <- new_lambda/n_bin
efunctions_scaled <- efunctions_new*sqrt(n_bin)

#### Generate test data ####
# generated testing data 
N_test <- 200
# xi_test <- matrix(rnorm(N_test*K), N_test, K)
xi_test <- rmvnorm(N_test, mean = rep(0, K), sigma = diag(lambda))
b_test <- xi_test %*% t(phi) # of size N by J
# dim(b_test)

# latent gaussian function
eta_test <- t(vapply(1:N_test, function(x){
  f_0(x) + b_test[x,]
}, numeric(J)))

# outcome binary function
Y_test <- matrix(rbinom(N_test*J, size=1, prob=plogis(eta_test)), 
                 N_test, J, byrow=FALSE)

# format into dataframe
test_df <- data.frame(id = factor(rep(1:N_test, each=J)),
                      t = rep(t, N_test), 
                      Y = as.vector(t(Y_test)),
                      eta_i = as.vector(t(eta_test)),
                      sind = rep(1:J, N_test))

test_df$id <- as.factor(test_df$id)
test_id <- unique(test_df$id)

#### Laplace Approximation ####
Model <- function(parm, Data){
  xi <- parm[Data$pos.xi]
  
  # log-prior
  xi.prior <- dmvnorm(xi, mean = rep(0, Data$K), sigma=Data$tao, log = TRUE)
  
  # log-posterior likelihood
  eta <- Data$f0+Data$X %*% xi
  p <- exp(eta)/(1+exp(eta))
  LL <- sum(dbern(x=Data$y, prob=p, log = TRUE)) # log likelihood of Y|xi
  LP <- LL+sum(xi.prior) # joint log likelihood of (Y, xi)
  
  # output
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LL,
                   yhat=Data$y, parm=parm)
  return(Modelout)
}

# parameter names
name_lst <- as.list(rep(0, K))
names(name_lst) <- paste("xi", 1:K, sep = "")
parm.names <- as.parm.names(name_lst)
pos.xi <- grep("xi", parm.names)
mon.names <- "LP"

#### For one single subject ####
this_id <- 43
this_t <- 0.4

# data
df_it <- test_df %>% filter(id==this_id & t <= this_t)
max_rid <- nrow(df_it)

# use Unscaled params 
MyData_unscaled <- list(K=K, 
                        X=efunctions_new[1:max_rid, ], 
                        mon.names=mon.names,
                        parm.names=parm.names, 
                        pos.xi=pos.xi, 
                        y=df_it$Y, 
                        tao=diag(new_lambda), f0=new_mu[1:max_rid])
Fit_unscaled <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData_unscaled,
                                     Method = "BFGS", Iterations = 1000,
                                     CovEst = "Identity")
Fit_unscaled$Converged

# Scaled params
MyData_scaled <- list(K=K, 
                      X=efunctions_scaled[1:max_rid, ], 
                      mon.names=mon.names,
                      parm.names=parm.names, 
                      pos.xi=pos.xi, 
                      y=df_it$Y, 
                      tao=diag(scaled_lambda), f0=new_mu[1:max_rid])
Fit_scaled <- LaplaceApproximation(Model, parm = rep(0, K), 
                                   Data=MyData_scaled, Method = "BFGS", 
                                   Iterations = 1000,
                                   CovEst = "Identity")

## iterative distribution, compared to true value
true_xi_it <- data.frame(name = colnames(Fit_unscaled$History),
                         true_xi = xi_test[id_diff$id[1], ])
Fit_scaled$Converged

# compare 
Fit_unscaled$Summary1
Fit_scaled$Summary1
xi_test[43, ]


#### For all subjects ####
new_mu <- new_mu[1:1000]
# new_mu[1:6, 1]
# new_mu[1001:1006, 1]


# prediction using unscaled parameter
# container
pred_list <- list()
converge_vec <- matrix(NA, N_test, 4)

t1 <- Sys.time()
# per subject
# seq_along(test_id)
for(i in seq_along(test_id)){
  df_i <- test_df %>% filter(id==test_id[i])
  
  # per max obs time
  for(tmax in c(0.2, 0.4, 0.6, 0.8)){
    df_it <- df_i %>% filter(t <= tmax)
    max_rid <- nrow(df_it)
    
    # into a list
    MyData <- list(K=K, 
                   X=efunctions_new[1:max_rid, ], 
                   mon.names=mon.names,
                   parm.names=parm.names, 
                   pos.xi=pos.xi, 
                   y=df_it$Y, 
                   tao=diag(new_lambda), f0=new_mu[1:max_rid])
    
    
    # fit laplace approximation
    Fit <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData, Method = "BFGS", Iterations = 1000,
                                CovEst = "Identity")
    converge_vec[i, which(c(0.2, 0.4, 0.6, 0.8)==tmax)] <- Fit$Converged
    # Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values = parm = rep(0, K),  Method = "BFGS")
    score <- Fit$Summary1[, "Mode"]
    
    # prediction
    eta_pred_out <- new_mu+efunctions_new%*%score
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
  }
  
  # df_i$pred0.2[df_i$t<=0.2] <- NA
  # df_i$pred0.4[df_i$t<=0.4] <- NA
  # df_i$pred0.6[df_i$t<=0.6] <- NA
  # df_i$pred0.8[df_i$t<=0.8] <- NA
  
  pred_list[[i]] <- df_i
}
t2 <- Sys.time()
t_pred <- t2-t1 # About 3.5 minutes

mean(converge_vec)

# scaled parameters
pred_list3 <- list()
converge_vec2 <- matrix(NA, N_test, 4)

# per subject
for(i in seq_along(test_id)){
  df_i <- test_df %>% filter(id==test_id[i])
  
  # per max obs time
  for(tmax in c(0.2, 0.4, 0.6, 0.8)){
    df_it <- df_i %>% filter(t <= tmax)
    max_rid <- nrow(df_it)
    
    # into a list
    MyData <- list(K=K,
                   X=efunctions_scaled[1: max_rid, ], 
                   mon.names=mon.names,
                   parm.names=parm.names, 
                   pos.xi=pos.xi, 
                   y=df_it$Y, 
                   tao=diag(scaled_lambda), f0=new_mu[1:max_rid])
    
    
    # fit laplace approximation
    Fit <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData, Method = "BFGS", Iterations = 1000,
                                CovEst = "Identity")
    converge_vec2[i, which(c(0.2, 0.4, 0.6, 0.8)==tmax)] <- Fit$Converged
    score <- Fit$Summary1[, "Mode"]
    
    # prediction
    eta_pred_out <- new_mu+efunctions_scaled%*%score
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
  }
  
  # df_i$pred0.2[df_i$t<=0.2] <- NA
  # df_i$pred0.4[df_i$t<=0.4] <- NA
  # df_i$pred0.6[df_i$t<=0.6] <- NA
  # df_i$pred0.8[df_i$t<=0.8] <- NA
  
  pred_list3[[i]] <- df_i
}


mean(converge_vec2)

## compare
df_pred <- bind_rows(pred_list)

df_pred$pred0.2[df_pred$t<=0.2] <- NA
df_pred$pred0.4[df_pred$t<=0.4] <- NA
df_pred$pred0.6[df_pred$t<=0.6] <- NA
df_pred$pred0.8[df_pred$t<=0.8] <- NA

rand_test_id <- sample(test_id, 4)

p1<- df_pred %>%
  filter(id %in% rand_test_id) %>%
  mutate_at(vars(eta_i, pred0.2, pred0.4, pred0.6, pred0.8), function(x){exp(x)/(1+exp(x))}) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), linetype="dashed", na.rm=T)+
  facet_wrap(~id)+
  labs(title = "Unscaled, debiased eigenvalues")

# scaled
df_pred3 <- bind_rows(pred_list3)

df_pred3$pred0.2[df_pred3$t<=0.2] <- NA
df_pred3$pred0.4[df_pred3$t<=0.4] <- NA
df_pred3$pred0.6[df_pred3$t<=0.6] <- NA
df_pred3$pred0.8[df_pred3$t<=0.8] <- NA

p3 <- df_pred3 %>%
  filter(id %in% rand_test_id) %>%
  mutate_at(vars(eta_i, pred0.2, pred0.4, pred0.6, pred0.8), function(x){exp(x)/(1+exp(x))}) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), linetype="dashed", na.rm=T)+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), linetype="dashed", na.rm=T)+
  facet_wrap(~id)+
  labs(title = "Scaled, unbiased eigenvalues")

# calcualte ISE
err1 <- df_pred %>%
  mutate(err1 = (pred0.2-eta_i)^2,
         err2 = (pred0.4-eta_i)^2,
         err3 = (pred0.6-eta_i)^2,
         err4 = (pred0.8-eta_i)^2) %>%
  select(id, t, starts_with("err"))

err1$window = cut(err1$t, breaks = seq(0, 1, by = 0.2), include.lowest = T)
# table(err1$window)
err1 <- split(err1, f = err1$window)

tb1 <- lapply(err1, function(x){
  x %>%  group_by(id) %>% 
    summarise_at(vars(err1, err2, err3, err4), sum) %>%
    summarise_at(vars(err1, err2, err3, err4), mean)
}) %>% bind_rows(.id = "Window")

# calcualte ISE
err3 <- df_pred3 %>%
  mutate(err1 = (pred0.2-eta_i)^2,
         err2 = (pred0.4-eta_i)^2,
         err3 = (pred0.6-eta_i)^2,
         err4 = (pred0.8-eta_i)^2) %>%
  select(id, t, starts_with("err"))

err3$window = cut(err3$t, breaks = seq(0, 1, by = 0.2), include.lowest = T)
# table(err1$window)
err3 <- split(err3, f = err3$window)

tb3 <- lapply(err3, function(x){
  x %>%  group_by(id) %>% 
    summarise_at(vars(err1, err2, err3, err4), sum) %>%
    summarise_at(vars(err1, err2, err3, err4), mean)
}) %>% bind_rows(.id = "Window")

tb1

tb3
