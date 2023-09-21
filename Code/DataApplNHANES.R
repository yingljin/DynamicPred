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
# ggsave(here("Images/NHANES/data.png"), width=7, height=5, bg="white")



##### Bin data #####

# bin data
bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

# bin labels
df$bin <- cut(df$sind, breaks = brks, include.lowest = T, labels = mid)
df$bin <- as.numeric(as.character(df$bin))

df <- df %>% rename(id = SEQN, Y=Z)


##### Data split #####

id_vec <- unique(df$id)
train_id <- sample(id_vec, size = 0.8*N, replace = FALSE) # 7010 subjects for training
test_id <- setdiff(id_vec, train_id) # 1753 subjects for testing

#### fGFPCA estimation ####
# fit model on the training set
df_train <- df %>% filter(id %in% train_id)
train_bin_lst <- split(df_train, f = df_train$bin)
# length(train_bin_lst)
# lapply(train_bin_lst, nrow)
# train_bin_lst[[1]] %>% View()

# local GLMM and estimate latent function
# near-unidentifiability issues 
t1=Sys.time()
df_est_latent <- lapply(train_bin_lst, function(x){pred_latent(x, n_node = 0)}) 
t2= Sys.time()
t_local_glmm <- t2-t1 # 3.3 minutes for model estimation
# no numeric warnings appeared

lapply(df_est_latent, function(x)range(x$eta_hat))
lapply(df_est_latent, function(x)unique(x$eta_hat)) %>% head()
df_est_latent[[1]] %>% View()

# on the binned grid, one unique value each bin
df_est_latent <- bind_rows(df_est_latent) %>% 
  select(-sind, -Y) %>% distinct(.) 

# visualization of the estimated latent function
rand_id <- sample(train_id, size = 4)
df %>% 
  filter(id %in% rand_id) %>%
  left_join(df_est_latent, by = c("id", "bin")) %>%
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_line(aes(x=sind, y=eta_hat, group = id))+
  geom_point(aes(x=sind, y = Y, group = id), size = 0.5)+
  facet_wrap(~id, scales = "free")+
  labs(x = "Time", y = "Estimated latent function (probablity scale)")

# fPCA
mat_est_unique <- matrix(df_est_latent$eta_hat,
                         nrow=length(train_id), 
                         ncol=n_bin, byrow = F) # row index subject, column binned time
dim(mat_est_unique)
fpca_mod <- fpca.face(mat_est_unique, argvals = mid, var=T)

dim(fpca_mod$efunctions) # 27 eigenfunctions total
fpca_mod$evalues
K <- 4 # use for eigen functions for prediction

#### Dynamic prediction ####
# Out-of-sample prediction
df_test <- df %>% filter(id %in% test_id)
df_test <- df_test %>% select(-sind, -Y) %>% distinct(.)

# observations time: 9am (540), 1pm (780),  5pm (1020)
df_test[, 'pred_t540'] <- df_test[, 'pred_t780'] <- df_test[, 'pred_t1020'] <- NA
# score_out_mat <- array(NA, dim = c(length(test_id), K, 3))

# predict all test subjects
skip_id <- c()

t1=Sys.time()
pb = txtProgressBar(min = 0, max = length(test_id), initial = 0, style = 3) 
# prediction for a single subject
for(i in seq_along(test_id)){
  df_i <- df %>% filter(id==test_id[i]) 
  
  # prediction 
  ## up to 540
  pred_t1 <- tryCatch({out_pred_laplace(fpca_mod, df_i %>% filter(sind<=540), kpc = K)},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t1))==0){
    df_test[df_test$id == test_id[i], 'pred_t540'] <- pred_t1$eta_pred
    # score_out_mat[i, ,1] <- pred_t1$score_out
  } else{
    skip_id <- append(skip_id, test_id[i])
  }
  
  ## up to 780
  pred_t2 <- tryCatch({out_pred_laplace(fpca_mod, df_i %>% filter(sind<=780), kpc = K)},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t2))==0){
    df_test[df_test$id == test_id[i], 'pred_t780'] <- pred_t2$eta_pred
    # score_out_mat[i, ,1] <- pred_t1$score_out
  } else{
    skip_id <- append(skip_id, test_id[i])
  }
  
  ## up to 1020
  pred_t3 <- tryCatch({out_pred_laplace(fpca_mod, df_i %>% filter(sind<=1020), kpc = K)},
                      error = function(e){return(NA)})
  if(sum(is.na(pred_t3))==0){
    df_test[df_test$id == test_id[i], 'pred_t1020'] <- pred_t3$eta_pred
    # score_out_mat[i, ,1] <- pred_t1$score_out
  } else{
    skip_id <- append(skip_id, test_id[i])
  }
  
  setTxtProgressBar(pb, i)
}
t2=Sys.time()

t_pred = t2-t1 # total 20min for out-of-sample prediction
t_pred/length(test_id) # about 0.01 mins per subject 


#### numeric problems ####
skip_id # 65876, 80414, 83271 had numeric problems

df_test[df_test$id==65876, ] %>% View() 
df_test[df_test$id==83271, ] %>% View() 
df_test[df_test$id==80414, ] %>% View() 

# for all three subjects, failed approximation happened with the shortest observer track (540)

df %>% filter(id %in% skip_id) %>% 
  ggplot()+
  geom_point(aes(x=sind, y=Y))+
  facet_wrap(~id)+
  geom_vline(xintercept = c(540, 780, 1020))

#### check results ####
# score
# score_out_mat[1,,]

# prediction on probability scale
df_test$pred_t540[df_test$bin<=540] <- NA
df_test$pred_t780[df_test$bin<=780] <- NA
df_test$pred_t1020[df_test$bin<=1020] <- NA

rand_id <- sample(test_id, 4)

df %>% 
  filter(id %in% rand_id) %>%
  left_join(df_test %>% select(id, bin, pred_t540, pred_t780, pred_t1020)) %>%
  mutate_at(vars(pred_t540, pred_t780, pred_t1020), function(x)exp(x)/(1+exp(x))) %>% 
  ggplot()+
    geom_line(aes(x=bin, y = pred_t540, col = "9am"), linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t780, col = "1pm"),linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t1020, col = "5pm"),linetype = "dashed")+
    geom_point(aes(x=bin, y = Y, col = "Outcome"), size = 0.2)+
    facet_wrap(~id)+
    labs(x = "Time", y = "Estimated latent function (probablity scale)")
# ggsave(here("Images/NHANES/pred_eta_prob.png"), width=7, height=5, bg="white")

# prediction on latent function scale
df_test %>%
  filter(id %in% rand_id) %>% 
  ggplot()+
    geom_line(aes(x=bin, y = pred_t540, col = "9am"), linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t780, col = "1pm"),linetype = "dashed")+
    geom_line(aes(x=bin, y = pred_t1020, col = "5pm"),linetype = "dashed")+
    # geom_line(aes(x=bin, y = eta_hat, col = "eta_hat"), linetype = "solid" )+
    facet_wrap(~id)+
    labs(x = "Time", y = "Estimated latent function")
# ggsave(here("Images/NHANES/pred_eta_latent.png"), width=7, height=5, bg="white")


# save results
save(df_test, skip_id, file = here("Data/ApplOutput_fGFPCA.RData"))


#### Reference method: GLMMadaptive ####

# fit GLMMadaptvie model on the training set
head(df_train)
t1 <- Sys.time()
adglmm_mod <- mixed_model(Y ~ sind, random = ~ 1 | id, data = df_train %>% mutate(sind=sind/J),
                      family = binomial())
# adglmm_mod <- mixed_model(Y ~ sind, random = ~ sind | id, data = df_train %>% mutate(sind=sind/J),
#                           family = binomial())
t2 <- Sys.time()
t_est_adglmm <- t2-t1 # model fitting took 20.72 mins
summary(adglmm_mod)

# use the fitted model for prediction
df_test <- df %>% filter(id %in% test_id)
df_test$sind <- df_test$sind/J
head(df_test)

## up to 540
t1 <- Sys.time()
adglmm_pred_t540 <- predict(adglmm_mod, 
        newdata = df_test %>% filter(sind <= 540/J),
        newdata2 =  df_test %>% filter(sind > 540/J), 
        type = "subject_specific", type_pred = "link",
        se.fit = TRUE, return_newdata = TRUE)

adglmm_pred_t780 <- predict(adglmm_mod, 
                            newdata = df_test %>% filter(sind <= 780/J),
                            newdata2 =  df_test %>% filter(sind > 780/J), 
                            type = "subject_specific", type_pred = "link",
                            se.fit = TRUE, return_newdata = TRUE)

adglmm_pred_t1020 <- predict(adglmm_mod, 
                            newdata = df_test %>% filter(sind <= 1020/J),
                            newdata2 =  df_test %>% filter(sind > 1020/J), 
                            type = "subject_specific", type_pred = "link",
                            se.fit = TRUE, return_newdata = TRUE)
t2 <- Sys.time()
t_pred_adglmm <- t2-t1 # 12.18 minutes spent on prediction

save(adglmm_pred_t540, adglmm_pred_t780, adglmm_pred_t1020, file = here("Data/ApplOutput_GLMMadaptive.RData"))
