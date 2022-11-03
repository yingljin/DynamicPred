
set.seed(1025)
##### generate data #####

library(tidyverse)
library(ggplot2)
library(fcr)
library(refund)
library(lme4)

#####  generate data for one subject ####
times <- seq(0, 1, by = 0.001)
gen_data <- function(id = 1){
  scores <- rnorm(4)
  pcs <- cbind(sin(times)^2, cos(times)^2, times^2, times)
  Z <- pcs %*% scores
  probs <- exp(Z)/(1+exp(Z))
  Y <- sapply(probs, rbinom, n = 1, size = 1) 
  
  return(data.frame(id = id, argvals = times, Y=Y, Z=Z))
}

df_train2 <- bind_rows(lapply(1:100, gen_data))
table(df_train2$id)


##### Local GLMMs #####

# bins
w <- 0.01 # bin width
brks <- seq(0, 01, by = w) # cutoff points
mid <- (brks+w/2)[1:100] # mid points
nb <- length(mid) # number of bins

# bin observations
df_train2$time_bin <- cut(df_train2$argvals, breaks = brks, include.lowest = T, labels = mid)
df_train2$time_bin <- as.numeric(as.character(df_train2$time_bin))
df_train2$id <- as.factor(df_train2$id)

table(df_train2$id, df_train2$time_bin)


# fit local linear mixed models in each bin
## df: function with binary outcome Y, with observations in one time bein
pred_latent <- function(df){
  this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial)
  Zhat <- predict(this_glm, type = "link")
  df$Zhat <- Zhat
  return(df)
}


## do that for all time bin
df_bin_lst <- split(df_train2, f = df_train2$time_bin)
df_pred_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
df_pred_latent <- bind_rows(df_pred_latent)



ggplot(df_pred_latent %>% filter(id %in% 15:18))+
  geom_line(aes(x = argvals, y=Zhat))+
  geom_line(aes(x=argvals, y = Z, col = "red"))+
  geom_point(aes(x=argvals, y = Y, col = "blue"), size = 0.01)+
  facet_wrap(~id)





