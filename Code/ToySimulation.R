
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

## four subjects
ggplot(df_train2 %>% filter(id %in% 1:4))+
  geom_point(aes(x=argvals, y=Y, col = "blue"),  size = 0.1,
             show.legend = F)+
  geom_line(aes(x=argvals, y=Z, col = "Red"), 
            show.legend = F)+
  facet_wrap(~id)+
  labs(title = "Simulated function for the 4 subjects",
       x = "time", y ="")

## all subjects

ggplot(df_train2)+
  geom_line(aes(x=argvals, y=Z, col = as.factor(id)), 
            show.legend = F)


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

# fit local linear mixed models in each bin
## df: function with binary outcome Y, id and time_bin
pred_latent <- function(df = df_train2, t = 0.005){
  this_df <- df[df$time_bin==t, ]
  this_glm <- glmer(Y ~ 1 + (1|id), data = this_df, family = binomial)
  coef(this_glm)
  Zhat <- predict(this_glm, type = "link")
  return(Zhat)
}

df_try <- df_train2 %>% select(id, Z, argvals, time_bin) %>% 
  filter(time_bin == 0.025) %>% 
  mutate(Zhat = pred_latent(t = 0.025))

ggplot(df_try %>% filter(id %in% 1:8))+
  geom_line(aes(x = argvals, y=Z, col = "red"))+
  geom_line(aes(x = argvals, y=Zhat, cold = "blue"))+
  facet_wrap(~id, )

allZhat <- unlist(sapply(df_train2_lst, pred_latent))
df_train_sm <- data.frame(subj = df_train2$id, 
                          argvals = df_train2$time_bin, 
                          y = allZhat) %>% 
  distinct()


ggplot(df_train_sm %>% filter(subj %in% 1:4))+
  geom_line(aes(x=argvals, y=y), 
            show.legend = F)+
  facet_wrap(~subj)+
  labs(title = "Estimated latent function tracks for the 4 subjects",
       x = "time", y ="")


rm(bin, bins, t, times, tnew)




