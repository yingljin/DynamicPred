
set.seed(1025)
##### generate data #####

library(tidyverse)
library(ggplot2)
library(fcr)

# generate data for one subject

GenData<-function(id = 1){
  ## number of observation
  mi <- sample(seq(15, 25, by=1), size = 1)
  ## observation time
  ti <- sample(seq(0, 500, by = 1), size = mi)
  ti <- sort(ti)
  ## coefficient
  ai <- runif(1, -10, 10)
  # function
  yi <- sqrt(ti)+10*sin(ti)+ai*cos(ti)+rnorm(mi, 0, 0.2)
  
  return(data.frame(id = i, t = ti, y=yi))
}

# generate data for N subjects
df <- data.frame()
for(i in 1:100){
  df <- bind_rows(df, GenData(id=i))
}

# check
ggplot(df %>% filter(id %in% 1:5), aes(x=t, y=y, col=as.factor(id)))+
  geom_line(show.legend = F)

##### Fit model #####

# functional domain
allt <- seq(0, 500, by = 1)

# FCR model
fit_fcr <- fcr(formula = y ~ s(t, k=30, bs="ps"),
               argvals = "t", subj = "id", argvals.new = allt,
               data = df, use_bam = T,
               face.args = list(knots = 30, pve = 0.99))

fit_fcr$face.object


##### Out of sample prediction #####

df_new <- data.frame()
for(i in 101:105){
  df_new <- bind_rows(df_new, GenData(id=i))
}

df_new <- subset(df_new, t<=200)
df_new <- rename(df_new, "subj" = "id")


# check
ggplot(df_new, aes(x=t, y=y, col=as.factor(subj)))+
  geom_line(show.legend = F)
pred <- data.frame()

for(i in unique(df_new$subj)){
  this_df <- subset(df_new, subj == i & t <= 300)
  pred_df <- data.frame(subj=i,
                        t = max(this_df$t):500, y=NA)
  pred_df <- bind_rows(this_df, pred_df)
  
  this_pred <- predict(fit_fcr, newdata = pred_df)
  
  pred <- bind_rows(pred,
                    data.frame(subj=i,
                               t = pred_df$t,
                               pred = this_pred$dynamic_predictions$fitted.values$y.pred,
                               se = this_pred$dynamic_predictions$fitted.values$se.fit))
  
  
  
}

subset(pred, subj==101) %>% head()
subset(df_new, subj==101) %>% head()

ggplot(pred)+
  geom_line(aes(x=t, y=pred), show.legend = F)+
  geom_line(aes(x = t, y = pred-se), linetype = "dashed")+
  geom_line(aes(x = t, y = pred+se), linetype = "dashed")+
  geom_point(data = df_new, aes(x=t, y=y))+
  facet_wrap(~subj)
