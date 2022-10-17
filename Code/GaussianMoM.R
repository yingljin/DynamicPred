library(fcr)
library(ggplot2)
library(tidyverse)

# https://cran.r-project.org/web/packages/fcr/vignettes/dynamic-prediction.html


#### Data ####
df = content
df_train = subset(content, include == 1)
df_test = subset(content, include == 0)

ggplot(df)+
  geom_line(aes(x = argvals, y = Y, col = as.factor(subj)), show.legend = F)



#### functional concurrent regression ####

k = 12

# functional domain
tnew = unique(df$argvals)
tnew = sort(tnew)

# smoothing time-varying covariate: measured with error
ggplot(df_train)+
  geom_line(aes(x = argvals, y = waz.true, col = as.factor(subj)), show.legend = F)

df_waz <- data.frame("y" = df_train$waz, "subj" = df_train$subj, "argvals" = df_train$argvals)
fit_waz <- face.sparse(df_waz, newdata = df_waz, knots = k, argvals.new = tnew) # not 
## mean function
plot(tnew, fit_waz$mu.new)
## covariance function
heatmap(fit_waz$Cor.new)
heatmap(fit_waz$Cor.raw.new)
## eigenfunctions
fit_waz$eigenvalues
plot(tnew, fit_waz$eigenfunctions[, 1])
plot(tnew, fit_waz$eigenfunctions[, 2])
plot(tnew, fit_waz$eigenfunctions[, 3])
plot(tnew, fit_waz$eigenfunctions[, 4])
## smoothed value 
df_train$sm_waz <- fit_waz$y.pred


# fit FCR 
K2 <- 15

fit_fcr <- fcr(formula = Y ~ s(argvals, k=K2, bs="ps") + s(argvals, by=Male, k=K2, bs="ps") + s(argvals, by=sm_waz, k=K2, bs="ps"), 
               argvals = "argvals", subj = "subj", data = df_train, use_bam = T, argvals.new = tnew,
               face.args = list(knots = K2, pve = 0.99))
## plot coefficients
par(mfrow = c(1, 1))
plot(fit_fcr,select=3)
## plot covariance function
plot(fit_fcr, plot.covariance=TRUE) # eigen functions different from tutorial

#### Prediction #####
# in-sample: for subject used in the model fitting
# dynamic: for subjects not used in the model fitting
# difference in whether subject ID is in the training set

# Example: subject 1
data_pred <- subset(df_train, subj == 1) # insample
# data_pred_dyn: dynamic prediction
# data_pred_dyn2: dynamic prediction but with all data
data_pred_dyn  <- data_pred_dyn2 <- data_pred
data_pred_dyn$subj <- data_pred_dyn2$subj <- data_pred$subj + 100000 # new subject ID

# set time-varying covariate and outcome beyond t as NA
data_pred_dyn$Y[data_pred_dyn$argvals >= 0.5] <- NA
data_pred_dyn$waz[data_pred_dyn$argvals >= 0.5] <- NA

# use fpca to predict time-varying covariate
data_dyn_waz <- data.frame(y = data_pred_dyn$waz, subj = data_pred_dyn$subj, argvals = data_pred_dyn$argvals)
data_pred_dyn$sm_waz <- predict(fit_waz, newdata = data_dyn_waz)$y.pred

data_pred_dyn$sm_waz[data_pred_dyn$argvals >= 0.5]
data_pred_dyn2$sm_waz[data_pred_dyn2$argvals >= 0.5]

# prediction
in_sample <- predict(fit_fcr, newdata=data_pred)
dyn <- predict(fit_fcr, newdata = data_pred_dyn)
dyn2 <- predict(fit_fcr, newdata = data_pred_dyn2)
dyn$dynamic_predictions$fitted.values
dyn2$dynamic_predictions$fitted.values

plot(data_pred$argvals, in_sample$insample_predictions, type = "l")
lines(data_pred$argvals, dyn$dynamic_predictions$fitted.values$y.pred, col = 2)
lines(data_pred$argvals, dyn2$dynamic_predictions$fitted.values$y.pred, col = 3)
# dynamic prediction with all data
# does not need to predict time-varying covariates


#### Prediction on the test set ####

# on the last 1/2 of the functional domain
data_dyn <- df_test
data_dyn$Y[data_dyn$argvals >= 0.5]   <- NA
data_dyn$waz[data_dyn$argvals >= 0.5] <- NA

# remove smoothed waz 
# -- dyanmic predictions calculated shortly
data_dyn$wazPred <- NULL


set.seed(1012341)
ids <- sample(unique(data_dyn$subj), 4, replace = FALSE)


data_plot <- c()
for(i in 1:length(ids)){
  tmp <- subset(data_dyn, subj %in% ids[i])
  ut_tmp <- tnew[!tnew %in% tmp$argvals & tnew > min(tmp$argvals)]
  n_ut <- length(ut_tmp)
  empty_dat <- data.frame("Y" = rep(NA,n_ut),
                          "Ytrue" = rep(NA,n_ut),
                          "X" = rep(NA,n_ut),
                          "waz.true"=rep(NA,n_ut),
                          "waz"=rep(NA,n_ut),
                          "Male" = rep(tmp$Male[1], n_ut),
                          "argvals"=ut_tmp,
                          "subj"=rep(tmp$subj[1],n_ut),
                          "include"=rep(0,n_ut))
  tmp <- rbind(tmp, empty_dat) # why combine the two data? why not just use tnew? 
  tmp <- tmp[order(tmp$argvals),]
  
  data_plot <- rbind(data_plot, tmp)
  rm(list=c("tmp","ut_tmp","n_ut","empty_dat"))
}

data_dyn <- data_plot
rm(list=c("data_plot"))
data_dyn[data_dyn$subj == 177,] %>% View()

# in fact, observations from one subject can be spread out
# across the time grid

# get dynamic waz predictions (from face.sparse)
data_dyn_waz <- data.frame("y" = data_dyn$waz, 
                           "subj" = data_dyn$subj, 
                           "argvals" = data_dyn$argvals)
data_dyn$sm_waz <- predict(fit_waz, newdata=data_dyn_waz)$y.pred
ggplot(data_dyn)+
  geom_line(aes(x=argvals, y=sm_waz, col=as.factor(subj)))

# predict the outcome
preds <- predict(fit_fcr, newdata=data_dyn)$dynamic_predictions
preds$fitted.values

# visualization
par(mfrow=c(2,2),las=1)
for(i in 1:length(ids)){
  inx <- which(data_dyn$subj==ids[i])
  inx2 <-which(df_test$subj==ids[i])
  
  yl  <- range(c(preds$fitted.values$y.pred[inx] + 
                   rep(c(1,-1),each=length(inx))*1.96*preds$fitted.values$se.fit.p[inx],
                 df_test$Y[inx2]))
  
  plot(preds$data$argvals[inx], preds$fitted.values$y.pred[inx],type='l',
       main = paste("Subject",i),xlab="Functional Domain", ylab="Response",
       xlim=c(0,1), ylim=yl)
  lines(preds$data$argvals[inx], lty=2,
        preds$fitted.values$y.pred[inx] - 1.96*preds$fitted.values$se.fit[inx],col='blue')
  lines(preds$data$argvals[inx], lty=2,
        preds$fitted.values$y.pred[inx] + 1.96*preds$fitted.values$se.fit[inx],col='blue')
  
  lines(preds$data$argvals[inx],  lty=2,
        preds$fitted.values$y.pred[inx] - 1.96*preds$fitted.values$se.fit.p[inx],col='red')
  lines(preds$data$argvals[inx],  lty=2,
        preds$fitted.values$y.pred[inx] + 1.96*preds$fitted.values$se.fit.p[inx],col='red')
  axis(1,0.5,labels=expression(t[m]))
  
  points(data_test$argvals[inx2],data_test$Y[inx2])
  
  abline(v=0.5,col='grey',lty=2)
  if(i == 4){
    legend("top",c("Observed Data","fcr() Prediction","95% Confidence Interval","95% Prediction Interval"),
           col=c("black","black","blue","red"),bty="n",pch=c(1,NA,NA,NA), lty = c(NA,1,2,2))
  }
}
