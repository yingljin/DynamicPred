##### Local GLMMs #####

bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

# bin observations
df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
df$bin <- as.numeric(as.character(df$bin))

head(df, 20)

# check overlap
# df %>% 
#   filter(id %in% 1:4) %>% 
#   ggplot()+
#   geom_line(aes(x=sind, y=eta_i, col=as.factor(bin)), show.legend = F)+
#   geom_point(aes(x=sind, y=Y, col=as.factor(bin)), size = 0.05, show.legend = F)+
#   facet_wrap(~id)

# fit local linear mixed models in each bin
## df: function with binary outcome Y, with observations in one time bein
pred_latent <- function(df){
  this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial)
  eta_hat <- predict(this_glm, type = "link")
  df$eta_hat <- eta_hat
  return(df)
}

## do that for all time bin
df_bin_lst <- split(df, f = df$bin)
df_pred_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
lapply(df_pred_latent, dim)
df_pred_latent <- bind_rows(df_pred_latent)

#### fPCA on predicted latent variables #####
rm(df_bin_lst)

# put predicted latent function to wide format
df_pred_latent <- df_pred_latent %>% arrange(id, sind_inx)
df_pred_unique <- df_pred_latent %>% select(id, bin, eta_hat) %>% distinct()

mat_pred_unique <- matrix(df_pred_unique$eta_hat, nrow=N, ncol=n_bin, byrow = T)
# dim(mat_pred_unique)

# fPCA
fpca_fit <- fpca.face(mat_pred_unique, pve = 0.95, argvals = unique(df_pred_unique$bin), knots=20, var=T)
# fpca_fit <- face.sparse(df_pred_unique, argvals.new = 1:J, knots = 20, pve=0.95)

# Things for prediction
tao <- diag(fpca_fit$evalues) ## egienvalues, also var-cov matrix of scores
allB <- fpca_fit$efunctions ## eigen functions on the entire binned grid
fpca_fit$mu ## mean function
fpca_fit$sigma2 ## residual variance



##### In-sample prediction on binned grid #####

# start with one obesrvation
df_id1 <- df_pred_latent %>% filter(id==1 & sind_inx <= 300) %>% 
  select(id, bin, eta_hat) %>% distinct(.)

obs_bin <- length(unique(df_id1$bin))

B <- allB[1:obs_bin, ] # eigenfunctions up to observed time point

# score
score_hat <- tao %*% t(B) %*% solve(B%*%tao%*%t(B)+fpca_fit$sigma2*diag(obs_bin)) %*% (df_id1$eta_hat-fpca_fit$mu[1:obs_bin])

fpca_fit$scores[1, ] ## scores from fPCA face

# prediction
df_pred_id1 <- df_pred_latent %>% filter(id==1 )

eta_pred <- fpca_fit$mu+allB%*%score_hat

df_pred_latent %>%
  filter(id==1) %>% 
  select(id, bin, eta_hat) %>% distinct(.) %>% 
  mutate(eta_pred=eta_pred) %>%
  ggplot()+
   geom_line(aes(x=bin, y = eta_hat))+
   geom_line(aes(x=bin, y = eta_pred), col="red")










#### Project eigen functions to a denser grid #####
# actually, used interpolation instead
keigen <- length(fpca_fit$evalues)

dense_efunc <- matrix(NA, nrow = J, ncol = keigen)
for(i in 1:keigen){
  dense_efunc[, i] <- approx(x=mid, y=fpca_fit$efunctions[ , i], xout=1:J)$y
}



