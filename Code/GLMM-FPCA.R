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
fpca_fit <- fpca.face(mat_pred_unique, pve = 0.95, argvals = unique(df_pred_unique$bin), knots=20)
# fpca_fit <- face.sparse(df_pred_unique, argvals.new = 1:J, knots = 20, pve=0.95)


#### Project eigen functions to a denser grid #####
# actually, used interpolation instead
keigen <- length(fpca_fit$evalues)

dense_efunc <- matrix(NA, nrow = J, ncol = keigen)
for(i in 1:keigen){
  dense_efunc[, i] <- approx(x=mid, y=fpca_fit$efunctions[ , i], xout=1:J)$y
}



##### Prediction #####

# data with partically observed curves
# observations up to 1
df1 <- df %>% filter(id==1 & sind_inx<=100)

# ggplot(df1)+
#   geom_line(aes(x=sind_inx, y = eta_i, col = "red"))+
#   geom_point(aes(x=sind_inx, y = Y, col = "blue"), size = 0.01)+
#   xlim(1, J)


# score with empirical bayes
# fpca_fit$scores[1, ] ## scores from fPCA face
# 
# tao <- diag(fpca_fit$evalues) ## egenvalues, also var-cov matrix of scores
# B <- dense_efunc[1:100, ]
# 
# uhat <- tao %*% t(B) %*% solve(B %*% tao %*% t(B)) %*% df1$Y





