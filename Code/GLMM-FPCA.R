##### Local GLMMs #####

bin_w <- 10 # bin width
n_bin <- J/bin_w # number of bins
brks <- seq(0, J, by = bin_w) # cutoff points
mid <- (brks+bin_w/2)[1:n_bin] # mid points

# bin observations
df$bin <- cut(df$sind_inx, breaks = brks, include.lowest = T, labels = mid)
df$bin <- as.numeric(as.character(df$bin))

# check overlap
# df %>% 
#   filter(id %in% 1:4) %>% 
#   ggplot()+
#   geom_line(aes(x=sind, y=eta_i, col=as.factor(bin)), show.legend = F)+
#   geom_point(aes(x=sind, y=Y, col=as.factor(bin)), size = 0.05, show.legend = F)+
#   facet_wrap(~id)

# fit local linear mixed models in each bin
# and extract subject level estimation
## df: function with binary outcome Y, with observations in one time bein
pred_latent <- function(df){
  this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial)
  eta_hat <- predict(this_glm, type = "link")
  df$eta_hat <- eta_hat
  return(df)
}

# additionally, get the average estimation across subject 
# by extacting the random intercepts from each bin
mean_latent <- function(df){
  this_glm <- glmer(Y ~ 1 + (1|id), data = df, family = binomial)
  beta_hat <- coef(summary(this_glm))[1]
  return(beta_hat)
}

## do that for all time bin
df_bin_lst <- split(df, f = df$bin)

###### subject specific
df_pred_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
# lapply(df_pred_latent, dim)
df_pred_latent <- bind_rows(df_pred_latent)

###### subject average
mean_latent <- sapply(df_bin_lst, function(x){mean_latent(x)}) 

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






