library(rstan)

#### true values ####

# fixed effect: mean function
true_mu <- rep(0, J)

#eigenfunctions
true_phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
                     sin(4*pi*t),cos(4*pi*t)) # orthogonal random functions

true_lambda = 0.5^(0:(K-1)) # eigenvalues/variance of scores



#### stan ####
test_df <- sim_data[[1]] %>% filter(id %in% 501:600)
test_id <- unique(test_df$id) 
tmax_vec <- c(0.2, 0.4, 0.6, 0.8)
df_pred <- list()

for(id in seq_along(test_id)){
  m <- test_id[[id]]
  df_i <- test_df %>% filter(id==m)
  
  scores <- matrix(NA, K, length(tmax_vec))
  sd_scores <- matrix(NA, K, length(tmax_vec))
  
  for(i in seq_along(tmax_vec)){
    
    tmax <- tmax_vec[i]
    df_it <- df_i %>% filter(t <= tmax)
  
    # into a list
    stanData <- list(
      J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
      efuncs = true_phi, 
      b0 = true_mu,
      lambda = true_lambda
    )
    
    fit <- stan(
      file = here("Code/prediction.stan"),  # Stan program
      data = stanData,    # named list of data
      chains = 2,             # number of Markov chains
      warmup = 500,          # number of warmup iterations per chain
      iter = 1000,            # total number of iterations per chain
      cores = 1,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )
    
    scores[ ,i] <- summary(fit)$summary[1:K, "mean"]
    sd_scores[ ,i] <- summary(fit)$summary[1:K, "sd"]
    
    # try a different algorithm
    # fit2 <- stan(
    #   file = here("Code/prediction.stan"),  # Stan program
    #   data = stanData,    # named list of data
    #   chains = 2,             # number of Markov chains
    #   warmup = 500,          # number of warmup iterations per chain
    #   iter = 1000,            # total number of iterations per chain
    #   cores = 1,              # number of cores (could use one per chain)
    #   refresh = 0,
    #   algorithm = "HMC"# no progress shown
    # )
    # 
    # scores2[ ,i] <- summary(fit2)$summary[1:K, "mean"]
    # sd_scores2[ ,i] <- summary(fit2)$summary[1:K, "sd"]
    
    
    # latent function predictions
    eta_pred_out <- true_mu+true_phi%*%scores[ ,i]
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
    
    # prediction interval
    sd_eta <- sqrt((true_phi^2) %*% sd_scores[ ,i]^2)
    df_i[ , paste0("pred", tmax, "_lb")] <- as.vector(eta_pred_out[, 1]-qnorm(0.975)*sd_eta)
    df_i[ , paste0("pred", tmax, "_ub")] <- as.vector(eta_pred_out[, 1]+qnorm(0.975)*sd_eta)
  }
  
  df_pred[[id]] <- df_i

}


plot(fit)
summary(fit)
traceplot(fit)
Fit$Summary1

scores
sd_scores

#### Prediction interval ####

View(df_pred[[1]])

# figure
rand_id <- sample(test_id, 4)

df_exp <- df_pred %>%
  bind_rows() %>%
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))})  
df_exp[df_exp$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
df_exp[df_exp$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
df_exp[df_exp$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
df_exp[df_exp$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA

df_exp %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.2)+
  geom_line(aes(x=t, y=eta_i, col = "True"))+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.2_lb, ymax = pred0.2_ub,
                  col = "0.2", fill = "0.2", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.4_lb, ymax = pred0.4_ub,
                  col = "0.4", fill = "0.4", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.6_lb, ymax = pred0.6_ub,
                  col = "0.6", fill = "0.6", alpha = 0.1),
              linetype="dashed")+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), na.rm = T)+
  geom_ribbon(aes(x=t, ymin=pred0.8_lb, ymax = pred0.8_ub,
                  col = "0.8", fill = "0.8", alpha = 0.1),
              linetype="dashed")+
  facet_grid(rows = vars(id))+
  guides(alpha = "none", col="none")
ggsave(here("Images/IntervalExp.jpeg"), height=, width = 12)

# coverate rate
list_cover <- list()

df_cover <- df_exp %>% 
    mutate(
      cover0.2 = pred0.2_lb<=eta_i & pred0.2_ub>=eta_i,
      cover0.4 = pred0.4_lb<=eta_i & pred0.4_ub>=eta_i,
      cover0.6 = pred0.6_lb<=eta_i & pred0.6_ub>=eta_i,
      cover0.8 = pred0.8_lb<=eta_i & pred0.8_ub>=eta_i
    ) %>% 
    group_by(t) %>% 
    summarize_at(vars(starts_with("cover")), mean)
}


df_cover %>% pivot_longer(starts_with("cover")) %>%
  ggplot()+
  geom_line(aes(x=t, y=value, col=name))+
  ylim(0, 1)
ggsave(here("Images/IntervalCov.jpeg"), height=4, width = 4)
