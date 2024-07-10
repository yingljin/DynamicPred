library(rstan)

test_id <- 501:504
tmax_vec <- c(0.2, 0.4, 0.6, 0.8)
df_pred <- list()

for(m in test_id){

  df_i <- test_df %>% filter(id==m)
  
  scores <- matrix(NA, K, length(tmax_vec))
  sd_scores <- matrix(NA, K, length(tmax_vec))
  
  for(i in seq_along(tmax_vec)){
    
    tmax <- tmax_vec[i]
    df_it <- df_i %>% filter(t <= tmax)
  
    # into a list
    stanData <- list(
      J = J, Ju = nrow(df_it), Y = df_it$Y, K = K,
      efuncs = new_phi, 
      b0 = new_mu,
      lambda = new_lambda
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
    eta_pred_out <- new_mu+new_phi%*%scores[ ,i]
    df_i[ , paste0("pred", tmax)] <- eta_pred_out[,1]
    
    # prediction interval
    sd_eta <- sqrt((new_phi^2) %*% sd_scores[ ,i]^2)
    df_i[ , paste0("pred", tmax, "_lb")] <- as.vector(eta_pred_out[, 1]-qnorm(0.975)*sd_eta)
    df_i[ , paste0("pred", tmax, "_ub")] <- as.vector(eta_pred_out[, 1]+qnorm(0.975)*sd_eta)
  }
  
  df_pred[[m-500]] <- df_i

}


plot(fit)
summary(fit)
traceplot(fit)
Fit$Summary1

scores
sd_scores

#### Prediction interval ####

View(df_pred[[1]])

# prediction results
df_exp <- bind_rows(df_pred) %>% 
  # filter(id %in% sample(test_id, 4)) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            function(x){exp(x)/(1+exp(x))}) 
df_exp[df_exp$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
df_exp[df_exp$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
df_exp[df_exp$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
df_exp[df_exp$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA

df_exp %>% filter(t>0.78) %>%select(starts_with("pred0.8"))  %>% View()
# df_exp <- df_exp %>% filter(id == 568)
# figure with prediction interval
ggarrange(
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.2), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.2_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.2_ub), linetype="dashed",na.rm=T, col = "red")+
    facet_wrap(~id, ncol = 1),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.4), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.4_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.4_ub), linetype="dashed",na.rm=T, col = "red")+
    facet_wrap(~id, ncol = 1),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.6), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.6_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.6_ub), linetype="dashed",na.rm=T, col = "red")+
    facet_wrap(~id, ncol = 1),
  
  df_exp %>% 
    ggplot()+
    geom_point(aes(x=sind, y=Y), size = 0.2)+
    geom_line(aes(x=sind, y=eta_i))+
    geom_line(aes(x=sind, y=pred0.8), na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.8_lb), linetype="dashed", na.rm=T, col = "red")+
    geom_line(aes(x=sind, y=pred0.8_ub), linetype="dashed",na.rm=T, col = "red")+
    facet_wrap(~id, ncol = 1),
  nrow = 1
)
ggsave(here("Images/IntervalExp.jpeg"), height=3, width = 12)
