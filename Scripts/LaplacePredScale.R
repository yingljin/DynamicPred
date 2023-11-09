
# I want to look at the posterior distribution from scaled and unscaled estiamte
# also figure out why scaled parameters would lead to such low convergence rate
# with scaled parameter, the posterior distribution of scores would be much, much more centered
# does it mean that sampling got stuck in a neighborhood of false values and was never able to get close to the "true" values? 

this_id <- 43
this_t <- 0.4

# data
df_it <- test_df %>% filter(id==this_id & t <= this_t)
max_rid <- nrow(df_it)

#### use Unscaled params ####
MyData <- list(K=K, 
               X=efunctions_new[1:max_rid, ], 
               mon.names=mon.names,
               parm.names=parm.names, 
               pos.xi=pos.xi, 
               y=df_it$Y, 
               tao=diag(new_lambda), f0=new_mu[1:max_rid])
Fit <- LaplaceApproximation(Model, parm = rep(0, K), Data=MyData, Method = "BFGS", Iterations = 1000,
                            CovEst = "Identity")
Fit$Converged # converged
## posterior samples
Fit$Posterior %>% as.data.frame() %>%
  pivot_longer(starts_with("xi")) %>%
  ggplot()+
  geom_histogram(aes(x=value), bins = 50)+
  facet_wrap(~name)

#### use Scaled params ####
MyData2 <- list(K=K, 
               X=efunctions_scaled[1:max_rid, ], 
               mon.names=mon.names,
               parm.names=parm.names, 
               pos.xi=pos.xi, 
               y=df_it$Y, 
               tao=diag(scaled_lambda), f0=new_mu[1:max_rid])
Fit_scaled <- LaplaceApproximation(Model, parm = rep(0, K), 
                                   Data=MyData2, Method = "BFGS", Iterations = 1000,
                            CovEst = "Identity")
Fit_scaled$Converged # not converged

## posterior samples
Fit$History %>% as.data.frame() %>%
  mutate(iter = 1:nrow(Fit$History)) %>%
  pivot_longer(starts_with("xi")) %>%
  ggplot()+
  geom_point(aes(x=iter, y=value))+
  facet_wrap(~name)

Fit_scaled$History %>% as.data.frame() %>% 
  mutate(iter = 1:nrow(Fit_scaled$History)) %>%
  filter(iter > 10) %>%
  pivot_longer(starts_with("xi")) %>%
  ggplot()+
  geom_point(aes(x=iter, y=value))+
  facet_wrap(~name, scales = "free")

xi_test[43, ] # true scorer: 1.1198069  0.4869535 -0.8083169  0.3692140


Fit_scaled$Posterior %>% as.data.frame() %>% 
  pivot_longer(starts_with("xi")) %>%
  ggplot()+
  geom_histogram(aes(x=value), bins = 50)+
  facet_wrap(~name)
  
plot(Fit, Data=MyData)
plot(Fit_scaled, Data=MyData2)

Fit$Summary1
Fit_scaled$Summary1

Fit$History %>% head()

