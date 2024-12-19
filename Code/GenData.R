
# This sciprt generateds 500 data sets for the simulation study

set.seed(1114)

library(mvtnorm)

N_train <- 500 # training sample size
N_test <- 100 # testing sample size
N <- N_train+N_test

J <- 1000 # number of observation points
t = seq(0,1,len=J) # observation grid 

# fixed effect: mean function
f_0 <- function(s){return(0)}

# random effects
K <- 4 # number of eigenfunctions

phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
                     sin(4*pi*t),cos(4*pi*t)) # orthogonal random functions

lambda = 0.5^(0:(K-1)) # eigenvalues/variance of scores

sind = seq(0,1,len=J) # observations points

# visualize eigenfunctions
par(mfrow=c(2, 2))
plot(t, phi[,1])
plot(t, phi[,2])
plot(t, phi[,3])
plot(t, phi[,4])



#### generated data ####

sim_data <- list() # container
M <- 500 # number of datasets to generate

for(m in 1:M){
  
  # generate score
  xi <- rmvnorm(N, mean = rep(0, K), sigma = diag(lambda))
  bi <- xi %*% t(phi) # random effect functions

  # latent gaussian function
  eta_i <- f_0(t)+bi
  
  this_df <- data.frame(id = factor(rep(1:N, each=J)),
                        t = rep(t, N), 
                        eta_i = as.vector(t(eta_i)),
                        sind = rep(1:J, N))
  
  this_df$Y <- rbinom(N*J, size=1, prob=plogis(this_df$eta_i))
  
  sim_data[[m]] <- this_df

}



# visualization

rand_id <- sample(N, size = 4)

sim_data[[357]] %>% filter(id %in% rand_id) %>% 
  ggplot()+
  geom_point(aes(x=t, y=Y))+
  geom_line(aes(x=t, y=plogis(eta_i)), col = "red")+
  facet_wrap(~id)

# save
save(sim_data, file = here("Data/sim_data.RData"))
 
