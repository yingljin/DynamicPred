
# This sciprt generateds 500 data sets

set.seed(314)
N <- 500 # sample size
J <- 1000 # number of observation points

sind = seq(0,1,len=J) # observations points


f_0 <- function(s) 0 # mean functions
f_0(5)


#### eigenfunctions #####
K <- 4 # number of eigenfunctions
phi <- sqrt(2)*cbind(sin(2*pi*sind),cos(2*pi*sind),
                     sin(4*pi*sind),cos(4*pi*sind))

par(mfrow=c(2, 2))
plot(sind, phi[,1])
plot(sind, phi[,2])
plot(sind, phi[,3])
plot(sind, phi[,4])


#### generated data ####

sim_data <- list()
M <- 500 # number of datasets to generate
lambda = 0.5^(0:(K-1)) # eigenvalues

for(m in 1:M){
  
  # generate score
  xi <- matrix(rnorm(N*K),N,K)
  xi <- xi %*% diag(sqrt(lambda))
  
  # subject-specific random effect
  b_i <- xi %*% t(phi); # of size N by J
  
  # latent gaussian function
  eta_i <- t(vapply(1:N, function(x){
    f_0(sind) + b_i[x,]
  }, numeric(J)))
  
  # outcome binary function
  Y_i <- matrix(rbinom(N*J, size=1, prob=plogis(eta_i)), 
                N, J, byrow=FALSE)
  
  # format into dataframe
  # id = subject identifier (factor variable)
  # sind = numeric value corresponding to the observed functional domain
  # Y = functional response (binary)
  # sind_inx = numeric value associated with order of "sind"
  #            this is not necessary, but may be of convenience when you implement the method
  df <- data.frame(id = factor(rep(1:N, each=J)),
                   sind = rep(sind, N), 
                   Y = as.vector(t(Y_i)),
                   eta_i = as.vector(t(eta_i)),
                   sind_inx = rep(1:J, N))
  
  sim_data[[m]] <- df

}

# save
save(sim_data, file = here("Data/sim_data.RData"))

# visualization

rand_id <- sample(N, size = 4)

sim_data[[1]] %>% filter(id %in% rand_id) %>% 
  ggplot()+
  geom_point(aes(x=sind_inx, y=Y))+
  geom_line(aes(x=sind_inx, y=plogis(eta_i)), col = "red")+
  #geom_line(aes(x=sind_inx, y=eta_i), col = "blue")+
  facet_wrap(~id)

 
