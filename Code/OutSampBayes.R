# This script writes functions for 
# calculating conditional expectation
# of fPCA scores 
# of an out-of-sample observation


## model
Model <- function(parm, Data){
  xi <- parm[Data$pos.xi]
  
  # log-prior
  xi.prior <- dmvnorm(xi, mean = rep(0, Data$K), sigma=Data$tao, log = TRUE)
  
  # log-posterior likelihood
  eta <- Data$f0+Data$X %*% xi
  p <- exp(eta)/(1+exp(eta))
  LL <- sum(dbern(x=Data$y, prob=p, log = TRUE)) # log likelihood of Y|xi
  LP <- LL+sum(xi.prior) # joint log likelihood of (Y, xi)
  
  # output
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LL,
                   yhat=Data$y, parm=parm)
  return(Modelout)
}

# parameter names
name_lst <- as.list(rep(0, K))
names(name_lst) <- paste("xi", 1:K, sep = "")
parm.names <- as.parm.names(name_lst)
pos.xi <- grep("xi", parm.names)
mon.names <- "LP"

# 
# ## a Bernoulli distribution at each time point
# Model <- function(parm, Data){
#   xi <- parm[Data$pos.xi]
#   
#   # log-prior
#   xi.prior <- dmvnorm(xi, mean = rep(0, Data$J), sigma=Data$tao, log = TRUE)
#   
#   # log-posterior likelihood
#   eta <- Data$f0+Data$X %*% xi
#   p <- exp(eta)/(1+exp(eta))
#   LL <- sum(dbern(x=Data$y, prob=p, log = TRUE)) # log likelihood of Y|xi
#   LP <- LL+sum(xi.prior) # joint log likelihood of (Y, xi)
#   
#   # output
#   Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LL,
#                    yhat=Data$y, parm=parm)
#   return(Modelout)
# }


# # data: 
# out_pred_laplace <- function(mu, evalues, phi_mat, df_new, kpc, max_t){
#   
#   # put data into correct format
#   # ns <- as.vector(table(df_new$bin)) # number of observations
#   # hs <- df_new %>% group_by(bin) %>% summarize_at("Y", sum) %>% select(Y) %>% unlist()# number of success
#   # nf <- ns-hs # number of failure
#   # max_bin <- length(unique(df_new$bin)) # assume new skipped bins  
#   # df_new2 <- data.frame(bin = unique(df_new$bin), ns, hs, nf)  
#   
#   # truncate data
#   df_new <- df_new[df_new$t<=max_t, ]
#   y <- df_new$Y[df_new$t<=max_t] 
#   
#   # into a list
#   tao <- diag(evalues)
#   f0 <- mu[1:max_t*J]
#   # N <- nrow(df_new2) # "sample size" which is in fact number of observed bins in this case for a new subject
#   # n <- df_new2$ns # number of experiements at each bin
#   y <- df_new$Y[df_new$t<=max_t] # outcome, number of 1
#   X <- phi_mat[1:max_bin, ] # "covariates", in this case PC functions
#   J <- kpc # number of parameters/scores
#   mon.names <- "LP"
#   
#   # parameter names
#   name_lst <- as.list(rep(0, J))
#   names(name_lst) <- paste("xi", 1:J, sep = "")
#   parm.names <- as.parm.names(name_lst)
#   pos.xi <- grep("xi", parm.names)
#   PGF <- function(Data) {
#     xi <- rnorm(Data$J)
#     return(xi)
#   }
#   
#   MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
#                  parm.names=parm.names, pos.xi=pos.xi, y=y, n=n, tao=tao, f0=f0)
#   
#   
#   # fit laplace approximation
#   Fit <- LaplaceApproximation(Model, parm = rep(0, J), Data=MyData, Method = "BFGS")
#   score <- Fit$Summary1[, "Mode"]
#   
#   # prediction
#   eta_pred_out <- mu+phi_mat%*%score
#   
#   return(list(eta_pred = eta_pred_out,
#                score_out = score))
# }



# fit LaplacesDemon model
# Fit <- LaplacesDemon(Model, Data=MyData, rep(0, 4),
#                      Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#                      Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

# test

# score_out <- out_pred_laplace(fpca_fit,  df %>% filter(id==rand_id[1] & sind<=720))$score_out
# eta_pred_out <- fpca_mod$mu+fpca_mod$efunctions%*%score_out

