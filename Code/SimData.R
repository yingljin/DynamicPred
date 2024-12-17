
# This scripts include codes to simulated data
# corresponding to Section 4 in the manuscript

# clear workspace
rm(list=ls())

# set seed
set.seed(1090)

# load and, if necessary, install packages used in simulating data 
# (some of these you'll need for model fitting)
pckgs <- c("tidyverse","mgcv","refund","lme4")
invisible(
        sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
                install.packages(x)
                require(x, character.only=TRUE,quietly=TRUE)
        })
)

#### Simulate data ####

## basic parameters
# number of subjects
N <- 500
# number of observations per subject (observations per function)
J <- 1000
# functional domain, equally spaced grid on [0,1]
sind <- seq(0,1,len=J)

## simulate outcomes g(E[Y_i(s)]) = \eta_i(s) = f_0(s) + f_1(s)X_i + b_i(s)
# f_0(s) = 0
# b_i(s) = \sum_{k=1}^4 \phi_k(s)\xi_{ik}
#        \xi_{ik} \sim N(0, \lambda_k)

## fixed effects f_0, f_1
f_0 <- function(s) 0

## random effects
# eigenfunctions \phi_k(s) evaluated on the observed grid
phi <- sqrt(2)*cbind(sin(2*pi*sind),cos(2*pi*sind),
                     sin(4*pi*sind),cos(4*pi*sind))
# eigenvalues \lambda
K <- 4
lambda <- 0.5^(0:(K-1))
# subject-specific weights/coefficients \xi_ik
xi <- matrix(rnorm(N*K),N,K);
xi <- xi %*% diag(sqrt(lambda))
b_i <- xi %*% t(phi); # of size N by J

## linear predictor \eta_i(s)
eta_i <- t(vapply(1:N, function(x){
                        f_0(sind) + b_i[x,]
                }, numeric(J)))
## response Y_i(s) \sim Binom(p = expit(eta_i))
Y_i <- matrix(rbinom(N*J, size=1, prob=plogis(eta_i)), 
              N, J, byrow=FALSE)

## create data matrix in long format for estimation 
# id = subject identifier (factor variable)
# sind = numeric value corresponding to the observed functional domain
# Y = functional response (binary)
# sind_inx = numeric value associated with order of "sind"
#            this is not necessary, but may be of convenience when you implement the method
df <- data.frame(id = factor(rep(1:N, each=J)),
                 sind = rep(sind, N), 
                 Y = as.vector(t(Y_i)),
                 sind_inx = rep(1:J, N))


