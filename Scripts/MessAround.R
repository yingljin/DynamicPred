library(GLMMadaptive)

set.seed(12345)

##### simulation #####

n <- 500 # number of subjects
K <- 10 # number of measurements per subject
t_max <- 4 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

lapply(DF, summary)

# design matrices for the fixed and random effects non-zero part
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ 1, data = DF)
# design matrices for the fixed and random effects zero part
X_zi <- model.matrix(~ sex, data = DF)
Z_zi <- model.matrix(~ 1, data = DF)

betas <- c(0.8, -0.5, 0.8, -0.5) # fixed effects coefficients non-zero part
shape <- 2 # shape/size parameter of the negative binomial distribution
gammas <- c(-2.5, 0.5) # fixed effects coefficients zero part
D11 <- 1.0 # variance of random intercepts non-zero part
D22 <- 0.8 # variance of random intercepts zero part

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor non-zero part
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, 1, drop = FALSE]))
# linear predictor zero part
eta_zi <- as.vector(X_zi %*% gammas + rowSums(Z_zi * b[DF$id, 2, drop = FALSE]))
# we simulate negative binomial longitudinal data
DF$y <- rnbinom(n * K, size = shape, mu = exp(eta_y))
# we set the extra zeros
DF$y[as.logical(rbinom(n * K, size = 1, prob = plogis(eta_zi)))] <- 0


# split data 
ids_train <- sample(500, 250)
DF_train <- DF[DF$id %in% ids_train, ]
DF_test <- DF[!DF$id %in% ids_train, ]


##### Model #####

fm1 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF_train,
                   family = poisson())

fm2 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF_train,
                   family = zi.negative.binomial(), 
                   zi_fixed = ~ sex, zi_random = ~ 1 | id)
preds_fm1 <- predict(fm1, newdata = DF_test[DF_test$time < 2, ],
                     newdata2 = DF_test[DF_test$time >= 2, ], 
                     type = "subject_specific",
                     se.fit = TRUE, return_newdata = TRUE)

preds_fm2 <- predict(fm2, newdata = DF_test[DF_test$time < 2, ],
                     newdata2 = DF_test[DF_test$time >= 2, ], 
                     type = "subject_specific",
                     se.fit = TRUE, return_newdata = TRUE)
