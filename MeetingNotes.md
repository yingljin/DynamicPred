# DynamicPred Meeting Notes

# 23/11/13
1. For Laplace approximation, try a few other optimizers: NM, Newton-Ralphson
2. Find different packages 
3. Set up overleaf document for manuscript 

# 23/11/20
1. Try GLMMadaptive() with ns() wrapped inside
2. Send Andrew overleaf document

# 23/11/30
For the small scale simulation:
1. Sample size N_train and N_test should at least be 100, and number of simulations should be 500. We can reduce the flexibility of of GLMMadaptive so that the model fits within 30 minutes.
2. fGFPCA should be fit on the original data
3. GLMMadaptive should be fit on a reduced density (e.g take every one observation every 10 observation). Prediction should be done on the original scale.  
4. Look into how GLMMadaptive makes prediction

Presentation notes: 
1. Move the exponential family formula later, right before the posterior likelihood
2. Step 4: change "projection" to "re-evaluation"
3. Add some data application results
3. Figure before 

# 23/12/14
1. Second reference method: Generalized Historical model
g(E[Y_i(s)|Y_i(t), t<s]) = \beta_0(s) + \beta_1 (s)* Y_i(t) for s > t
2. Run 500 iterations for the small-scale simulation
3. Fit on the who data set on NHANES (30%-40% for testing maybe?). The debias model may take more than a day to fit

Qeustions:
1. Essentially, what is the difference between treating Y as a series of binary values or a function?





