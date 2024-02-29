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

# 24/1/11
Data application
1. Run on the whole data and see what happens with eigenvalues
2. Percentage variance explained by each PC and time domain

--- Ask for the paper where conclusion lambda = 1/sp
--- What happened to the fGFPCA paper estiamtes? 

2nd reference method
1. Actually, for each interval, use the observation immediately before the start of the interval as a time-fixed observation. All predictions will be made based on this single observation. This is a function-on-scalar regression


# 24/1/18
1. Send Andrew the NHANES document, with summary statistics and histograms of score distribution
- The score for PC3 is very skewed (negatively) on binned data. It may have drawn variation from the night time. 
- Also, report variantion explained by each PC
- Ask Suni about the conclusion linking eigenvalues and smoothing parameters. 

# 24/2/1

1. Flip coordinate of tables (landscape) and figures. Let make each row the same subject, and compared across methods
2. Run the data application with all four methods
3. For the reversed eigenvalue problem: generate data using estimates from fGFPCA (eigenvalues, eigenfunctions, mean, etc). Can we reproduce the reverse eigenvalue problem? 

# 24/2/15
1. Add citation about preprocessing of binary indicator
2. Change figure line thickness and color, condense space
3. Figure 3: yellow line for fGFPCA looks strange. Double check
4. Write an appendix explaining the reverse eigenvalues problem

# 24/2/22
1. NHANES example: get distribution of scores, population average probability curve
2. Reverse eigenvalue problem:
- get the probability plot of mu(t)+2*\lambda_2*phi(t) for 2nd and 2rd PCs. Bascially, we wanna see how this probability affect the observed outcome.

# 24/2/29
1. Individual scores: calculated the mean
2. Do iterative simulation and summarize the eigenfunctions and eigenvalues, see when they'd flip
3. A few experienments:
- Put a lower bound on population mean, so that the nighttime observations are about half active and inactive. Use the same PC, simulate data, do the fGFPCA, see if eigenvalues still flip
- Upsample with replacement (e.g. 5N), do fGFPCA and see if PC2 and PC3 flip and if the PC stablize
- Maybe try more than 4 eigenfunctions? 




