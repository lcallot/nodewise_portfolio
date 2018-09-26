# This file contains the code to replicate the simulations in A Nodewise Regression Approach to Estimating Large Portfolios
# The script saves an object with statistics which are processed in a .Rnw file to produce the output. 
# Laurent Callot (l.callot@gmail.com), 29/04/2016


library('parallel') # Getting mclapply
library('MASS') # for multivariate normal
library('mvtnorm') # for multivariate t
library('glmnet') # for the Lasso
library('msgps') # for the Lasso
library('rPython') # to call the python sparse_spd_matrix method

# Number of CPUs to use
ncores <- 64

# Seed
set.seed(1453)

# Number of iterations
niter <-1000

# setting directory to script directory
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

# Get the library
source('../lib/lib_sim.R')

# Get the DGP functions
source('../lib/lib_DGP.R')

# Number of observations (sample length)
nobs <- 252

# target vector (yearly,monthly,daily)
tgvec <- c(0.000378)
names(tgvec) <- c('daily')

# Exposure vector
expvec <- seq(1,2,length.out = 3)
exprep <- 2

# Sequence of model sizes
# pvec <- seq(from = 50, to = 300, by = 50)
pvec <- c(seq(from = 25, to = 100, by = 25), seq(from = 150, to = 300, by = 50), seq(from = 400, to =500, by = 100))

# Do the work
mc <- mclapply(1:niter,mciter,nobs,pvec,tgvec,expvec,exprep,ic='GIC',mc.cores = ncores)
# mc <- lapply(1:niter,mciter,nobs,pvec,tgvec,expvec,exprep,ic='GIC')

# save 
saveRDS(mc,file = 'sim_results_2018_GIC_1000.rds')
