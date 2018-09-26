# The script saves the portfolio estimation results which are processed to produce the output in a .Rnw file. 

library('parallel') # Getting mclapply
library('MASS') # for multivariate normal
library('glmnet') # for the Lasso
library('zoo') # for dates
library('lubridate') # for date class
library('scales')

# Number of CPUs to use
ncores <- 2

# Seed
set.seed(1453)

# setting directory to script directory
setwd("../application")

# Get the constants
source('../code/constant.R')

# Get the library
source('../code/roll.R')

# Do the work
gpmc <- mciter(nvec,snp.excess,snp.returns,returns,rf,target,Rvec,ic=NULL)

# save 
save(gpmc,file ='Ndw.Rda')


