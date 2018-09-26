# packages
library('parallel') # Getting mclapply
library('glmnet') # for the Lasso
library('knitr') # knitting
library('Quandl') # data download
library('reshape2') # melting

# horizon setup
Rvec <- c(252)   # estimation window size 
#names(Rvec) <- c('one-year horizon')

# Load the data
source('../data')

# return matrix 
# monthly sample
# snp.returns <- readRDS('snp_returns_monthly_2006_2018.rds')
# snp.returns <- readRDS('snp_returns_monthly_1994_2018.rds')
# daily sample
 snp.returns <- readRDS('snp_returns_daily_2013_2018.rds')
# snp.returns  <- readRDS('snp_returns_daily_2016_2018.rds')
# snp.returns  <- readRDS('snp_returns_daily_2017_2018.rds')

# excess return matrix
# monthly sample
# snp.excess  <- readRDS('snp_excess_monthly_2006_2018.rds')
# snp.excess  <- readRDS('snp_excess_monthly_1994_2018.rds')
# daily sample
 snp.excess  <- readRDS('snp_excess_daily_2013_2018.rds')
# snp.excess  <- readRDS('snp_excess_daily_2016_2018.rds')
# snp.excess  <- readRDS('snp_excess_daily_2017_2018.rds')

#save as matrix
returns <- as.matrix(snp.excess)

# return target
#target <- 0.007974 # monthly expected return target (10% per year)
target <- 0.000378 # daily expected return target (10% per year)

# dimensions
nvec   <- nrow(returns)

