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

# save as matrix
returns <- as.matrix(snp.excess)

# risk free rate
# rf  <- (snp.returns[,1]/100)/12   # convert to monthly return
 rf  <- (snp.returns[,1]/100)/252  # convert to daily return

# dimensions
nvec   <- nrow(returns)


