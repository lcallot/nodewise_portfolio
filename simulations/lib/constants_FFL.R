# This file contains the constants used in Fan, Fan, Lv 2008 Journal of Econometrics
# These values are sourced by the scripts and knitr files. 
# Laurent Callot (l.callot@gmail.com), 14/01/2016


# Moments for the factors
muf  <- c(0.023558,0.012989,0.020714)
covf <- matrix(c(1.2507,-0.034999,-0.20419,
                 -0.034999,0.31564,-0.0022526,
                 -0.20419,-0.0022526,0.19303),
               3,3)

# Moments for the loadings
mub  <- c(0.78282,0.51803,0.41003)
covb <- matrix(c(0.029145,0.023873,0.010184,
                 0.023873,0.053951,-0.006967,
                 0.010184,-0.006967,0.086856),
               3,3)


# Parameters of the Gamma distribution
alpha <- 3.3586
beta  <- 0.1876

