# This file contains the functions used to estimate large covariance matrices and out-of-sample forecasting. 

mciter <- function(nvec,snp.excess,snp.returns,returns,rf,target,Rvec,ic=NULL){
  
# Creating storage
sest.names <- c('POET','NodeWise','Ledoit-Wolf')
stat.names <- c('without transaction cost','with transaction cost')
rol.names <- c('Return','Variance','Sharpe','Turnover')
ptf.names <- c('Global','Markowitz')

ptf.stats  <- array(NA,dim = c(length(Rvec),length(rol.names),length(ptf.names),length(stat.names),length(sest.names)),
                    dimnames = list('R'=Rvec,'performance'=rol.names,'Portfolio'=ptf.names,'TC'=stat.names,'Estimator'=sest.names))

# sourcing constants
source('../code/constant.R')

# initial storage
gmv_smp <- mkw_smp <- gmv_poe <- mkw_poe <- gmv_ndw <- mkw_ndw <- gmv_lw <- mkw_lw <- NULL

# Looping over R dimension
ip <- 1 #counter

for(R in Rvec) {
  pred <- nvec-R         # pred is forecast horizon
for(j in 1:R)  {
  source('../code/data_generation.R')

    Y   <- returns[j:(pred+j-1),]
    n <- nrow(Y)
    p <- ncol(Y)
    
    # Nodewise information criterion (WIC, BIC, AIC)
    if(is.null(ic)) ic <- 'GIC'
    
    # COVARIANCE MATRIX COMPUTATIONS
    # Comptue the sample cov mat (and inverse)    
    ssmp <- est_smpcov(Y)
    # Compute the covariance matrix (and inverse) with UNKNWON FACTORS 
    spoe <- est_poetcov(ssmp[[1]],kmax = 7,n,p)
    # Compute the NODEWISE inverse covariance matrix
    sndw <- est_ndwcov(Y,ic)
    sndw[[1]] <- ssmp[[1]] # sample covariance matrix
    # Compute the LedoitWolf inverse covariance matrix
    slw  <- est_lwcov(Y)
    
    # PORTFOLIO COMPUTATIONS
    # expected return vectors
    mu_ndw <- matrix(colMeans(Y),ncol=1) 
    mu_smp <- mu_ndw
    mu_poe <- mu_ndw
    mu_lw  <- mu_ndw
   
    # portfolio variances  
    # sample
    if(class(ssmp[[2]])!='try-error') {
    smp_ptf <- portfolios(ssmp,target,mu_smp) 
    gmv_smp  <- cbind(gmv_smp,smp_ptf[[1]])
    mkw_smp  <- cbind(mkw_smp,smp_ptf[[2]])
    }
    # POET
    poe_ptf  <- portfolios(spoe,target,mu_poe)
    gmv_poe  <- cbind(gmv_poe,poe_ptf[[1]])
    mkw_poe  <- cbind(mkw_poe,poe_ptf[[2]])
    # nodewise
    ndw_ptf  <- portfolios(sndw,target,mu_ndw)
    gmv_ndw  <- cbind(gmv_ndw,ndw_ptf[[1]])
    mkw_ndw  <- cbind(mkw_ndw,ndw_ptf[[2]])
    # ledoitwolf
    lw_ptf  <- portfolios(slw,target,mu_lw)
    gmv_lw  <- cbind(gmv_lw,lw_ptf[[1]])
    mkw_lw  <- cbind(mkw_lw,lw_ptf[[2]])
  } 
  
  # global min and markowitz portfolios
  # performance for poet
  ptf_gmv_poe <- rolling(returns,gmv_poe,pred,R,rf,c=0.001)
  ptf_mkw_poe <- rolling(returns,mkw_poe,pred,R,rf,c=0.001)
  ptf.stats[ip,,1,1,'POET'] <- ptf_gmv_poe[[1]]
  ptf.stats[ip,,1,2,'POET'] <- ptf_gmv_poe[[2]]
  ptf.stats[ip,,2,1,'POET'] <- ptf_mkw_poe[[1]]
  ptf.stats[ip,,2,2,'POET'] <- ptf_mkw_poe[[2]]
  # performance for ndw
  ptf_gmv_ndw <- rolling(returns,gmv_ndw,pred,R,rf,c=0.001)
  ptf_mkw_ndw <- rolling(returns,mkw_ndw,pred,R,rf,c=0.001)
  ptf.stats[ip,,1,1,'NodeWise'] <- ptf_gmv_ndw[[1]]
  ptf.stats[ip,,1,2,'NodeWise'] <- ptf_gmv_ndw[[2]]
  ptf.stats[ip,,2,1,'NodeWise'] <- ptf_mkw_ndw[[1]]
  ptf.stats[ip,,2,2,'NodeWise'] <- ptf_mkw_ndw[[2]]
  # performance for lw
  ptf_gmv_lw <- rolling(returns,gmv_lw,pred,R,rf,c=0.001)
  ptf_mkw_lw <- rolling(returns,mkw_lw,pred,R,rf,c=0.001)
  ptf.stats[ip,,1,1,'Ledoit-Wolf'] <- ptf_gmv_lw[[1]]
  ptf.stats[ip,,1,2,'Ledoit-Wolf'] <- ptf_gmv_lw[[2]]
  ptf.stats[ip,,2,1,'Ledoit-Wolf'] <- ptf_mkw_lw[[1]]
  ptf.stats[ip,,2,2,'Ledoit-Wolf'] <- ptf_mkw_lw[[2]]
}
# End of R-loop, update counter
ip <- ip + 1

return(list('performance'=ptf.stats)) 
}


# A small standard error function
se <- function(x) sqrt(var(x))
# A small RMSE function
rmse <- function(x) sqrt(mean(x^2))
# A small MSE function
mse <- function(x) mean(x^2)
# A small Mean Abolute Error function
mnae <- function(x) mean(abs(x))
# A small Median Absolute Error function
mdae <- function(x) median(abs(x))
# A small trace function
tr <- function(x)sum(diag(x))

# Nodewise estimation of the covariance matrix
est_ndwcov <- function(Y,ic){
  
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  Y <- Y- t((apply(Y,2,mean))%*%matrix(1,1,n)) # Y is de-meaned
  
  C <- matrix(0,p,p)
  diag(C) <- 1
  tau <- NULL
  
  # Loop over the assets
  for(j in 1:p){
    # Estimate the Lasso
    jlas <- glmnet(x=Y[,-j],y=Y[,j],family = 'gaussian',intercept = FALSE)
    # Get fit
    jfit <- predict(jlas, newx=Y[,-j], type="response")    
    # residuals
    jres <- matrix(Y[,j],n,length(jlas$lambda)) - jfit
    # std err
    jsig <- colSums(jres^2)/n
    # Computing information criterion
    if(ic=='WIC') jbic  <- log(jsig) + jlas$df *log(n)/n * log(log(p)) # BIC (Wang,2010)
    if(ic=='GIC') jbic  <- log(jsig) + jlas$df *log(p)/n * log(log(n)) # GIC
    if(ic=='BIC') jbic  <- log(jsig) + jlas$df *log(n)/n  #BIC
    if(ic=='MIC') jbic  <- jsig + jlas$df *log(n) * log(log(p))/n # MC's IC
    if(ic=='AIC') jbic  <- log(jsig) + 2 * jlas$df # AIC
    # Index of selected model 
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[,jind]
    # Computing tau squared
    jtau <- sum(jres[,jind]^2)/n + (1/2)*jlas$lambda[jind]*sum(abs(jpar)) # using (10)
    # Storing the parameters
    C[j,-j] <- -jpar
    tau <- c(tau,jtau)
  }
  
  # Construct T-squared inverse
  T2inv <- diag(1/tau)
  
  # Construct Theta-hat
  Theta <- T2inv %*% C
  
  # sparsity
  sp <- sum(Theta==0)/(p^2)
  
  return(list(NULL,Theta,sp))
}


# ledoit wolf version of the squared frobenius norm
lwfrob <- function(x){sum(diag(tcrossprod(x)))/ncol(x)}

# Estimate Ledoit-Wolf covariance matrix
# From Ledoit Wolf JMA 2004
est_lwcov <- function(Y){
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  # sample cov mat
  S <- crossprod(Y)/n
  
  # from lemma 3.2
  m <- sum(diag(tcrossprod(S,diag(p))))/p
  # from lemma 3.3
  d2 <- lwfrob(S-m*diag(p))
  # from lemma 3.4
  b2bar <- sum(apply(Y,1,function(x,S){lwfrob(x%*%t(x)-S)},S=S))/n^2
  b2 <- min(d2,b2bar)
  # from lemma 3.5
  a2 <- d2-b2
  
  # Computing the LW covariance matrix
  siglw <- (b2/d2)*m*diag(p) + (a2/d2)*S
  # computing the inverse
  siglwinv <- chol2inv(chol(siglw))
  
  return(list(siglw,siglwinv,b2/d2,m))
}

# Sample covariance matrix estiamtion with inversion
est_smpcov <- function(Y){

  sigsmp <- var(Y)
  x <- try(chol(sigsmp),silent = TRUE)
  if(class(x)!='try-error')  sigsmpinv <- chol2inv(x)
  else sigsmpinv <- x
  
  return(list(sigsmp,sigsmpinv))  
}

# construct the partial matrix reconstruction for POET
sum_ev_back <- function(eval,evec,kmax){
  #init
  evsum <- list()
  # loop
  for(i in 1:(kmax+1)){
    tmp <- matrix(0,length(eval),length(eval))
    for(j in i:length(eval)){
      tmp <- tmp + eval[j]*(evec[,j]%*%t(evec[,j]))
    }
    evsum[[i]] <- tmp
  }
  return(evsum)
}

# construct the partial forward matrix reconstruction for POET with given khat
sum_ev_back_khat <- function(eval,evec,khat){
  #init
  tmp <- matrix(0,length(eval),length(eval))
  for(j in (khat+1):length(eval)){
    tmp <- tmp + eval[j]*(evec[,j]%*%t(evec[,j]))
  }
  return(tmp)
}

# construct the partial forward matrix reconstruction for POET with given khat
sum_ev_forward_khat <- function(eval,evec,khat){
  #init
  tmp <- matrix(0,length(eval),length(eval))
  for(j in 1:khat){
    tmp <- tmp + eval[j]*(evec[,j]%*%t(evec[,j]))
  }
  return(tmp)
}

# POET covariance matrix estimation
est_poetcov <- function(ssmp,kmax,n,p,Khat=NULL){
  
  # Eigen value/vector decomposition
  ev <- eigen(ssmp)
  eval <- ev$values
  evec <- ev$vectors
  
  if(is.null(Khat)){
    # Determining K-hat
    evsum_back <- sum_ev_back(eval,evec,kmax)
    Kvec <- (1/p) * sapply(evsum_back,tr) + (0:kmax)*(p+n)*log((p*n)/(p+n))/(p*n)
    Khat <- which.min(Kvec)-1
    # Construct Omega-hat
    Omegahat <- evsum_back[[Khat+1]]
  }
  if(!is.null(Khat)) Omegahat <- sum_ev_back_khat(eval,evec,Khat)
  
  # Construct 1st part of Sigma POET
  if(Khat==0)Spoet <- matrix(0,p,p)
  
  if(Khat!=0)Spoet <- sum_ev_forward_khat(eval,evec,Khat)  
  
  
  # Looping over C until covariance matrix invertible: 
  pinv <- TRUE
  Clev <- 0
  while(pinv){
    # Construct threshold matrix
    Cthres <- Clev*(sqrt(log(p)/n)+sqrt(1/p))
    domega <- matrix(diag(Omegahat),ncol=1)
    momega <- Cthres * sqrt(domega%*%t(domega))
    diag(momega) <- 0
    
    # hard thresholding
    # Omegahat[abs(Omegahat)<momega] <- 0
    # soft thresholding
    Omegahat[abs(Omegahat)<momega] <- 0
    Omegahat[abs(Omegahat)>=momega] <- sign(abs(Omegahat[abs(Omegahat)>=momega]))*(abs(Omegahat[abs(Omegahat)>=momega])-momega[abs(Omegahat)>=momega])
    
    
    # Constructing the poet estimator 
    Stmp <- Spoet + Omegahat
    
    # checking if positive definite
    try.inv <- try(chol(Stmp),silent = TRUE)
    pinv <- class(try.inv)=='try-error'
    if(pinv) Clev <- Clev + 0.1
    if(!pinv) {
      invSpoet <- try(chol2inv(try.inv))
      if(class(invSpoet)=='try-error') Clev <- Clev + 0.1
    }
    
  }
  
  # Storing final POET estimator
  Spoet <- Stmp
  
  return(list(Spoet,invSpoet,Clev,Khat))
}



# Compute the portfolio weights and variance
# target: a scalar with the target return
# mu: column vector of expected returns 
portfolios <- function(estsigm,target,mu){
  
  # Compute preliminary values
  p <- length(mu)
  one <- matrix(1,p,1)
  phi <- as.single(t(one) %*% estsigm[[2]] %*% one)
  psi <- as.single(t(one) %*% estsigm[[2]] %*% mu)
  Phi <- as.single(t(mu) %*% estsigm[[2]] %*% mu)
  
  # Global min variance
  gmv <- (estsigm[[2]] %*% one) / phi 
  gmvvar <- 1 / phi
  
  # Markowitz 
  denom  <- phi*Phi - psi^2
  mkw    <- (Phi-target*psi)/denom * estsigm[[2]] %*% one + (target*phi-psi)/denom * estsigm[[2]]%*%mu 
  mkwvar <- (phi*target^2 - 2*psi*target + Phi)/denom
  
  ptfret <- list(cbind(gmv),cbind(mkw),c(gmvvar),c(mkwvar))

  return(ptfret)
}


rolling <- function(excess,weight,pred,R,rf,c=0.001) {
  rets    <- t(excess[(pred+1):nrow(excess),])
  
  # compute vector for the portfolio return
  rewet = NULL
  # weight                                    # originally (p*R) dimension
  for (i in 1:R) {                            # (R*p) dimensional weight matrix 
    rewet[i] = t(weight)[i,] %*% rets[,i]     # (p*R) dimensional return matrix
  }
  # rewet = diag(t(weight) %*% rets)          
  
  # out-of-sample return of the porfolio 
  mur     <- mean(rewet)
  # out-of-sample variance 
  sigma2  <- var(rewet) 
  # out-of-sample sharpe ratio 
  sharpe  <- mur/sqrt(sigma2)
  
  # compute portfolio total return
  wexret  <- t(t(rets) + rf[(pred+1):nrow(excess)])              
  portret <- diag(t(weight) %*% wexret) 
  # portret <- rewet + rf[(pred+1):nrow(excess)]
  
  # asset return ratio in total portfolio return 
  retra <- sweep((1+rets), 2, (1+portret), `/`)

  # portfolio weight before rebalancing 
  wmins <- weight * retra
  # portfolio weight after rebalancing 
  wtpl  <- weight[,2:R]
  wmin  <- wmins[,1:(R-1)]
  
  # init storage
  turn = matrix((rep(NA, (R-1)*1)), nrow=(R-1))
  # compute turnover rate vector
  for (i in 1:(R-1)) { 
    turn[i] = sum(abs(wtpl[,i] - wmin[,i]))   # colSums(abs(wtpl - wmin))
     }
  # turnover rate for all assets 
  turnover <- mean(turn)
  

# return with transaction cost
  reventc        <- matrix((rep(NA, (R-1))))
  for (i in 1:length(reventc)){ 
  reventc[i]     <- rewet[i] - (c * (1+rewet[i])) * (sum(abs(wtpl[,i] - wmin[,i])))
  }
  # oos return with transaction cost   
  mutc = mean(reventc) 
  # oos variance with transaction cost  
  sigtc = var(reventc)
  # oos sharpe with transaction cost  
  shartc = mutc/sqrt(sigtc)
  
  perest <- list(c(mur,sigma2,sharpe,turnover),c(mutc,sigtc,shartc,NA))

  return(perest)
  
}

