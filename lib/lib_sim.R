# This file contains the functions used to estimate large covariance matrices. 
# Laurent Callot (l.callot@gmail.com), 29/04/2016



# Main simulation function
# i: iteration index (for compatibility with mclapply)
# nobs: integer with sample length
# pvec: a vector with the number of assets
# tgvec: the vector of return targets
# expvec: the vector of exposures for the random weights
# exprep: the number of random portfolios to generate for each iteration
# ic: the information criterion for the penalty parameter selection, default to WIC. Other options: AIC, BIC, MIC.
mciter <- function(i,nobs,pvec,tgvec,expvec,exprep,ic=NULL){
  
  # timer
  tstart <- proc.time()
  
  # Creating storage
  sest.names <- c('Sample','Factor','POET','NodeWise','Ledoit-Wolf')
  dgp.names <- c('FFL 2008','S&P Student')
  ptf.names <- c(paste0('Optimal, ',names(tgvec),' target'),'Global min variance')
  stat.names <- c('Risk','Weights','Exposure','Risk_diff','Risk_ratio')
  norm.names <- c('Frobenius','Sigma norm','Entropy loss','Mean AE','Median AE','Frobenius inverse')
  misc.names <- c('POET C','POET Khat','NodeWise Lem 1','Sparsity','LW ratio','LW scale')
  comp.names <- c('Estimation','Exposure')
  expo.names <- c(paste0('Exposure: ',expvec))
  exst.names <- c('Risk','Weights')
  
  norm.stats <- array(NA,dim = c(length(pvec),length(norm.names),length(sest.names),length(dgp.names)),
                      dimnames = list('p'=pvec,'norm'=norm.names,'sest'=sest.names,'dgp'=dgp.names))
  
  expo.stats <- array(NA,dim = c(length(pvec),length(ptf.names),length(expo.names),length(sest.names),length(dgp.names),length(exst.names)),
                      dimnames = list('p'=pvec,'portfolio'=ptf.names,'exposure'=expo.names,'sest'=sest.names,'dgp'=dgp.names,'stats'=exst.names))
  
  ptf.stats <- array(NA,dim = c(length(pvec),length(ptf.names),length(sest.names)+1,length(dgp.names),length(stat.names)),
                     dimnames = list('p'=pvec,'portfolio'=ptf.names,'sest'=c('True',sest.names),'dgp'=dgp.names,'stats'=stat.names))
  
  misc.stats <- array(NA,dim = c(length(pvec),length(misc.names),length(dgp.names)),
                      dimnames = list('p'=pvec,'misc'=misc.names,'dgp'=dgp.names))
  
  time.stats <- array(NA,dim = c(length(pvec),length(sest.names)+2,length(comp.names),length(dgp.names)),
                      dimnames = list('p'=pvec,'times'=c('True',sest.names,'data generation'),'comp'=comp.names,'dgp'=dgp.names))
  
  # Generate the factors for FFL: 
  # loading the DGP moments
  source('../lib/constants_FFL.R')
  fac_FFL <- mvrnorm(n=nobs,mu=muf,Sigma=covf)
  
  # Looping over the p dimension
  ip <- 1 #counter
  for(p in pvec)
  {
    for(dgp in dgp.names)
    {
      
      if(dgp=='FFL 2008'){
        # Generating data from FFL, 
        time.stats[ip,'data generation',1,dgp] <- system.time(ffl <- FFL_DGP(fac_FFL,p))[1]
        Y <- ffl[[1]]
        stru <- ffl[2:3]
        mu_tru <- ffl[[4]]
      }
      
      if(dgp=='S&P Block Diag'){
        # Generating data from the S&P moments, with block diagonal matrix 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- Sparse_block_diag_DGP(p,nobs,ndist = 'norm',sparsity = 0.9))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
      }
      
      if(dgp=='LW Gaussian'){
        # Generating data from the S&P moments, with block diagonal matrix 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- LW_2004_DGP(p,nobs,ndist = 'norm'))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
      }
      
      if(dgp=='LW Student'){
        # Generating data from the S&P moments, with block diagonal matrix 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- LW_2004_DGP(p,nobs,ndist = 't', tdf = 9))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
      }
  
      if(dgp=='S&P Gaussian'){
        # Generating data from the S&P moments, 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- SnP_DGP(p,nobs,ndist = 'norm'))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
        evc <- snp[[6]]
      }
      
      if(dgp=='S&P Student'){
        # Generating data from the S&P moments, 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- SnP_DGP(p,nobs,ndist = 't',tdf = 9))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
        evc <- snp[[6]]
      }
      
      if(dgp=='Toeplitz'){
        # Generating data from the S&P moments, 
        time.stats[ip,'data generation',1,dgp] <- system.time(snp <- Toeplitz_DGP(p, nobs, rho = 0.5, ndist = 'norm'))[1]
        Y <- snp[[1]]
        stru <- snp[2:3]
        mu_tru <- snp[[4]]
      }
      
      if(dgp=='Sparse Inv'){
        # Generating data from sparse inverse cov matrix, 
        time.stats[ip,'data generation',1,dgp] <- system.time(spd <- Sparse_inv_DGP(p,nobs,sparsity = 0.9))[1]
        Y <- spd[[1]]
        stru <- spd[2:3]
        mu_tru <- spd[[4]]
      }
      
       if(length(grep('SPD Chol',dgp))){
        # Generating data from sparse SPD matrix, cholesky based.
        ods <- 10^-4
        if(length(grep('small od',dgp))) ods <- 10^-4
        time.stats[ip,'data generation',1,dgp] <- system.time(spd <- Sparse_SPD_DGP_chol(p,nobs,od_scale = ods,sparsity = 0.9))[1]
        Y <- spd[[1]]
        stru <- spd[2:3]
        mu_tru <- spd[[4]]
      }
          
      # Nodewise information criterion (WIC, BIC, AIC)
      if(is.null(ic)) ic <- 'WIC'
      
      # COVARIANCE MATRIX COMPUTATIONS
  
      # Comptue the sample cov mat (and inverse)    
      smp_time <- system.time(ssmp <- est_smpcov(Y))[1]
      time.stats[ip,'Sample',1,dgp] <- smp_time
      
      # Compute the covariance matrix (and inverse) with UNKNWON FACTORS 
      fac_time <- system.time(sfac <- est_faccov(Y,nfac = 3))[1]
      time.stats[ip,'Factor',1,dgp] <- fac_time
      
      # Compute the covariance matrix (and inverse) with POET 
      poet_time <- system.time(spoe <- est_poetcov(ssmp[[1]],kmax = 7,nobs,p))[1]
      time.stats[ip,'POET',1,dgp] <- poet_time
      misc.stats[ip,1:2,dgp] <- unlist(spoe[3:4])
      
      # Compute the NODEWISE inverse covariance matrix
      ndw_time <- system.time(sndw <- est_ndwcov(Y,ic))[1]
      sndw[[1]] <- ssmp[[1]] # sample covariance matrix
      time.stats[ip,'NodeWise',1,dgp] <- ndw_time
      misc.stats[ip,'Sparsity',dgp] <- sndw[[3]]
   
      # Compute the LedoitWolf inverse covariance matrix
      lw_time <- system.time(slw <- est_lwcov(Y))[1]
      time.stats[ip,'Ledoit-Wolf',1,dgp] <- lw_time
      misc.stats[ip,5:6,dgp] <- unlist(slw[3:4])
      
      # NORM and EXPOSURE COMPUTATIONS
      
      # Compute the norms for the sample 
      smp_norm <- err_norms(ssmp,stru,p)
      norm.stats[ip,,'Sample',dgp] <- smp_norm
      # Compute norms for the factor estimator 
      fac_norm <- err_norms(sfac,stru,p)
      norm.stats[ip,,'Factor',dgp] <- fac_norm
      # Compute norms for POET 
      poe_norm <- err_norms(spoe,stru,p)
      norm.stats[ip,,'POET',dgp] <- poe_norm
      # Compute the norms for the NodeWise 
      ndw_norm <- err_norms(sndw,stru,p)
      norm.stats[ip,,'NodeWise',dgp] <- ndw_norm
      misc.stats[ip,'NodeWise Lem 1',dgp] <- max(abs(sndw[[2]]%*%sndw[[1]]-diag(p)))# Computing Lemma 1 stat
      # Compute the norms for the Ledoit Wolf 
      lw_norm <- err_norms(slw,stru,p)
      norm.stats[ip,,'Ledoit-Wolf',dgp] <- lw_norm
   
      # PORTFOLIO COMPUTATIONS
      # Estimated expected return vectors
      mu_ndw <- matrix(colMeans(Y),ncol=1) 
      mu_smp <- mu_ndw
      mu_poe <- mu_ndw
      mu_lw <- mu_ndw
      mu_fac <- mu_ndw

      # true ptf
      tru_ptf <- portfolios(stru,tgvec,mu_tru,stru)
      ptf.stats[ip,,'True',dgp,'Exposure'] <- tru_ptf[[3]]
    
      # sample
      if(class(ssmp[[2]])!='try-error'){
        smp_ptf <- portfolios(ssmp,tgvec,mu_smp,stru)
        ptf.stats[ip,,'Sample',dgp,'Risk'] <- smp_ptf[[1]]
        ptf.stats[ip,,'Sample',dgp,'Risk_ratio'] <- abs(smp_ptf[[1]]/tru_ptf[[1]]-1)
        ptf.stats[ip,,'Sample',dgp,'Weights'] <- colSums(abs(smp_ptf[[2]] - tru_ptf[[2]]))
        ptf.stats[ip,,'Sample',dgp,'Exposure'] <- smp_ptf[[3]]
        ptf.stats[ip,,'Sample',dgp,'Risk_diff'] <- smp_ptf[[4]]
        time.stats[ip,'Sample',2,dgp] <- system.time(rw <- rand_ptf_expo(stru[[1]],smp_ptf[[1]],smp_ptf[[2]],expvec = expvec,rep = exprep))[1]
        expo.stats[ip,,,'Sample',dgp,'Risk'] <- rw[[1]]
        expo.stats[ip,,,'Sample',dgp,'Weights'] <- rw[[2]]
      }
      # Factor
      fac_ptf <- portfolios(sfac,tgvec,mu_fac,stru)
      ptf.stats[ip,,'Factor',dgp,'Risk'] <- fac_ptf[[1]]
      ptf.stats[ip,,'Factor',dgp,'Risk_ratio'] <- abs(fac_ptf[[1]]/tru_ptf[[1]]-1)
      ptf.stats[ip,,'Factor',dgp,'Weights'] <- colSums(abs(fac_ptf[[2]] - tru_ptf[[2]])) 
      ptf.stats[ip,,'Factor',dgp,'Exposure'] <- fac_ptf[[3]]
      ptf.stats[ip,,'Factor',dgp,'Risk_diff'] <- fac_ptf[[4]]
      time.stats[ip,'Factor',2,dgp] <- system.time(rw <- rand_ptf_expo(stru[[1]],fac_ptf[[1]],fac_ptf[[2]],expvec = expvec,rep = exprep))[1]
      expo.stats[ip,,,'Factor',dgp,'Risk'] <- rw[[1]]
      expo.stats[ip,,,'Factor',dgp,'Weights'] <- rw[[2]]
      
      # POET
      poe_ptf <- portfolios(spoe,tgvec,mu_poe,stru)
      ptf.stats[ip,,'POET',dgp,'Risk'] <- poe_ptf[[1]]
      ptf.stats[ip,,'POET',dgp,'Risk_ratio'] <- abs(poe_ptf[[1]]/tru_ptf[[1]]-1)
      ptf.stats[ip,,'POET',dgp,'Weights'] <- colSums(abs(poe_ptf[[2]] - tru_ptf[[2]])) 
      ptf.stats[ip,,'POET',dgp,'Exposure'] <- poe_ptf[[3]]
      ptf.stats[ip,,'POET',dgp,'Risk_diff'] <- poe_ptf[[4]]
      time.stats[ip,'POET',2,dgp] <- system.time(rw <- rand_ptf_expo(stru[[1]],poe_ptf[[1]],poe_ptf[[2]],expvec = expvec,rep = exprep))[1]
      expo.stats[ip,,,'POET',dgp,'Risk'] <- rw[[1]]
      expo.stats[ip,,,'POET',dgp,'Weights'] <- rw[[2]]
      # nodewise
      ndw_ptf <- portfolios(sndw,tgvec,mu_ndw,stru)
      ptf.stats[ip,,'NodeWise',dgp,'Risk'] <- ndw_ptf[[1]]
      ptf.stats[ip,,'NodeWise',dgp,'Risk_ratio'] <- abs(ndw_ptf[[1]]/tru_ptf[[1]]-1)
      ptf.stats[ip,,'NodeWise',dgp,'Weights'] <- colSums(abs(ndw_ptf[[2]] - tru_ptf[[2]]))  
      ptf.stats[ip,,'NodeWise',dgp,'Exposure'] <- ndw_ptf[[3]]
      ptf.stats[ip,,'NodeWise',dgp,'Risk_diff'] <- ndw_ptf[[4]]
      time.stats[ip,'NodeWise',2,dgp] <- system.time(rw <- rand_ptf_expo(stru[[1]],ndw_ptf[[1]],ndw_ptf[[2]],expvec = expvec,rep = exprep))[1]
      expo.stats[ip,,,'NodeWise',dgp,'Risk'] <- rw[[1]]
      expo.stats[ip,,,'NodeWise',dgp,'Weights'] <- rw[[2]]
      
      # ledoitwolf
      lw_ptf <- portfolios(slw,tgvec,mu_lw,stru)
      ptf.stats[ip,,'Ledoit-Wolf',dgp,'Risk'] <- lw_ptf[[1]]
      ptf.stats[ip,,'Ledoit-Wolf',dgp,'Risk_ratio'] <- abs(lw_ptf[[1]]/tru_ptf[[1]]-1)
      ptf.stats[ip,,'Ledoit-Wolf',dgp,'Weights'] <- colSums(abs(lw_ptf[[2]] - tru_ptf[[2]]))  
      ptf.stats[ip,,'Ledoit-Wolf',dgp,'Exposure'] <- lw_ptf[[3]]
      ptf.stats[ip,,'Ledoit-Wolf',dgp,'Risk_diff'] <- lw_ptf[[4]]
      time.stats[ip,'Ledoit-Wolf',2,dgp] <- system.time(rw <- rand_ptf_expo(stru[[1]],lw_ptf[[1]],lw_ptf[[2]],expvec = expvec,rep = exprep))[1]
      expo.stats[ip,,,'Ledoit-Wolf',dgp,'Risk'] <- rw[[1]]
      expo.stats[ip,,,'Ledoit-Wolf',dgp,'Weights'] <- rw[[2]]
      
    }
    
    # End of p-loop, update counter
    ip <- ip + 1
  }
  
  titer <- (proc.time() - tstart)[3]
  
  cat(paste0('\n Iteration ',i,' completed in ',round(titer,0),' sec.'))
  
  return(list('norms' = norm.stats,'portfolios' = ptf.stats, 'times' = time.stats,
              'titer' = titer, 'misc' = misc.stats, 'exposures' = expo.stats))
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

# A function that takes the list returned by mclapply and bind the arrays for post processing 
postmc <- function(mc){
  
  # bind the arrays
  norm.list <- abind(lapply(mc,'[[','norms'),along = 5)
  expo.list <- abind(lapply(mc,'[[','exposures'),along = 7)
  ptf.list <- abind(lapply(mc,'[[','portfolios'),along = 6)
  time.list <- abind(lapply(mc,'[[','times'),along = 5)
  titer <- sapply(mc,'[[','titer')
  misc.list <- abind(lapply(mc,'[[','misc'),along = 4)
    
    
  # compute avg error and error std dev for the norms
  # keeping the dimension names is a pain in the ...
  norm.stats <- array(NA,dim = c(dim(norm.list)[1:4],2),
                      dimnames = c('p'=dimnames(norm.list)[1],
                                   'Norm'=dimnames(norm.list)[2],
                                   'Estimator'=dimnames(norm.list)[3],
                                   'dgp' = dimnames(norm.list)[4],
                                   'stats' = list(c('Mean','Standard deviation'))))
  norm.stats[,,,,1] <- apply(norm.list,1:4,mean)
  norm.stats[,,,,2] <- apply(norm.list,1:4,se)
  
  # Aggregation for the exposures
  expo.stats <- array(NA,dim = c(dim(expo.list)[1:6],2),
                      dimnames = c('p'=dimnames(expo.list)[1],
                                   'Portfolio'=dimnames(ptf.list)[2],
                                   'Exposure'=dimnames(expo.list)[3],
                                   'Estimator'=dimnames(expo.list)[4],
                                   'dgp' = dimnames(expo.list)[5],
                                   'measure' = dimnames(expo.list)[6],
                                   'stats' = list(c('Median','Standard deviation'))))
  expo.stats[,,,,,,1] <- apply(expo.list,1:6,median)
  expo.stats[,,,,,,2] <- apply(expo.list,1:6,se)
  
  
  # compute stats for portfolio variance estimation error
  ptf.stats <- array(NA,dim = c(dim(ptf.list)[1:5],3),
                     dimnames = c('p'=dimnames(ptf.list)[1],
                                  'Portfolio'=dimnames(ptf.list)[2],
                                  'Estimator'=dimnames(ptf.list)[3],
                                  'dgp' = dimnames(ptf.list)[4],
                                  'measure' = dimnames(ptf.list)[5],
                                  'stats' = list(c('Root mean square','Mean','Median'))))
  
  
  ptf.stats[,,,,,1] <- apply(ptf.list,1:5,rmse)
  ptf.stats[,,,,,2] <- apply(ptf.list,1:5,mnae)
  ptf.stats[,,,,,3] <- apply(ptf.list,1:5,mdae)
  
  # compute average times 
  time.stats <- apply(time.list,1:4,mean)
  names(attributes(time.stats)$dimnames) <- c('p','Estimator','Computation','dgp')
  
  
  # compute average iter time
  niter <- length(titer)
  time.iter <- mean(titer)
  

  # compute average poet stats
  misc.stats <- apply(misc.list,1:3,mean)
  names(attributes(misc.stats)$dimnames) <- c('p','Parameter','dgp')
  
  return(list('norms'=norm.stats,'ptf'=ptf.stats,'time'=time.stats
              ,'iter'=time.iter,'misc'=misc.stats,'exposure'=expo.stats,'niter'=niter))
}






# Nodewise estimation of the covariance matrix
est_ndwcov <- function(Y,ic){
  
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  C <- matrix(0,p,p)
  diag(C) <- 1
  tau <- NULL
  
  # Loop over the assets
  for(j in 1:p){
    # Estimate the Lasso
    jlas <- glmnet(x=Y[,-j],y=Y[,j],family = 'gaussian')
    # Get fit
    jfit <- predict(jlas, newx=Y[,-j], type="response")    
    # residuals
    jres <- matrix(Y[,j],n,length(jlas$lambda)) - jfit
    # std err
    jsig <- colSums(jres^2)/n
    # Computing information criterion
    if(ic=='WIC') jbic  <- log(jsig) + jlas$df * log(n)/n * log(log(p)) # BIC (Wang,2010)
    if(ic=='BIC') jbic  <- log(jsig) + jlas$df * log(n)/n  #BIC
    if(ic=='GIC') jbic  <- log(jsig) + jlas$df * log(p) * log(log(n))/n # Fan & Tang JRSS-B 2004
    if(ic=='AIC') jbic  <- log(jsig) + 2 * jlas$df # AIC 
    # Index of selected model 
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[,jind]
    # Computing tau squared. Two formulas in text
    jtau <- sum(jres[,jind]^2)/n + jlas$lambda[jind]*sum(abs(jpar)) # using (12)
    
    # using the msgps package
    # Not used because it's very slow!
    #jlas <- msgps(Y[,-j],as.vector(Y[,j]))
    #jpar <- jlas$dfbic_result$coef[-1]
    #jtun <- jlas$dfbic_result$tuning
    #jres <- Y[,j]-predict(jlas,Y[,-j],tuning = jlas$dfbic_result$tuning)
    #jtau <- sum(jres^2) + jtun*sum(abs(jpar))

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



# factor covariance matrix estimation with inversion
est_faccov <- function(Y,nfac = 3){
  
  # Estimating factors
  fac <- prcomp(Y)$x[,1:nfac]
  
  # OLS estimation
  fmod <- lm(Y~fac)
  bhat <- t(fmod$coefficients[-1,])
  vres <- apply(fmod$residuals,2,var)
  
  # cov mat est
  sigfac <- bhat %*% var(fac) %*% t(bhat) + diag(vres)
  
  # inversion (several options)
  #sigfacinv <- solve(sigfac)
  sigfacinv <- chol2inv(chol(sigfac))
  
  return(list(sigfac,sigfacinv,bhat,fac))
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
est_poetcov <- function(ssmp,kmax,nobs,p,Khat=NULL){
  
  # Eigen value/vector decomposition
  ev <- eigen(ssmp)
  eval <- ev$values
  evec <- ev$vectors
  
  if(is.null(Khat)){
    # Determining K-hat
    evsum_back <- sum_ev_back(eval,evec,kmax)
    Kvec <- (1/p) * sapply(evsum_back,tr) + (0:kmax)*(p+nobs)*log((p*nobs)/(p+nobs))/(p*nobs)
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
    Cthres <- Clev*(sqrt(log(p)/nobs)+sqrt(1/p))
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


# Compute the covariance matrix estimation error norms 
err_norms <- function(estsigma, trusigma, p){
  
  # Estimation error
  serr <- estsigma[[1]] - trusigma[[1]]
  if(class(estsigma[[2]])!='try-error') {
    serr_inv <- estsigma[[2]] - trusigma[[2]]
    }
  else {
    serr_inv <- NA
  }
   
  # Frobenius norms
  frob <- norm(serr,type = 'F')
  frob_inv <- ifelse(length(serr_inv)==1,NA,norm(serr_inv,type = 'F'))
  
  # est * inverse true 
  sprod <- estsigma[[1]] %*% trusigma[[2]]
  # Entropy norm
  entropy <- tr(sprod) - log(det(sprod)) - p
  # Sigma norm
  signorm <- (1/sqrt(p)) * tr(sprod - diag(p))^2
  # MnAE
  mnae <- mean(abs(serr))
  # MdAE
  mdae <- median(abs(serr))
  
  return(c(frob,signorm,entropy,mnae,mdae,frob_inv))
}


# Compute the portfolio weights and variance
# target: a scalar with the target return
# mu: column vector of expected returns 
portfolios <- function(estsigm,tgvec,mu,trusigma){
  
  # Compute preliminary values
  p <- length(mu)
  one <- matrix(1,p,1)
  phi <- as.single(t(one) %*% estsigm[[2]] %*% one)
  psi <- as.single(t(one) %*% estsigm[[2]] %*% mu)
  Phi <- as.single(t(mu) %*% estsigm[[2]] %*% mu)
  
  # init storage
  mkwtg <- NULL
  mkwvartg <- NULL
  mkwexptg <- NULL
  mkwrsktg <- NULL
  # looping over the targets
  for(target in tgvec){
    # Markowitz 
    denom <- phi*Phi - psi^2
    mkw <- (Phi-target*psi)/denom * estsigm[[2]] %*% one + (target*phi-psi)/denom * estsigm[[2]]%*%mu 
    mkwvar <- (phi*target^2 - 2*psi*target + Phi)/denom
    rsk <- t(mkw) %*% (estsigm[[1]] - trusigma[[1]]) %*% mkw
    
    # storing 
    mkwtg <- cbind(mkwtg,mkw)
    mkwvartg <- c(mkwvartg,mkwvar)
    mkwexptg <- c(mkwexptg,sum(abs(mkw)))
    mkwrsktg <- c(mkwrsktg,rsk)
  }
    
  # Global min variance
  gmv <- (estsigm[[2]] %*% one) / phi 
  gmvvar <- 1 / phi
  gmvexp <- sum(abs(gmv))
  gmvrsk <- t(gmv) %*% (estsigm[[1]] - trusigma[[1]]) %*% gmv 

  return(list(c(mkwvartg,gmvvar),cbind(mkwtg,gmv),c(mkwexptg,gmvexp),c(mkwrsktg,gmvrsk)))
}


# Compute the risk of a portfolio with random weights and fixed exposure.
rand_ptf_expo <- function(tru_sig,est_rsk,est_w,expvec,rep){
  
  # dimensions
  p <- ncol(tru_sig)
  # prep storage
  exprsk <- matrix(NA,nrow=length(est_rsk),ncol=length(expvec))
  expw <- matrix(NA,nrow=length(est_rsk),ncol=length(expvec))
  
  # looping over the estimated portfolios
  for(cp in 1:length(est_rsk)){
    # counter
    expcount <- 1
    # looping over the exposures
    for(expo in expvec){
      wrsk <- c()
      ww <- c()
      for(r in 1:rep){
        # generate weights using FFL 2015 method
        w <- mk_ptf_weights(p,expo)
        # compute w delta
        dw <- sum(abs(w-est_w[,cp]))
        # compute risk delta
        drsk <- (est_rsk[cp]-t(w)%*%tru_sig%*%w)
        # aggregate
        wrsk <- c(wrsk,drsk)
        ww <- c(ww,dw)
      }
      #storing
      exprsk[cp,expcount] <- mean(wrsk)
      expw[cp,expcount] <- mean(ww)
      # update counter
      expcount <- expcount + 1
    }
  }
  # return
  return(list(exprsk,expw))
}


# Generate portfolio weights for p assets and expo exposure (l1 norm)
mk_ptf_weights <- function(p,expo){
  # storage
  w <- rep(0,p)
  # draw from a binomial, ensure k<p
  k <- p
  if(expo>1)while(k==p)  k <- rbinom(1,p,(expo+1)/(2*expo))
  # indices for pos and neg weights
  rp <- sample(1:p,k)
  rn <- (1:p)[-rp]
  # positive weights
  ep <- rexp(k,rate = 1)
  wpos <- (1+expo)*ep/sum(2*ep)
  # negative weights
  en <- rexp(p-k,rate = 1)
  wneg <- (1-expo)*en/sum(2*en)
  # construct the weight vector
  w[rp] <- wpos
  w[rn] <- wneg
  
  return(matrix(w,ncol=1))
}
