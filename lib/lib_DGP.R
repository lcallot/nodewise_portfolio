# A set of functions to generate data for the large covariance matrix estimation project. 
# Laurent Callot (l.callot@gmail.com), 29/04/2016

# Data generating process from Fan, Fan, Liao 2008 JoE. 
# The moments and other values are sourced from constant_FFL.R
FFL_DGP <- function(fac,p){
  
  # number of obs
  nobs <- nrow(fac)
  
  # loading the DGP moments
  source('../lib/constants_FFL.R')
  
  # generate the factor loadings
  bfac <- mvrnorm(n=p,mu=mub,Sigma=covb)
  
  # generate the erros std dev and lower bounding
  sdev <- rgamma(p,shape=alpha,rate=1/beta)
  sdev[sdev<0.1950] <- 0.1950
  
  # generate the errors
  errs <- mvrnorm(n = nobs, mu = rep(0,p), Sigma = diag(sdev^2))
  
  # generate the data    
  Y <- fac %*% t(bfac) + errs
  
  # cov mat est
  sigtru <- bfac %*% var(fac) %*% t(bfac) + diag(sdev^2)
  # inversion 
  sigtruinv <- chol2inv(chol(sigtru))
  
  # expected return vector
  mu_tru <- bfac %*% matrix(colMeans(fac),ncol=1)
  
  return(list(Y,sigtru,sigtruinv,mu_tru))
}


SnP_DGP <- function(p,nobs,ndist='norm',tdf=9,sparsity = 1){
  
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  
  not.ok <- TRUE
  ntry <- 0
  evc <- 0
  while(not.ok){
    # create the grid of matrix entries
    gr <- expand.grid(1:p,1:p)
    gr <- as.matrix(gr[gr[,1]<gr[,2],])
    # selecting non-zero entries
    nz <- sample(1:nrow(gr),floor(sparsity*nrow(gr)))
    # enforcing symmetry
    fgr <- rbind(gr[nz,],gr[nz,2:1])
    
    # Generating the standard deviations
    sdev <- 100 * rgamma(p,shape=ver['shape'],rate=ver['rate'])
    # Generating the mean
    mu_tru <- rnorm(p,mean = mer['mean'],sd=mer['sd'])
    # Generating the off-diagonal elements
    oder <- 100 * rnorm(length(nz),mean = oer['mean'],sd=oer['sd'])
    
  
    # creating the matrix
    sigtru <- Diagonal(x = sdev^2)
    # populating off diagonal
    sigtru[fgr] <- rep(oder,2)

    evs <- eigen(sigtru)
    ev_threshold <- 1e-10
    # eigenvalue cleaning if some values too low
    if(sum(evs$values<ev_threshold)>0){
      evc <- 1
      minev <- min(evs$values[evs$values>ev_threshold])
      evs$values[evs$values<ev_threshold] <- minev
      sigtru <- t(evs$vectors)%*%diag(evs$values)%*%evs$vectors
    }
      
    # checking invertibility
    trychol <- try(chol(sigtru),silent = TRUE)
    not.ok <- (class(trychol)=='try-error')
    if(!not.ok) sigtruinv <- chol2inv(trychol)
    ntry <- ntry + 1
    #if(not.ok) {print(eigen(sigtru)$values)}
    if(ntry>99) break('Too many attemps at generating and invertible covariance matrix')
  }
  
  # de-sparsifying
  sigtru <- as.matrix(sigtru) 
  sigtruinv <- as.matrix(sigtruinv)
 
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sigtru)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sigtru, df= tdf)
  
  return(list(Y,sigtru,sigtruinv,mu_tru,ntry,evc))
}



Sparse_SPD_DGP_chol <- function(p,nobs,sparsity=0.8,od_scale=0.1,ndist='norm',tdf=9){
  
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  # Generating the mean of the return vectors
  mu_tru <- rnorm(p,mean = mer['mean'],sd=mer['sd'])
  
  # create the grid of matrix entries
  gr <- expand.grid(1:p,1:p)
  # taking the upper diagonal
  gr <- as.matrix(gr[gr[,1]<gr[,2],])
  # selecting non-zero entries
  nz <- sample(1:nrow(gr),floor((1-sparsity)*nrow(gr)))
  # generating the cholesky factors 
  cfac <- matrix(0,p,p)
  cfac[nz] <- od_scale*runif(length(nz),-1,1)
  diag(cfac) <- runif(p,0,1)
  # creating the cov matrix
  sigtru <- cfac%*%t(cfac)
  
  # the inverse covariance matrix
  sigtruinv <- chol2inv(chol(sigtru))
  
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sigtru)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sigtru, df= tdf)
  
  return(list(Y,sigtru,sigtruinv,mu_tru))
}


Sparse_inv_DGP <- function(p,nobs,sparsity=0.8,ndist='norm',tdf=9){
  
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  # Generating the mean of the return vectors
  mu_tru <- 100*rnorm(p,mean = mer['mean'],sd=mer['sd'])
  
  # random state for pyhton
  rs <- sample(1:100000, 1)
  
  # Generating sparse sym positive definite matrix with sklearn make_sparse_spd_matrix
  # This will be the inverse covariance matrix
  python.load('../lib/sparse_spd.py')
  sigtruinv <- python.call('mk_spd',dim=p,alpha=sparsity,maxc=0.009,minc=0.001,rs=rs)
  sigtruinv <- 100*do.call(cbind,sigtruinv)
  # the covariance matrix
  sigtru <- chol2inv(chol(sigtruinv))
  
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sigtru)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sigtru, df= tdf)
  
  return(list(Y,sigtru,sigtruinv,mu_tru))
}



# Sparse block diag DGP
# A fraction 1-sparsity of the matrix is a non-sparse SPD matrix, the remainder is 0 
Sparse_block_diag_DGP <- function(p,nobs,sparsity = 0.5, ndist = 'norm',tdf = 9){
  
  # number of blocks computed to ensure nbr non-zero correlations per entry consistent with sparsity   
  nbr_blocks <- round(1/(1-sparsity),0)
  size_block <- floor(p/nbr_blocks)
  # to be filled...
  sparse_cov <- matrix(0,p,p)
  
  for(i in 1:nbr_blocks){
    if(i==nbr_blocks) size_block <- p-(nbr_blocks-1)*size_block # ensures the last block fills the matrix
    # index for block i in full matrix
    if(i<nbr_blocks) block_ind <-  ((i-1)*size_block+1):(i*size_block) 
    if(i==nbr_blocks) block_ind <-  (p-size_block+1):p 
    # generate block i
    non_sparse_block <- SnP_DGP(size_block,nobs,ndist=ndist,tdf=tdf,sparsity = 1)[[2]]
    sparse_cov[block_ind,block_ind] <- non_sparse_block
  }
  
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  # Generating the mean of the return vectors
  mu_tru <- rnorm(p,mean = mer['mean'],sd=mer['sd'])
  
  # inverse
  # no need to try, the blocks are checked for invertibility so the block-diag mat will be invertible
  sigtruinv <- solve(sparse_cov)
  
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sparse_cov)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sparse_cov, df= tdf)
  
  return(list(Y,sparse_cov,sigtruinv,mu_tru))
}


LW_2004_DGP <- function(p, nobs, ndist = 'norm',tdf = 9){
  
  # draw the variances
  dg <- rlnorm(p, meanlog = 0, sdlog = 1)
  # generate the covariance matrix
  sigtru <- diag(dg)
  # generate the true inverse
  sigtruinv <- diag(1/dg)
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  # Generating the mean of the return vectors
  mu_tru <- rnorm(p,mean = mer['mean'],sd=mer['sd'])
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sigtru)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sigtru, df= tdf)

  return(list(Y,sigtru,sigtruinv,mu_tru))
}


Toeplitz_DGP <- function(p, nobs, rho = 0.5, ndist = 'norm', tdf = 9){
  
  # loading the DGP moments
  source('../lib/constants_SnP.R')
  
  # Generating the standard deviation vector
  sdev <- rgamma(p,shape=ver['shape'],rate=ver['rate'])
  smat <- matrix(sdev,p ,1) %*% matrix(sdev, 1, p)
  
  # Generating Toeplitz matrix
  trow <- rho^seq.int(0, p-1)
  tmat <- toeplitz(trow)
  
  # Generating the covariance matrix 
  sigtru <- smat * tmat
  # sigtru <- tmat

  # the inverse covariance matrix
  sigtruinv <- chol2inv(chol(sigtru))
  
  # Generating the means
  mu_tru <- rnorm(p,mean = mer['mean'],sd=mer['sd'])
  
  # generating the data 
  if(ndist=='norm') Y <- mvrnorm(n = nobs, mu = mu_tru, Sigma = sigtru)
  if(ndist=='t') Y <- mu_tru + rmvt(n = nobs,sigma = sigtru, df= tdf)
  
  
  return(list(Y,sigtru,sigtruinv,mu_tru))  
}