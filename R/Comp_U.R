# the first derivative of SCAD function
Dpfunc = function(u,lambda,a)
{
  if(u<=lambda) Dpval = lambda
  else if(u<a*lambda) Dpval = -(u-a*lambda)/(a-1)
  else Dpval = 0
  Dpval
}

# Compute TT_j
slos.compute.weights = function(basis){
  L       = basis$nbasis
  rng     = fda::getbasisrange(basis)
  breaks  = c(rng[1],basis$params,rng[2])
  M       = length(breaks) - 1
  norder  = L-M+1
  TT       = array(0,dim=c(L,L,M))

  for (j in 1:M)
  {
    temp = fda::inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    # TT[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
    TT[,,j] <- temp
  }
  TT
}

# Compute U
Comp_U <- function(basis, gamma_est, lambda, absTol){

  L <- basis$nbasis

  if(lambda == 0){
    U <- matrix(0, ncol = L, nrow = L)
    id <- 1:L
  }else{

    rng <- fda::getbasisrange(basis)
    length_t <- rng[2] - rng[1]
    breaks <- c(rng[1],basis$params,rng[2])
    M <- length(breaks) - 1
    d <- L - M
    TT <- slos.compute.weights(basis)

    bZeroMat <- rep(FALSE,L)
    bZeroMat[c(1, L)] <- TRUE

    U <- matrix(0, ncol = L, nrow = L)
    beta_norm <- NULL
    for(j in 1:M){
      beta_norm[j] <- sqrt(t(gamma_est) %*% TT[,,j] %*% gamma_est)
      if(beta_norm[j] <= absTol){
        bZeroMat[j:(j + d)] <- TRUE
      }else{
        Uj <- sqrt(M/length_t) * Dpfunc(sqrt(M/length_t) * beta_norm[j], lambda = lambda,
                                        a = 3.7)/beta_norm[j] * TT[,,j]
        U <- U + Uj
      }

    }

    U <- U/2
    id <- which(bZeroMat == FALSE)

  }

  res <- list()
  res$U <- U
  res$id <- id

  return(res)

}


