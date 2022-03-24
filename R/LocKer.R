#' @title Locally sparse estimator of generalized varying coefficient model for asynchronous longitudinal data.
#'
#' @description Locally sparse estimator of generalized varying coefficient model for asynchronous longitudinal data by kernel-weighted estimating equation. The function is suitable for generalized varying coefficient model with one covariate.
#'
#' @param X A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of the covariate for each subject at the observation time correspond to \code{X_obser}.
#' @param Y A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of the response for each subject at the observation time correspond to \code{Y_obser}.
#' @param family A \code{character} string representing the distribution family of the response. The value can be "Gaussian", "binomial", "poisson".
#' @param X_obser_num A \code{vector} denoting the observation size of the covariate for each subject.
#' @param Y_obser_num A \code{vector} denoting the observation size of the response for each subject.
#' @param X_obser A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation times of the covariate for each subject.
#' @param Y_obser A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation times of the response for each subject.
#' @param timeint A \code{vector} of length two denoting the supporting interval.
#' @param L_list A \code{vector} denoting the candidates for the number of B-spline basis functions. The best \code{L} is chosen by cross-validation.
#' @param roupen_para_list A \code{vector} denoting the candidates for the roughness parameters. The best roughness parameter is chosen by \code{EBIC} together with sparseness parameter.
#' @param lambda_list A \code{vector} denoting the candidates for the sparseness parameter. The best sparseness parameter is chosen by \code{EBIC} together with roughness parameter.
#' @param absTol_list A \code{vector} denoting the threshold of the norm for coefficient function on each sub-interval. The \code{vector} is related to \code{L_list}, with the same length as \code{L_list}.
#' @param nfold An \code{integer} denoting the number of fold for the selection of \code{L} by cross-validation. (default: 5)
#' @param d An \code{integer} denoting the degree of B-spline basis functions. (default: 3)
#'
#' @return A \code{list} containing the following components:
#' \item{beta0fd_est}{A functional data object denoting the estimated intercept function.}
#' \item{betafd_est}{A functional data object denoting the estimated coefficient function.}
#' \item{time}{A \code{scalar} denoting the computation time.}
#' \item{L}{An \code{integer} denoting the selected number of B-spline basis function.}
#' \item{roupen_select}{A \code{scalar} denoting the selected roughness parameter.}
#' \item{lambda_select}{A \code{scalar} denoting the selected sparseness parameter.}
#' \item{EBIC}{A \code{matrix} denoting the \code{EBIC} scores for various roughness parameters and sparseness parameters belongs to the candidates when using the selected \code{L}.}
#' @export
#'
#' @examples
#' ####Generate data
#' n <- 200
#' beta0 <- function(x){cos(2 * pi * x)}
#' beta <- function(x){sin(2 * pi * x)}
#' Y_rate <- 15
#' X_rate <- 15
#' Y_obser_num <- NULL
#' X_obser_num <- NULL
#' Y_obser <- list()
#' X_obser <- list()
#' for(i in 1:n){
#' Y_obser_num[i] <- stats::rpois(1, Y_rate) + 1
#' Y_obser[[i]] <- stats::runif(Y_obser_num[i], 0, 1)
#' X_obser_num[i] <- stats::rpois(1, X_rate) + 1
#' X_obser[[i]] <- stats::runif(X_obser_num[i], 0, 1)
#' }
#' ## The covariate functions Xi(t)
#' X_basis <- fda::create.bspline.basis(c(0, 1), nbasis = 74, norder = 5,
#' breaks = seq(0, 1, length.out = 71))
#' a <- matrix(0, nrow = n, ncol = 74)
#' X <- list()
#' XY <- list() #X at the observation time of Y
#' muY <- list()
#' for(i in 1:n){
#' a[i,] <- stats::rnorm(74)
#' Xi_B <- splines::bs(X_obser[[i]], knots = seq(0, 1, length.out = 71)[-c(1, 71)],
#' degree = 4, intercept = TRUE)
#' X[[i]] <- Xi_B %*% a[i,]
#' Yi_B <- splines::bs(Y_obser[[i]], knots = seq(0, 1, length.out = 71)[-c(1, 71)],
#' degree = 4, intercept = TRUE)
#' XY[[i]] <- Yi_B %*% a[i,]
#' muY[[i]] <- beta0(Y_obser[[i]]) + XY[[i]] * beta(Y_obser[[i]])
#' }
#' Y <- list()
#' errY <- list()
#' for(i in 1:n){
#' errY[[i]] <- stats::rnorm(Y_obser_num[[i]], mean = 0, sd = 1)
#' Y[[i]] <- muY[[i]] + errY[[i]]
#' }
#' L_list <- 20
#' absTol_list <- 10^(-3)
#' roupen_para_list <- 1.5 * 10^(-3)
#' lambda_list <- c(0, 0.001, 0.002)
#' LocKer_list <- LocKer(X, Y, family = "Gaussian", X_obser_num, Y_obser_num, X_obser,
#' Y_obser, timeint = c(0, 1), L_list, roupen_para_list, lambda_list, absTol_list)
LocKer <- function(X, Y, family, X_obser_num, Y_obser_num, X_obser, Y_obser, timeint,
                   L_list, roupen_para_list, lambda_list, absTol_list, nfold = 5, d = 3){

  if(family == "Gaussian"){
    linkfun <- function(x){x}
    linkinv <- function(x){x}
    bfun <- function(x){x^2/2}
    linkinvder <- function(x){x - x + 1}
  }else if(family == "binomial"){
    linkfun <- function(x){(1 + exp(-x))^(-1)}
    linkinv <- function(x){log(x/(1 - x))}
    bfun <- function(x){log(1 + exp(x))}
    linkinvder <- function(x){1/(x * (1 -x))}
  }else if(family == "poisson"){
    linkfun <- function(x){exp(x)}
    linkinv <- function(x){log(x)}
    bfun <- function(x){exp(x)}
    linkinvder <- function(x){1/x}
  }else{
    print("The setting of family is not appropriate!")
  }

  kernfun <- function(x){
    max(0, 0.75 * (1 - x^2))
  }

  start_time <- Sys.time()

  n <- length(X)
  nl <- length(L_list)
  if(nl == 1){
    L_id <- 1
  }else{
    #######CV############
    #subsample
    sample_id <- rep(nfold, n)
    nid_cand <- 1:n
    n_cv <- round(n/nfold)
    for(i in 1:(nfold - 1)){

      id_select <- sample(nid_cand, size = n_cv, replace = F)
      sample_id[id_select] <- i

      nid_cand <- nid_cand[-which(nid_cand %in% id_select)]

    }

    #comoute CV score
    CV_score <- NULL # length = nl
    for(lid in 1:nl){

      L <- L_list[lid]
      absTol <- absTol_list[lid]

      CV_score_ind <- NULL # length = nfold

      for(i in 1:nfold){

        pre_id <- which(sample_id == i)
        gamma_cv_list <- LocKer_ind(X[-pre_id], Y[-pre_id], family, X_obser_num[-pre_id],
                                    Y_obser_num[-pre_id], X_obser[-pre_id], Y_obser[-pre_id],
                                    timeint, L, d, roupen_para_list,lambda_list, absTol,
                                    kernfun, linkfun, linkinvder)
        CV_score_ind[i] <- CV_ind(X[pre_id], Y[pre_id], family, X_obser[pre_id],
                                  Y_obser[pre_id], gamma_cv_list$gamma_est, timeint, L,
                                  d, gamma_cv_list$h, kernfun, linkfun)

      }

      CV_score[lid] <- mean(CV_score_ind)


    }

    L_id <- which.min(CV_score)
  }

  L <- L_list[L_id]
  absTol <- absTol_list[L_id]

  beta_basis <- fda::create.bspline.basis(timeint, nbasis = L, norder = d + 1,
                                     breaks = seq(timeint[1], timeint[2], length.out = L - d + 1))

  gamma_cv_list <- LocKer_ind(X, Y, family, X_obser_num, Y_obser_num, X_obser, Y_obser,
                              timeint, L, d, roupen_para_list,lambda_list, absTol, kernfun,
                              linkfun, linkinvder)

  B_est <- splines::bs(seq(timeint[1], timeint[2], length.out = 100), knots = beta_basis$params,
              degree = beta_basis$nbasis - length(beta_basis$params) - 1, intercept = T)
  beta_est <- B_est %*% gamma_cv_list$gamma_est[(L + 1):(2 * L)]
  beta_estfd <- fda::smooth.basis(seq(timeint[1], timeint[2], length.out = 100), beta_est, beta_basis)$fd
  beta0_est <- B_est %*% gamma_cv_list$gamma_est[1:L]
  beta0_estfd <- fda::smooth.basis(seq(timeint[1], timeint[2], length.out = 100), beta0_est, beta_basis)$fd

  end_time <- Sys.time()

  time_locker <- difftime(end_time, start_time, units = "sec")

  results <- list()
  results$beta0fd_est <- beta0_estfd
  results$betafd_est <- beta_estfd
  results$time <- time_locker
  results$L <- L
  results$roupen_select <- gamma_cv_list$roupen_select
  results$lambda_select <- gamma_cv_list$lambda_select
  results$EBIC <- gamma_cv_list$EBIC

  return(results)

}

EBIC <- function(X_tilde, Y, family, gamma_est, W, XWX, B_der, roupen_para, n_real,
                 Y_obser_num, linkfun){

  n <- length(Y)

  # sse
  if(family == "Gaussian"){
    SSE <- NULL
    for(i in 1:n){
      for(j in 1:Y_obser_num[i]){
        Y_hat_ik <- linkfun(X_tilde[[i]] %*% gamma_est)
        resi_ij <- (Y[[i]][j] - Y_hat_ik)^2/2
        SSE_ij <- 2 * t(resi_ij) %*% diag(W[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }else if(family == "binomial"){
    SSE <- NULL
    for(i in 1:n){
      for(j in 1:Y_obser_num[i]){
        Y_hat_ik <- linkfun(X_tilde[[i]] %*% gamma_est)
        if(Y[[i]][j] == 1){
          resi_ij <- Y[[i]][j] * log(Y[[i]][j]/Y_hat_ik)
        }else{
          resi_ij <- (1 - Y[[i]][j]) * log((1 - Y[[i]][j])/(1 - Y_hat_ik))
        }
        SSE_ij <- 2 * t(resi_ij) %*% diag(W[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }else if(family == "poisson"){
    SSE <- NULL
    for(i in 1:n){
      for(j in 1:Y_obser_num[i]){
        Y_hat_ik <- linkfun(X_tilde[[i]] %*% gamma_est)
        resi_ij <- Y_hat_ik - Y[[i]][j] * log(Y_hat_ik)
        SSE_ij <- 2 * t(resi_ij) %*% diag(W[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }

  # nonzero_num <- length(which(gamma_est != 0))
  # df
  id <- which(gamma_est != 0)
  XWX_sub <- XWX[id, id]
  B_der_sub <- B_der[id, id]
  df <- psych::tr(solve(XWX_sub + n * roupen_para * B_der_sub) %*% XWX_sub)
  nonzero_num <- df

  EBIC_score <- log(sum(SSE)) + nonzero_num * log(n_real)/n_real +
    0.5 * nonzero_num * log(length(gamma_est))/n_real

  EBIC_score

}

LocKer_ind <- function(X, Y, family, X_obser_num, Y_obser_num, X_obser, Y_obser, timeint,
                       L, d, roupen_para_list, lambda_list, absTol, kernfun, linkfun,
                       linkinvder){

  n <- length(X)
  nLM <- sum(X_obser_num * Y_obser_num)
  beta_basis <- fda::create.bspline.basis(timeint, nbasis = L, norder = d + 1,
                                     breaks = seq(timeint[1], timeint[2], length.out = L - d + 1))
  B_der <- fda::bsplinepen(beta_basis)
  p <- L
  V_der <- as.matrix(Matrix::bdiag(B_der, B_der))

  timecom <- NULL
  timecom_ind <- list()
  timediff_ind <- list()
  for(i in 1:n){
    timecom_ind[[i]] <- expand.grid(X_obser[[i]], Y_obser[[i]])
    timecom <- rbind(timecom, timecom_ind[[i]])
    timediff_ind[[i]] <- timecom_ind[[i]][,1] - timecom_ind[[i]][,2]
  }
  timediff <- timecom[,1] - timecom[,2]

  h <- max(stats::quantile(unlist(lapply(timediff_ind, function(x){min(abs(x))})), probs = 0.95),
           0.01)

  X_tilde <- list()
  beta_B <- list()
  for(i in 1:n){
    beta_B[[i]] <- splines::bs(X_obser[[i]],
                      knots = seq(timeint[1], timeint[2], length.out = L - d + 1)[-c(1, L - d + 1)],
                      degree = d, intercept = T)
    X_tilde[[i]] <- cbind(beta_B[[i]], as.vector(X[[i]]) * beta_B[[i]])
  }

  XWX <- matrix(0, ncol = 2 * p, nrow = 2 * p)
  XW <- list()
  omega <- NULL
  W <- list()
  for(i in 1:n){
    XW[[i]] <- list()
    W[[i]] <- list()
    for(j in 1:Y_obser_num[i]){
      timediff_ij <- timediff_ind[[i]][which(timecom_ind[[i]][,2] == Y_obser[[i]][j])]
      omega_ij <- sapply(timediff_ij/h, kernfun)/h
      omega <- c(omega, omega_ij)
      W[[i]][[j]] <- diag(omega_ij, ncol = length(timediff_ij), nrow = length(timediff_ij))
      XW[[i]][[j]] <- t(X_tilde[[i]]) %*% W[[i]][[j]]
      XWX <- XWX + XW[[i]][[j]] %*% X_tilde[[i]]
    }
  }
  n_real <- length(which(omega != 0))

  nii <- length(roupen_para_list)
  njj <- length(lambda_list)
  EBIC_score <- matrix(Inf, ncol = njj, nrow = nii)
  gamma_res <- list()
  for(ii in 1:nii){
    gamma_res[[ii]] <- list()
    for(jj in 1:njj){

      roupen_para <- roupen_para_list[ii]
      lambda <- lambda_list[jj]

      #####iteration
      # initial value
      # gamma_ini <- matrix(1, nr = p, nc = 1)
      Z <- list()
      for(i in 1:n){
        Z[[i]] <- list()
        for(j in 1:Y_obser_num[i]){
          Z[[i]][[j]] <- rep(Y[[i]][j], X_obser_num[i])
        }
      }
      XWZ <- matrix(0, ncol = 1, nrow = 2 * p)
      for(i in 1:n){
        for(j in 1:Y_obser_num[[i]]){
          XWZ <- XWZ + XW[[i]][[j]] %*% Z[[i]][[j]]
        }
      }
      gamma_ini <- solve(XWX + n * roupen_para * V_der) %*% XWZ # LSE as initial estimate

      gamma_est <- list()
      gamma_est[[1]] <- gamma_ini
      max_iter <- 100
      iter <- 1
      while(iter <= max_iter){

        eta <- list()
        link_eta <- list()
        wei_mu <- list()
        for(i in 1:n){
          eta[[i]] <- X_tilde[[i]] %*% gamma_est[[iter]]
          link_eta[[i]] <- sapply(eta[[i]], linkfun)
          wei_mu[[i]] <- sapply(link_eta[[i]], linkinvder)
        }

        Z <- list()
        for(i in 1:n){
          Z[[i]] <- list()
          for(j in 1:Y_obser_num[i]){
            Z[[i]][[j]] <- eta[[i]] + (Y[[i]][j] - link_eta[[i]]) * wei_mu[[i]]
          }
        }

        XW_wZ <- matrix(0, ncol = 1, nrow = 2 * p)
        for(i in 1:n){
          for(j in 1:Y_obser_num[[i]]){
            XW_wZ <- XW_wZ + XW[[i]][[j]] %*% diag(1/wei_mu[[i]],
                                                   nrow = length(wei_mu[[i]]),
                                                   ncol = length(wei_mu[[i]])) %*% Z[[i]][[j]]
          }
        }

        XW_wX <- matrix(0, ncol = 2 * p, nrow = 2 * p)
        for(i in 1:n){
          for(j in 1:Y_obser_num[[i]]){
            XW_wX <- XW_wX + XW[[i]][[j]] %*% diag(1/wei_mu[[i]],
                                                   nrow = length(wei_mu[[i]]),
                                                   ncol = length(wei_mu[[i]])) %*% X_tilde[[i]]
          }
        }

        ##add the computation of U
        U_sub <- Comp_U(beta_basis, gamma_est[[iter]][(p+1):(2 * p)],
                        lambda, absTol = absTol)
        comp_id <- c(1:p, U_sub$id + p)
        U_res <- as.matrix(Matrix::bdiag(diag(0, p), U_sub$U[U_sub$id, U_sub$id]))
        XW_wX_inv <- tryCatch(solve(XW_wX[comp_id, comp_id] + n * roupen_para * V_der[comp_id, comp_id] +
                                      2 * n * U_res),
                              error = function(err){return(Inf)})
        if(is.infinite(XW_wX_inv[1])){
          print(paste("smooth", roupen_para, ", sparse", lambda, "break"))
          break
        }
        gamma_est[[iter + 1]] <- matrix(0, nrow = 2 * p, ncol = 1)
        gamma_est[[iter + 1]][comp_id] <- XW_wX_inv %*% XW_wZ[comp_id]

        # print(iter)
        iter <- iter + 1

        if(norm(gamma_est[[iter]] - gamma_est[[iter - 1]], type = "F") <= 0.00001){
          break
        }

      }

      gamma_res[[ii]][[jj]] <- gamma_est[[iter]]
      if(is.infinite(XW_wX_inv[1])){
        EBIC_score[ii, jj] <- Inf
      }else{
        EBIC_score[ii, jj] <- EBIC(X_tilde, Y, family, gamma_est = gamma_est[[iter]], W,
                                   XW_wX, V_der, roupen_para, n_real, Y_obser_num, linkfun)
      }

    }

    # print(ii)
  }

  ii <- which(EBIC_score == min(EBIC_score), arr.ind = T)[1, 1]
  jj <- which(EBIC_score == min(EBIC_score), arr.ind = T)[1, 2]

  results <- list()
  results$gamma_est <- gamma_res[[ii]][[jj]]
  results$roupen_select <- roupen_para_list[ii]
  results$lambda_select <- lambda_list[jj]
  results$EBIC <- EBIC_score
  results$h <- h

  return(results)

}

CV_ind <- function(X_pre, Y_pre, family, X_obser_pre, Y_obser_pre, gamma_est, timeint,
                   L, d, h, kernfun, linkfun){

  n_pre <- length(X_pre)
  obser_num_Y <- unlist(lapply(Y_pre, function(x){length(x)}))

  timecom <- NULL
  timecom_ind <- list()
  timediff_ind <- list()
  for(i in 1:n_pre){
    timecom_ind[[i]] <- expand.grid(X_obser_pre[[i]], Y_obser_pre[[i]])
    timecom <- rbind(timecom, timecom_ind[[i]])
    timediff_ind[[i]] <- timecom_ind[[i]][,1] - timecom_ind[[i]][,2]
  }
  timediff <- timecom[,1] - timecom[,2]

  X_pre_tilde <- list()
  beta_B_pre <- list()
  for(i in 1:n_pre){
    beta_B_pre[[i]] <- splines::bs(X_obser_pre[[i]],
                          knots = seq(timeint[1], timeint[2], length.out = L - d + 1)[-c(1, L - d + 1)],
                          degree = d, intercept = T)
    X_pre_tilde[[i]] <- cbind(beta_B_pre[[i]], as.vector(X_pre[[i]]) * beta_B_pre[[i]])
  }

  W_pre <- list()
  for(i in 1:n_pre){
    W_pre[[i]] <- list()
    for(j in 1:obser_num_Y[i]){
      timediff_ij <- timediff_ind[[i]][which(timecom_ind[[i]][,2] == Y_obser_pre[[i]][j])]
      omega_ij <- sapply(timediff_ij/h, kernfun)/h
      W_pre[[i]][[j]] <- diag(omega_ij, ncol = length(timediff_ij), nrow = length(timediff_ij))
    }
  }

  # sse
  if(family == "Gaussian"){
    SSE <- NULL
    for(i in 1:n_pre){
      for(j in 1:obser_num_Y[i]){
        Y_hat_ik <- linkfun(X_pre_tilde[[i]] %*% gamma_est)
        resi_ij <- (Y_pre[[i]][j] - Y_hat_ik)^2/2
        SSE_ij <- 2 * t(resi_ij) %*% diag(W_pre[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }else if(family == "binomial"){
    SSE <- NULL
    for(i in 1:n_pre){
      for(j in 1:obser_num_Y[i]){
        Y_hat_ik <- linkfun(X_pre_tilde[[i]] %*% gamma_est)
        if(Y_pre[[i]][j] == 1){
          resi_ij <- Y_pre[[i]][j] * log(Y_pre[[i]][j]/Y_hat_ik)
        }else{
          resi_ij <- (1 - Y_pre[[i]][j]) * log((1 - Y_pre[[i]][j])/(1 - Y_hat_ik))
        }
        SSE_ij <- 2 * t(resi_ij) %*% diag(W_pre[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }else if(family == "poisson"){
    SSE <- NULL
    for(i in 1:n_pre){
      for(j in 1:obser_num_Y[i]){
        Y_hat_ik <- linkfun(X_pre_tilde[[i]] %*% gamma_est)
        resi_ij <- Y_hat_ik - Y_pre[[i]][j] * log(Y_hat_ik)
        SSE_ij <- 2 * t(resi_ij) %*% diag(W_pre[[i]][[j]])
        SSE <- c(SSE, SSE_ij)
      }
    }
  }

  MSE <- mean(SSE)
  return(MSE)

}

