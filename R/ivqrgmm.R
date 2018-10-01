#' @title GMM Estimation for Instrumental Variable Quantile Regression Model 
#' (IVQR_GMM)
#' 
#' @param y vector of outcomes.
#' @param w (n by k) matrix of the covariate dataset.
#' @param z (n by p ) matrix of the instrument variable dataset.
#' @param tau quantile index
#' @param intercept FALSE by default, and the function will NOT add intercept
#'                  term automatically, include it in w and z if needed.
#' @param tlimit scalar, 0 by default. If tlimit>0, then tlimit is the
#'               time limit specified for early termination of the MIO solver. 
#'               Otherwise, the MIO solver keeps running until convergence.
#' @param abgap scalar, 0 by default. The absolute gap specified 
#'              for early termination of the MIO solver.
#' @param bnd (k by 2) matrix, NA by default. The first and second columns of 
#'            the matrix respectively store the lower and upper bounds of the 
#'            unknown coefficients.
#' @import gurobi
#' @import matlib
#' 
#' @return a list with 6 elements:
#' \item{theta_hat}{the vector of the coefficient estimates}
#' \item{s_hat}{the estimated asymptotic standard errors}
#' \item{obj_v}{the value of the GMM objective function}
#' \item{gap}{the MIO optimization gap value in case of early termination
#             gap = 0 ==> optimal solution is found within the time limit}
#' \item{rtime}{the time used by the MIO solver in the estimation procedure}
#' \item{ncount}{the number of nodes already explored by the MIO solver}
#' 
#' @description
#' This is a package for Instrumental Variable Quantile Regression (IVQR)
#' using Generalized Methods of Moments (GMM) estimation.
#' It is based on the paper: Chen, Le-Yu and Lee, Sokbae (September 2017),
#' "Exact computation of GMM estimators for instrumental variable quantile 
#' regression models".
#' The paper has been published at Journal of Applied Econometrics.
#' @author Zizhe Xia <xiazizhejordan@@gmail.com>
#' @references \emph{Exact computation of GMM estimators for instrumental 
#' variable quantile regression models}
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/jae.2619}
#' by Chen, Le-Yu and Lee, Sokbae (September 2017).
#' The R codes are based on the MATLAB codes available at:
#' \url{https://github.com/LeyuChen/IVQR-GMM-computation-codes}
#' @export


library(matlib)
library(gurobi)

Two_Stage_LS <- function(y, x, z, robust=1) {
  # This function computes the coefficient estimates and the estimated asymptotic variance
  # for the two-stage least square (2SLS) regression of y on x using z as instrument
  #
  # Args:
  #   y     : the outcome column vector as a n by 1 matrix/array
  #   x     : (n by k) matrix of covariate data where k is the number of covariates
  #   z     : (n by p) matrix of instrument data where p is the number of instruments
  #   robust: set robust = 1 for the estimated heteroskedasticity robust asymptotic variance
  #
  # Returns:
  #   bhat  : the vector of 2SLS regression coefficient estimates
  #   avar  : the estimated 2SLS asymptotic variance
  n <- nrow(y)
  k <- ncol(x)
  P = z %*% solve(t(z) %*% z,t(z))
  xhat <- P %*% x
  bhat <- solve(t(xhat) %*% xhat, t(xhat) %*% y)
  uhat <- y - x %*% bhat
  inv_xhat_xhat <- solve(t(xhat) %*% xhat, diag(1, k))
  
  if (robust == 1){
    var <- inv_xhat_xhat %*% (t(xhat) %*% diag(uhat*uhat) %*% xhat) %*% inv_xhat_xhat
  } else {
    avar= ((t(uhat) %*% uhat) / n) %*% inv_xhat_xhat
  }
  list(bhat=bhat, avar=avar)
}

miobnd_fn <- function(y, x, bnd) {
  # Given (y,x), this function solves the following maximization problem
  # for each i, max |y(i)-x(i,:)*b| over b confined to the space described by bnd
  #
  # Args:
  # y     : (n by 1) matrix of outcomes
  # x     : (n by k) matrix of covariate data
  # bnd   : (k by 2) matrix where the first and second columns respectively store
  #         the lower and upper bounds of the unknown coefficients
  #
  # Returns:
  # value : the value of the maximized objective function
  n <- nrow(y)
  k <- ncol(x)
  
  model <- list()
  params <- list()
  
  model$modelsense <- 'max'
  model$sense <- '>'
  model$lb <- bnd[, 1]
  model$ub <- bnd[, 2]
  
  tol <- 1e-6
  
  params$outputflag <- 0
  params$OptimalityTol <- tol
  params$FeasibilityTol <- tol
  params$IntFeasTol <- tol
  
  value <- matrix(0, n, 1)
  
  for(i in 1:n){
    v <- matrix(0, 2, 1) 
    
    model$obj <- -x[i, ]
    model$objcon <- y[i]
    
    tryCatch({
      model$A <- -x[i, ]
      model$rhs <- -y[i]
      result <- gurobi(model, params)
      # %disp(result.x);
      v[1] <- result$objval
      # %rtime=result.runtime;
      # %disp(v(1));
      # %fprintf('Optimization returned status: %s\n', result.status);
    }, warning=function(gurobiError) {
      cat('Error reported\n')
      # %fprintf('Error reported\n');disp(v(1));
    })
    
    model$obj <- x[i, ]
    model$objcon <- -y[i]
    
    tryCatch({
      model$A <- x[i, ]
      model$rhs <- y[i]
      result <- gurobi(model, params)
      # %disp(result.x)
      v[2] <- result$objval
      # %rtime=result.runtime
      # %disp(v(2))
      # %fprintf('Optimization returned status: %s\n', result.status);
    }, warning=function(gurobiError) {
      cat('Error reported\n')
      # %fprintf('Error reported\n');disp(v(2));
    })
    
    value[i] <- max(v)
  }
  value
}

IVQR_MIO <- function(y, x, Q, tau, tlimit, abgap, bnd, method=1) {
  # Args:
  #   y     : vector of outcomes
  #   x     : (n by k) matrix of the covariate data
  #   Q     : (n by n) matrix equal to (G*Q_hat*G') stated in the MIQP formulation 
  #   tau   : quantile index
  #   tlimit     : the time limit specified for early termination of the MIO solver
  #   abgap : the absolute gap specified for early termination of the MIO solver
  #   bnd   : (k by 2) matrix where the first and second columns  
  #   respectively store the lower and upper bounds 
  #   of the unknown coefficients
  #   method: set method=1 for solving the MIQP formulation (3.3)
  #   set method=2 for solving the MIQP formulation (C.1)
  #   set method=3 for solving the MILP formulation (C.10)
  # 
  # Returns:
  #   bhat  : the vector of the coefficient estimates
  #   obj_v : the value of the GMM objective function
  #   gap   : the MIO optimization gap value in case of early termination
  #           gap = 0 ==> optimal solution is found
  #   rtime : the time used by the MIO solver in the estimation procedure
  #   ncount: the number of nodes already explored by the MIO solver 
  n <- nrow(y)
  k <- ncol(x)
  bhat <- matrix (0, k, 1)
  
  gap <- 0
  rtime <- 0
  ncount <- 0
  
  tau_vec <- matrix(tau, n, 1)
  model$objcon <- t(tau_vec) %*% Q %*% tau_vec
  
  model$modelsense <- 'min'
  tol <- 1e-6
  
  if (method==1) { # MIQP formulation (3.3)
    print('solving MIQP formulation (3.3)')
    
    model$sense <- '<';
    model$lb <- c(matrix(0, n, 1), bnd[ ,1])
    model$ub <- c(matrix(1, n, 1), bnd[ ,2])
    # 'B' : int code 66
    # 'C' : int code 67
    model$vtype <- chr(c(matrix(66, 1, n), matrix(67, 1, k)))
    
    miobnd <- miobnd_fn(y, x, bnd)  # miobnd_fn computes the values M(i) defined in (3.6)
    miobnd_bar <- miobnd + tol
    
    model$obj <- c(-2 %*% Q %*% tau_vec, matrix(0, k, 1))
    model$Q <- rbind(cbind(Q, matrix(0, n, k)), matrix(0, k, n+k))
    model$A <- rbind(cbind(diag(c(miobnd)), -x), cbind(-diag(c(miobnd_bar)), x))
    model$rhs <- c(miobnd * (1 - tol) - y, y - tol * miobnd_bar)
  } else if (method==2) { # MIQP formulation (C.1)
    print('solving MIQP formulation (C.1)')
    
    model$sense <- c(rep('=', times=2*n), rep('<', times=n))
    model$lb <- c(matrix(0, n, 1), bnd[ ,1], matrix(0, 3*n, 1))
    model$ub <- c(matrix(1, n, 1), bnd[ ,2], matrix(1, n, 1), matrix(1/eps, 2*n, 1))

    # 'B' : int code 66
    # 'C' : int code 67
    model$vtype <- chr(c(matrix(66, 1, n), matrix(67, 1, k), matrix(66, 1, n), matrix(67, 1, 2*n)))

    model$obj <- c(-2 * Q %*% tau_vec, matrix(0, 3*n+k, 1))
    Q.matrix <- Matrix(0, nrow=4*n+k, ncol=4*n+k, sparse=TRUE)
    Q.matrix[c(1:n), c(1:n)] <- Q
    model$Q <- Q.matrix
    # Matrix(cbind(rbind(Q, matrix(0, n, 3*n+k)), matrix(0, 3*n+k, 4*n+k))) # can be sparse
    
    model$A <- Matrix(rbind(cbind(diag(1, n), matrix(0, n, k), diag(1, n), matrix(0, n, 2*n)), 
                            cbind(matrix(0, n, n), x, matrix(0, n, n), diag(1, n), diag(-1, n)), 
                            cbind(matrix(0, n, 2*n+k), diag(-1, n), diag(-1, n))), sparse=TRUE)
    model$rhs <- c(matrix(1, n, 1), y, matrix(-1e-5, n, 1))

    params$PreSOS1BigM <- 0

    # Specification of the SOS-1 constraints
    sos.all <- list()
    for (j in 1:n) {
      sos.1 <- list()
      sos.1$type <- 1
      sos.1$index <- c(j, 2*n+k+j) # (r,e): SOS-1
      sos.all[[j]] <- sos.1
      sos.2 <- list()
      sos.2$type <- 1
      sos.2$index <- c(n+k+j, 3*n+k+j) # (s,1-e): SOS-1
      sos.all[[n+j]] <- sos.2
    }
    model$sos <- sos.all
  } else if (method == 3) { # MILP formulation (C.10)
    print('solving MILP formulation (C.10)') 

    aux_num <- n*(n-1)/2

    model$sense <- '<'
    model$lb <- c(matrix(0, n, 1), bnd[ ,1], matrix(0, aux_num, 1))
    model$ub <- c(matrix(1, n, 1), bnd[ ,2], matrix(1, aux_num, 1))

    # 'B' : int code 66
    # 'C' : int code 67
    model$vtype <- chr(c(matrix(66, 1, n), matrix(67, 1, k), matrix(66, 1, aux_num)))

    miobnd <- miobnd_fn(y, x, bnd)  # miobnd_fn computes the values M(i) defined in (3.6)
    miobnd_bar <- miobnd + tol

    Q_vecl <- matrix()
    aux_constr1 <- matrix(0, aux_num, n+k+aux_num)
    aux_constr2 <- matrix(0, aux_num, n+k+aux_num)
    aux_constr3 <- matrix(0, aux_num, n+k+aux_num)

    s <- 0
    for (i in 1:n-1) {
      if (i == 1) {
        Q_vecl <- Q[i+1:n,i]
      } else {
        Q_vecl <- rbind(Q_vecl, Q[i+1:n,i])
      }
  
      # -e_i+x_ij <= 0
      aux_constr1[s+1:s+n-i, i] <- matrix(-1, n-i, 1)
      aux_constr1[s+1:s+n-i, n+k+s+1:n+k+s+n-i] <- diag(1, n-i)
  
      # -e_j+x_ij <= 0
      aux_constr2[s+1:s+n-i, i+1:n] <- diag(-1, n-i)
      aux_constr2[s+1:s+n-i, n+k+s+1:n+k+s+n-i] <- diag(-1, n-i)
  
      # e_i+e_j-x_ij <= 1
      aux_constr3[s+1:s+n-i, i]<- matrix(1, n-i, 1)
      aux_constr3[s+1:s+n-i, i+1:n]<- diag(1, n-i)
      aux_constr3[s+1:s+n-i, n+k+s+1:n+k+s+n-i] <- diag(-1, n-i)
  
      s <- n-i+s
    }

    model$obj <- c(diag(Q)-2*Q %*% tau_vec, matrix(0, k, 1), 2*Q_vecl)
    model$A <- Matrix(rbind(cbind(diag(c(miobnd)), -x, matrix(0, n, aux_num)), 
                            cbind(-diag(c(miobnd_bar)), x, matrix(0, n, aux_num)), 
                            aux_constr1, aux_constr2, aux_constr3), sparse=TRUE)
    model$rhs <- c(miobnd*(1 - tol) - y, y - tol*miobnd_bar, matrix(0, 2*aux_num, 1), matrix(1, aux_num, 1))
  } else {
    print('error in input arguments')
    return
  }
  
  params$outputflag <- 0
  params$OptimalityTol <- tol
  params$FeasibilityTol <- tol
  params$IntFeasTol <- tol
  
  if (tlimit > 0) {
    params$TimeLimit <- tlimit
  }
  
  if (abgap > 0) {
    params$MIPGapAbs <- abgap
  }
  
  tryCatch({
    result <- gurobi(model, params)
    bhat <- result$x[n+1:n+k]
    obj_v <- result$objval
    gap <- (obj_v - result$objbound)
    rtime <- result$runtime
    ncount <- result$nodecount
    # cat(sprintf(('Optimization returned status: %s\n', result$status)))
    print('Optimization returned status: ')
    print(result$status)
  }, warning=function(gurobiError) {
    cat('Error reported\n')
  })
  list(bhat=bhat, obj_v=obj_v, gap=gap, rtime=rtime, ncount=ncount)
}

IVQR_GMM <- function(y, w, z, tau, intercept=FALSE, tlimit=0, abgap=0, bnd=NA, method=1) {
  # The IVQR_GMM function computes the exact GMM estimator of the IVQR model
  # via the MIO approach as described in Chen and Lee (2017).
  # 
  # Args:
  #   y        : vector of outcomes
  #   w        : (n by k) matrix of the covariate dataset
  #   z        : (n by p ) matrix of the instrument variable dataset
  #   tau      : quantile index
  #   intercept: False ==> The function will NOT add intercept term automatically, 
  #                        include it in w and z if needed
  #              True  ==> The function will ADD intercept term to w and z automatically
  #   tlimit   : scalar. If tlimit>0, then tlimit is the time limit specified for early termination
  #              of the MIO solver. Otherwise, the MIO solver keeps running until convergence.
  #   abgap    : the absolute gap specified for early termination of the MIO solver
  #   bnd      : (k by 2) matrix where the first and second columns  
  #              respectively store the lower and upper bounds of the unknown coefficients         
  #   method   : 1 set method = 1 for computing the GMM estimator based on the formulation (3.3) 
  #                               of Chen and Lee (2017)
  #              2 set method = 2 for using the formulation (C.1)
  #              3 set method = 3 for using the formulation (C.10)
  # 
  # The arguments tlimit, abgap and bnd are optional. When they are not specified,
  # the following default values are used:
  #   tlimit   : set T = 0  ==> solve the MIO problem until convergence
  #   intercept: set intercept=False, the function will NOT add intercept term by itself
  #   abgap    : set abgap = 0  ==> solve the MIO problem until convergence
  #   bnd      : Calculate the parameter bounds based on the two-stage least 
  #              square regression results as used in Chen and Lee (2017)
  # 
  # Returns:
  #   theta_hat: the vector of the coefficient estimates
  #   s_hat    : the estimated asymptotic standard errors
  #   obj_v    : the value of the GMM objective function
  #   gap      : the MIO optimization gap value in case of early termination
  #              gap = 0 ==> optimal solution is found within the time limit
  #   rtime    : the time used by the MIO solver in the estimation procedure
  #   ncount   : the number of nodes already explored by the MIO solver 
  n <- nrow(y)
  
  if (intercept) {
    w <- cbind(1, w)
    z <- cbind(1, z)
  }
  
  # Use of the Hall-Sheath bandwidth choice (Koenker 1994)
  q_tau <- qnorm(tau)
  H_S_ratio <- (dnorm(q_tau) ^ 2) / (1 + 2 * q_tau * q_tau)
  h_Hall_Sheath <- (1.96 ^ (2 / 3)) * ((1.5 / n * H_S_ratio) ^ (1 / 3))
  
  Q <- z * inv(t(z) %*% z / n) * (t(z) / (tau * (1 - tau))) 
  # Q is the matrix G %*% Q_hat %*% G' stated in the MIQP formulation of the GMM
  # estimation problem
  
  k <- ncol(w)
  theta_hat <- matrix(0, k, 1)
  
  if (is.na(bnd)) {
    TSLS_result <- Two_stage_LS(y,w,z,1)
    b <- TSLS_result$bhat
    var <- TSLS_result$avar
    bnd <- cbind(b - 10 * sqrt(diag(var)), b + 10 * sqrt(diag(var)))
  }
  
  mio_result <- IVQR_MIO(y,w,Q,tau,T,abgap,bnd,method)
  theta_hat <- mio_result$bhat
  obj_v <- mio_result$obj_v
  gap <- mio_result$gap
  rtime <- mio_result$rtime
  ncount <- mio_result$ncount
  
  # compute the estimated asymptotic standard errors based on the method of
  # Powell (1986) using Gaussian kernal and the Hall-Sheath bandwidth choice
  e_hat <- y - w %*% theta_hat
  kern <- dnorm(e_hat / h_Hall_Sheath) / h_Hall_Sheath
  k_x <- kronecker(matrix(1, 1, k), kern) * w
  s_hat <- sqrt(diag(inv(t(k_x) %*% z %*% inv(t(z) %*% z) %*% t(z) %*% k_x) * tau * (1-tau)))
  list(theta_hat=theta_hat, s_hat=s_hat, obj_v=obj_v, gap=gap, rtime=rtime, ncount=ncount)
}
