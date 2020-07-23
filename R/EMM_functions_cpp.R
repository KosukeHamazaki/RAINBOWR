#' Perform spectral decomposition (inplemented by Rcpp)
#'
#' @description Perform spectral decomposition for \eqn{G = ZKZ'} or \eqn{SGS} where \eqn{S = I - X(X'X)^{-1}X}.
#'
#' @param ZETA A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param ZWs A list of additional linear kernels other than genomic relationship matrix (GRM).
#' We utilize this argument in RGWAS.multisnp function, so you can ignore this.
#' @param X \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param weights If the length of ZETA >= 2, you should assign the ratio of variance components to this argument.
#' @param return.G If thie argument is TRUE, spectral decomposition results of G will be returned.
#' (\eqn{G = ZKZ'})
#' @param return.SGS If this argument is TRUE, spectral decomposition results of SGS will be returned.
#' (\eqn{S = I - X(X'X)^{-1}X}, \eqn{G = ZKZ'})
#' @param spectral.method The method of spectral decomposition.
#' In this function, "eigen" : eigen decomposition and "cholesky" : cholesky and singular value decomposition are offered.
#' If this argument is NULL, either method will be chosen accorsing to the dimension of Z and X.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param df.H The degree of freedom of K matrix. If this argument is NULL, min(n, sum(nrow(K1), nrow(K2), ...)) will be assigned.
#'
#' @return
#' \describe{
#' \item{$spectral.G}{The spectral decomposition results of G.}
#' \item{$U}{Eigen vectors of G.}
#' \item{$delta}{Eigen values of G.}
#' \item{$spectral.SGS}{Estimator for \eqn{\sigma^2_e}}
#' \item{$Q}{Eigen vectors of SGS.}
#' \item{$theta}{Eigen values of SGS.}
#' }
#'
#'
#'
spectralG.cpp <- function(ZETA, ZWs = NULL, X = NULL, weights = 1, return.G = TRUE,
                          return.SGS = FALSE, spectral.method = NULL,
                          tol = NULL, df.H = NULL){
  if((!is.null(tol)) & (length(tol) == 1)){
    tol <- rep(tol, 2)
  }

  lines.name.pheno <- rownames(ZETA[[1]]$Z)
  n <- nrow(ZETA[[1]]$Z)
  ms.ZETA <- unlist(lapply(ZETA, function(x) ncol(x$Z)))
  m.ZETA <- sum(ms.ZETA)
  lz <- length(ZETA)
  if(!is.null(ZWs)){
    lw <- length(ZWs)
    ms.ZWs <- unlist(lapply(ZWs, function(x) ncol(x$Z)))
    m.ZWs <- sum(ms.ZWs)
  }else{
    lw <- 0
    ms.ZWs <- NULL
    m.ZWs <- 0
  }

  m <- m.ZETA + m.ZWs

  if((lz + lw) != length(weights)){
    stop("Weights should have the same length as ZETA!")
  }
  stopifnot(all(weights >= 0))

  if(is.null(X)){
    X <- as.matrix(rep(1, n))
    rownames(X) <- lines.name.pheno
    colnames(X) <- "Intercept"
  }
  p <- ncol(X)

  if(is.null(spectral.method)){
    if(n <= m + p){
      spectral.method <- "eigen"
    }else{
      spectral.method <- "cholesky"
    }
  }

  if(is.null(df.H)){
    df.H <- n
  }

  if(spectral.method == "cholesky"){
    ZBt <- NULL
    for(i in 1:lz){
      Z.now <- ZETA[[i]]$Z
      K.now <- ZETA[[i]]$K

      diag(K.now) <- diag(K.now) + 1e-06
      B.now <- try(chol(K.now), silent = TRUE)
      if ("try-error" %in% class(B.now)) {
        stop("K not positive semi-definite.")
      }

      ZBt.now <- tcrossprod(Z.now, B.now)
      ZBt <- cbind(ZBt, sqrt(weights[i]) * ZBt.now)
    }

    if(!is.null(ZWs)){
      for(j in 1:lw){
        Z.now <- ZWs[[j]]$Z
        Bt.now <- ZWs[[j]]$W %*% ZWs[[j]]$Gamma

        ZBt.now <- Z.now %*% Bt.now
        ZBt <- cbind(ZBt, sqrt(weights[lz + j]) * ZBt.now)
      }
    }
    spectral_G.res <- try(spectralG_cholesky(zbt = ZBt, x = X, return_G = return.G,
                                             return_SGS = return.SGS),
                          silent = TRUE)

    if(!("try-error" %in% class(spectral_G.res))){

      if(return.G){
        U0 <- spectral_G.res$U
        delta0 <- c(spectral_G.res$D)

        if(!is.null(tol)){
          G.tol.over <- which(delta0 > tol[2])
        }else{
          G.tol.over <- 1:df.H
        }
        delta <- delta0[G.tol.over]
        U <- U0[, G.tol.over]
      }

      if(return.SGS){
        Q0 <- spectral_G.res$Q
        theta0 <- c(spectral_G.res$theta)

        if(!is.null(tol)){
          SGS.tol.over <- which(theta0 > tol[2])
        }else{
          SGS.tol.over <- 1:(df.H - p)
        }
        theta <- theta0[SGS.tol.over]
        Q <- Q0[, SGS.tol.over]
      }
    }else{
      spectral.method <- "eigen"
    }
  }


  if(spectral.method == "eigen"){
    ZKZt <- matrix(0, nrow = n, ncol = n)
    for(i in 1:lz){
      Z.now <- ZETA[[i]]$Z
      K.now <- ZETA[[i]]$K

      ZKZt.now <- tcrossprod(Z.now %*% K.now, Z.now)
      ZKZt <- ZKZt + weights[i] * ZKZt.now
    }

    if(!is.null(ZWs)){
      for(j in 1:lw){
        Z.now <- ZWs[[j]]$Z
        K.now <- tcrossprod(ZWs[[j]]$W %*% ZWs[[j]]$Gamma, ZWs[[j]]$W)

        ZKZt.now <- tcrossprod(Z.now %*% K.now, Z.now)
        ZKZt <- ZKZt + weights[lz + j] * ZKZt.now
      }
    }
    spectral_G.res <- spectralG_eigen(zkzt = ZKZt, x = X, return_G = return.G,
                                      return_SGS = return.SGS)

    if(return.G){
      U0 <- spectral_G.res$U
      delta0 <- c(spectral_G.res$D)

      if(!is.null(tol)){
        G.tol.over <- which(delta0 > tol[2])
      }else{
        G.tol.over <- 1:df.H
      }
      delta <- delta0[G.tol.over]
      U <- U0[, G.tol.over]
    }

    if(return.SGS){
      Q0 <- spectral_G.res$Q
      theta0 <- c(spectral_G.res$theta)

      if(!is.null(tol)){
        SGS.tol.over <- which(theta0 > tol[2])
      }else{
        SGS.tol.over <- 1:(df.H - p)
      }
      theta <- theta0[SGS.tol.over]
      Q <- Q0[, SGS.tol.over]
    }
  }


  if(return.G){
    spectral.G <- list(U = U, delta = delta)
  }else{
    spectral.G <- NULL
  }

  if(return.SGS){
    spectral.SGS <- list(Q = Q, theta = theta)
  }else{
    spectral.SGS <- NULL
  }


  return(list(spectral.G = spectral.G,
              spectral.SGS = spectral.SGS))
}










#' Equation of mixed model for one kernel, GEMMA-based method (inplemented by Rcpp)
#'
#' @description This function solves the single-kernel linear mixed effects model by GEMMA
#' (genome wide efficient mixed model association; Zhou et al., 2012) approach.
#'
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param ZETA A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param lam.len The number of initial values you set. If this number is large, the estimation will be more accurate,
#' but computational cost will be large. We recommend setting this value 3 <= lam.len <= 6.
#' @param init.range The range of the initial parameters. For example, if lam.len = 5 and init.range = c(1e-06, 1e02),
#' corresponding initial heritabilities will be calculated as seq(1e-06, 1 - 1e-02, length = 5),
#' and then initial lambdas will be set.
#' @param init.one The initial parameter if lam.len = 1.
#' @param conv.param The convergence parameter. If the diffrence of log-likelihood by updating the parameter "lambda"
#' is smaller than this conv.param, the iteration steps will be stopped.
#' @param count.max Sometimes algorithms won't converge for some initial parameters.
#' So if the iteration steps reache to this argument, you can stop the calculation even if algorithm doesn't converge.
#' @param bounds Lower and Upper bounds of the parameter 1 / lambda. If the updated parameter goes out of this range,
#' the parameter is reset to the value in this range.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param REML You can choose which method you will use, "REML" or "ML".
#' If REML = TRUE, you will perform "REML", and if REML = FALSE, you will perform "ML".
#' @param silent If this argument is TRUE, warning messages will be shown when estimation is not accurate.
#' @param plot.l If you want to plot log-likelihood, please set plot.l = TRUE.
#' We don't recommend plot.l = TRUE when lam.len >= 2.
#' @param SE If TRUE, standard errors are calculated.
#' @param return.Hinv If TRUE, the function returns the inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}. This is useful for GWAS.
#'
#' @return
#' \describe{
#' \item{$Vu}{Estimator for \eqn{\sigma^2_u}}
#' \item{$Ve}{Estimator for \eqn{\sigma^2_e}}
#' \item{$beta}{BLUE(\eqn{\beta})}
#' \item{$u}{BLUP(\eqn{u})}
#' \item{$LL}{Maximized log-likelihood (full or restricted, depending on method)}
#' \item{$beta.SE}{Standard error for \eqn{\beta} (If SE = TRUE)}
#' \item{$u.SE}{Standard error for \eqn{u^*-u} (If SE = TRUE)}
#' \item{$Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I} (If return.Hinv = TRUE)}
#' \item{$Hinv2}{The inverse of \eqn{H2 = ZKZ'/\lambda + I} (If return.Hinv = TRUE)}
#' \item{$lambda}{Estimators for \eqn{\lambda = \sigma^2_e / \sigma^2_u}}
#' \item{$lambdas}{Lambdas for each initial values}
#' \item{$reest}{If parameter estimation may not be accurate, reest = 1, else reest = 0}
#' \item{$counts}{The number of iterations until convergence for each initial values}
#' }
#'
#' @references Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis
#'  for association studies. Nat Genet. 44(7): 821-824.
#'
#'
#'
#'
EMM1.cpp <- function(y, X = NULL, ZETA, eigen.G = NULL, lam.len = 4, init.range = c(1e-04, 1e02),
                     init.one = 0.5, conv.param = 1e-06, count.max = 15, bounds = c(1e-06, 1e06),
                     tol = NULL, REML = TRUE, silent = TRUE, plot.l = FALSE,
                     SE = FALSE, return.Hinv = TRUE){

  #### The start of EMM1 (GEMMA-based) ####
  if(length(ZETA) != 1){
    stop("If you want to solve multikernel mixed-model equation, you should use EM3.cpp function.")
  }

  ### Some definitions ###
  y0 <- as.matrix(y)
  n0 <- length(y0)

  if (is.null(X)) {
    p <- 1
    X <- matrix(rep(1, n0), n0, 1)
  }
  p <- ncol(X)

  Z <- ZETA[[1]]$Z
  if (is.null(Z)) {
    Z <- diag(n0)
  }
  m <- ncol(Z)
  if (is.null(m)) {
    m <- 1
    Z <- matrix(Z, length(Z), 1)
  }


  stopifnot(nrow(Z) == n0)
  stopifnot(nrow(X) == n0)

  K <- ZETA[[1]]$K
  if (!is.null(K)) {
    stopifnot(nrow(K) == m)
    stopifnot(ncol(K) == m)
  }
  not.NA <- which(!is.na(y0))
  Z <- Z[not.NA, , drop = FALSE]
  X <- X[not.NA, , drop = FALSE]
  n <- length(not.NA)
  y <- matrix(y0[not.NA], n, 1)

  spI <- diag(n)
  logXtX <- sum(log(eigen(crossprod(X), symmetric = TRUE)$values))

  ### Decide initial parameter and calculate H ###
  if(lam.len >= 2){
    h2.0 <- seq(init.range[1], 1 - 1 / init.range[2], length = lam.len)
    #lambda.0 <- h2.0 / (1 - h2.0)
    lambda.0 <- 10 ^ seq(log10(init.range[1]), log10(init.range[2]), length = lam.len)
  }else{
    lambda.0 <- init.one
  }


  ### Eigen decomposition of ZKZt ###
  if(is.null(eigen.G)){
    eigen.G <- spectralG.cpp(ZETA = lapply(ZETA, function(x) list(Z = x$Z[not.NA, , drop = FALSE], K = x$K)),
                             X = X, return.G = TRUE, return.SGS = FALSE, tol = tol, df.H = NULL)[[1]]
  }
  U <- eigen.G$U
  delta <- eigen.G$delta


  if(lam.len >= 2){
    it <- split(lambda.0, factor(1:lam.len))
    mclap.res <- unlist(parallel::mclapply(it, ml_est_out, y = as.matrix(y), x = as.matrix(X),
                                           u = U, delta = as.matrix(delta), z = Z, k = K,
                                           logXtX = logXtX, conv_param = conv.param, count_max = count.max,
                                           bounds1 = bounds[1], bounds2 = bounds[2], REML = REML, mc.cores = lam.len))
    counts <- mclap.res[seq(3, 3 * lam.len, by = 3)]
    l.maxs <- mclap.res[seq(2, 3 * lam.len - 1, by = 3)]
    lambda.opts <- mclap.res[seq(1, 3 * lam.len - 2, by = 3)]
    lambda.opt <- lambda.opts[which.max(l.maxs)]
    names(lambda.opt) <- "lambda.opt"
  }else{
    ml.res <- ml_est_out(lambda.0[1], y = as.matrix(y), x = as.matrix(X),
                         u = U, delta = as.matrix(delta), z = Z, k = K,
                         logXtX = logXtX, conv_param = conv.param, count_max = count.max,
                         bounds1 = bounds[1], bounds2 = bounds[2], REML = REML)
    lambda.opts <- lambda.opt <- ml.res$lambda
    l.maxs <- ml.res$max.l
    counts <- ml.res$count
    names(lambda.opt) <- "lambda.opt"
  }


  if(all(counts[-1] == count.max)){
    if(!silent){
      warning("Parameter estimation may not be accurate. Please reestimate with EMM2.cpp function.")
    }
    reest <- 1
  }else{
    reest <- 0
  }

  res.list <- ml_one_step(lambda.opt, y = as.matrix(y), x = as.matrix(X),
                          u = U, delta = as.matrix(delta), z = Z, k = K,
                          logXtX = logXtX, conv_param = conv.param, count_max = count.max,
                          bounds1 = bounds[1], bounds2 = bounds[2], REML = REML, SE = SE, return_Hinv = return.Hinv)

  return(c(res.list, list(lambda = 1 / lambda.opt), list(lambdas = 1 / lambda.opts),
           list(reest = reest), list(counts = counts)))
}







#' Equation of mixed model for one kernel, EMMA-based method (inplemented by Rcpp)
#'
#' @description This function solves single-kernel linear mixed model by EMMA
#'  (efficient mixed model association; Kang et al., 2008) approach.
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param ZETA A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param bounds Lower and Upper bounds of the parameter lambda. If the updated parameter goes out of this range,
#' the parameter is reset to the value in this range.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param traceInside Perform trace for the optimzation if traceInside >= 1, and this argument shows the frequency of reports.
#' @param REML You can choose which method you will use, "REML" or "ML".
#' If REML = TRUE, you will perform "REML", and if REML = FALSE, you will perform "ML".
#' @param SE If TRUE, standard errors are calculated.
#' @param return.Hinv If TRUE, the function returns the inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}. This is useful for GWAS.
#'
#' @return
#' \describe{
#' \item{$Vu}{Estimator for \eqn{\sigma^2_u}}
#' \item{$Ve}{Estimator for \eqn{\sigma^2_e}}
#' \item{$beta}{BLUE(\eqn{\beta})}
#' \item{$u}{BLUP(\eqn{u})}
#' \item{$LL}{Maximized log-likelihood (full or restricted, depending on method)}
#' \item{$beta.SE}{Standard error for \eqn{\beta} (If SE = TRUE)}
#' \item{$u.SE}{Standard error for \eqn{u^*-u} (If SE = TRUE)}
#' \item{$Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I} (If return.Hinv = TRUE)}
#' }
#'
#' @references Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#'
#'
EMM2.cpp <- function(y, X = NULL, ZETA, eigen.G = NULL, eigen.SGS = NULL, tol = NULL, optimizer = "nlminb",
                     traceInside = 0, REML = TRUE, bounds = c(1e-09, 1e+09), SE = FALSE, return.Hinv = FALSE){

  #### The start of EMM2 (EMMA-based) ####
  if(length(ZETA) != 1){
    stop("If you want to solve multikernel mixed-model equation, you should use EM3.cpp function.")
  }

  pi <- 3.14159
  y0 <- as.matrix(y)
  n0 <- length(y0)

  if (is.null(X)) {
    p <- 1
    X <- matrix(rep(1, n0), n0, 1)
  }
  p <- ncol(X)

  Z <- ZETA[[1]]$Z
  if (is.null(Z)) {
    Z <- diag(n0)
  }
  m <- ncol(Z)
  if (is.null(m)) {
    m <- 1
    Z <- matrix(Z, length(Z), 1)
  }


  stopifnot(nrow(Z) == n0)
  stopifnot(nrow(X) == n0)

  K <- ZETA[[1]]$K
  if (!is.null(K)) {
    stopifnot(nrow(K) == m)
    stopifnot(ncol(K) == m)
  }
  not.NA <- which(!is.na(y0))
  Z <- Z[not.NA, , drop = FALSE]
  X <- X[not.NA, , drop = FALSE]
  n <- length(not.NA)
  y <- matrix(y0[not.NA], n, 1)

  spI <- diag(n)

  ### Eigen decomposition of ZKZt ###
  if(is.null(eigen.G)){
    eigen.G <- spectralG.cpp(ZETA = lapply(ZETA, function(x) list(Z = x$Z[not.NA, , drop = FALSE], K = x$K)),
                             X = X, return.G = TRUE, return.SGS = FALSE, tol = tol, df.H = NULL)[[1]]
  }
  U <- eigen.G$U
  delta <- eigen.G$delta

  if(is.null(eigen.SGS)){
    eigen.SGS <- spectralG.cpp(ZETA = lapply(ZETA, function(x) list(Z = x$Z[not.NA, , drop = FALSE], K = x$K)),
                               X = X, return.G = FALSE, return.SGS = TRUE, tol = tol, df.H = NULL)[[2]]
  }
  Q <- eigen.SGS$Q
  theta <- eigen.SGS$theta

  omega <- crossprod(Q, y)
  omega.sq <- omega ^ 2

  lambda.0 <- 1
  n <- nrow(X)
  p <- ncol(X)
  phi <- delta
  traceNo <- ifelse(traceInside > 0, 3, 0)
  traceREPORT <- ifelse(traceInside > 0, traceInside, 1)

  if (!REML) {
    f.ML <- function(lambda, n, theta, omega.sq, phi) {
      n * log(sum(omega.sq / (theta + lambda))) + sum(log(phi +
                                                            lambda))
    }
    # gr.ML <- function(lambda, n, theta, omega.sq, phi) {
    #   n * sum(omega.sq / (theta + lambda) ^ 2) / sum(omega.sq / (theta + lambda)) -
    #     sum(1 / (phi + lambda))
    # }

    if (optimizer == "optim") {
      soln <- optimize(f.ML, interval = bounds, n, theta,
                       omega.sq, phi)
      lambda.opt <- soln$minimum
      maxval <- soln$objective
    } else if (optimizer == "optimx") {
      soln <- optimx::optimx(par = lambda.0, fn = f.ML, gr = NULL, hess = NULL,
                             lower = bounds[1], upper = bounds[2],
                             method = "L-BFGS-B", n = n, theta = theta,
                             omega.sq = omega.sq, phi = phi,
                             control = list(trace = traceNo, starttests = FALSE, maximize = F,
                                            REPORT = traceREPORT, kkt = FALSE))
      lambda.opt <- soln[, 1]
      maxval <- soln[, 2]
    } else if (optimizer == "nlminb") {
      soln <- nlminb(start = lambda.0, objective = f.ML, gradient = NULL, hessian = NULL,
                     lower = bounds[1], upper = bounds[2], n = n, theta = theta,
                     omega.sq = omega.sq, phi = phi, control = list(trace = traceInside))
      lambda.opt <- soln$par
      maxval <- soln$objective
    } else {
      warning("We offer 'optim', 'optimx', and 'nlminb' as optimzers. Here we use 'optim' instead.")
      soln <- optimize(f.ML, interval = bounds, n, theta,
                       omega.sq, phi)
      lambda.opt <- soln$minimum
      maxval <- soln$objective
    }

    df <- n
  } else {
    f.REML <- function(lambda, n.p, theta, omega.sq) {
      n.p * log(sum(omega.sq / (theta + lambda))) + sum(log(theta + lambda))
    }
    # gr.REML <- function(lambda, n.p, theta, omega.sq) {
    #   n.p * sum(omega.sq / (theta + lambda) ^ 2) / sum(omega.sq / (theta + lambda)) -
    #     sum(1 / (theta + lambda))
    # }

    if (optimizer == "optim") {
      soln <- optimize(f.REML, interval = bounds, n - p, theta,
                       omega.sq)
      lambda.opt <- soln$minimum
      maxval <- soln$objective
    } else if (optimizer == "optimx") {
      soln <- optimx::optimx(par = lambda.0, fn = f.REML, gr = NULL, hess = NULL,
                             lower = bounds[1], upper = bounds[2],
                             method = "L-BFGS-B", n.p = n - p, theta = theta,
                             omega.sq = omega.sq, control = list(trace = traceNo, starttests = FALSE,
                                                                 maximize = F, REPORT = traceREPORT, kkt = FALSE))
      lambda.opt <- soln[, 1]
      maxval <- soln[, 2]
    } else if (optimizer == "nlminb") {
      soln <- nlminb(start = lambda.0, objective = f.REML, gradient = NULL, hessian = NULL,
                     lower = bounds[1], upper = bounds[2], n.p = n - p, theta = theta,
                     omega.sq = omega.sq, control = list(trace = traceInside))
      lambda.opt <- soln$par
      maxval <- soln$objective
    } else {
      warning("We offer 'optim', 'optimx', and 'nlminb' as optimzers. Here we use 'optim' instead.")
      soln <- optimize(f.REML, interval = bounds, n - p, theta,
                       omega.sq)
      lambda.opt <- soln$minimum
      maxval <- soln$objective
    }

    df <- n - p
  }


  EMM2.final.res <- EMM2_last_step(lambda.opt, as.matrix(y), X, U, as.matrix(delta),
                                   Z, K, as.matrix(theta), as.matrix(omega.sq),
                                   as.matrix(phi), maxval, n, p, df, SE, return.Hinv)
  return(EMM2.final.res)
}






#' Equation of mixed model for one kernel, a wrapper of two methods
#'
#' @description This function estimates maximum-likelihood (ML/REML; resticted maximum likelihood) solutions for the following mixed model.
#'
#' \deqn{y = X \beta + Z u + \epsilon}
#'
#' where \eqn{\beta} is a vector of fixed effects and \eqn{u} is a vector of random effects with
#' \eqn{Var[u] = K \sigma^2_u}. The residual variance is \eqn{Var[\epsilon] = I \sigma^2_e}.
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param ZETA A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param n.thres If \eqn{n >= n.thres}, perform EMM1.cpp. Else perform EMM2.cpp.
#' @param reestimation If TRUE, EMM2.cpp is performed when the estimation by EMM1.cpp may not be accurate.
#' @param lam.len The number of initial values you set. If this number is large, the estimation will be more accurate,
#' but computational cost will be large. We recommend setting this value 3 <= lam.len <= 6.
#' @param init.range The range of the initial parameters. For example, if lam.len = 5 and init.range = c(1e-06, 1e02),
#' corresponding initial heritabilities will be calculated as seq(1e-06, 1 - 1e-02, length = 5),
#' and then initial lambdas will be set.
#' @param init.one The initial parameter if lam.len = 1.
#' @param conv.param The convergence parameter. If the diffrence of log-likelihood by updating the parameter "lambda"
#' is smaller than this conv.param, the iteration steps will be stopped.
#' @param count.max Sometimes algorithms won't converge for some initial parameters.
#' So if the iteration steps reache to this argument, you can stop the calculation even if algorithm doesn't converge.
#' @param bounds Lower and Upper bounds of the parameter lambda. If the updated parameter goes out of this range,
#' the parameter is reset to the value in this range.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param traceInside Perform trace for the optimzation if traceInside >= 1, and this argument shows the frequency of reports.
#' @param REML You can choose which method you will use, "REML" or "ML".
#' If REML = TRUE, you will perform "REML", and if REML = FALSE, you will perform "ML".
#' @param silent If this argument is TRUE, warning messages will be shown when estimation is not accurate.
#' @param plot.l If you want to plot log-likelihood, please set plot.l = TRUE.
#' We don't recommend plot.l = TRUE when lam.len >= 2.
#' @param SE If TRUE, standard errors are calculated.
#' @param return.Hinv If TRUE, the function returns the inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}. This is useful for GWAS.
#'
#' @return
#' \describe{
#' \item{$Vu}{Estimator for \eqn{\sigma^2_u}}
#' \item{$Ve}{Estimator for \eqn{\sigma^2_e}}
#' \item{$beta}{BLUE(\eqn{\beta})}
#' \item{$u}{BLUP(\eqn{u})}
#' \item{$LL}{Maximized log-likelihood (full or restricted, depending on method)}
#' \item{$beta.SE}{Standard error for \eqn{\beta} (If SE = TRUE)}
#' \item{$u.SE}{Standard error for \eqn{u^*-u} (If SE = TRUE)}
#' \item{$Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I} (If return.Hinv = TRUE)}
#' \item{$Hinv2}{The inverse of \eqn{H2 = ZKZ'/\lambda + I} (If return.Hinv = TRUE)}
#' \item{$lambda}{Estimators for \eqn{\lambda = \sigma^2_e / \sigma^2_u} (If \eqn{n >= n.thres})}
#' \item{$lambdas}{Lambdas for each initial values (If \eqn{n >= n.thres})}
#' \item{$reest}{If parameter estimation may not be accurate, reest = 1, else reest = 0 (If \eqn{n >= n.thres})}
#' \item{$counts}{The number of iterations until convergence for each initial values (If \eqn{n >= n.thres})}
#' }
#'
#'
#' @references Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis
#'  for association studies. Nat Genet. 44(7): 821-824.
#'
#'
#' @example R/examples/EMM.cpp_example.R
#'
#'
#'
#'
EMM.cpp <- function(y, X = NULL, ZETA, eigen.G = NULL, eigen.SGS = NULL, n.thres = 450, reestimation = FALSE,
                    lam.len = 4, init.range = c(1e-06, 1e02), init.one = 0.5, conv.param = 1e-06,
                    count.max = 20, bounds = c(1e-06, 1e06), tol = NULL, optimizer = "nlminb", traceInside = 0, REML = TRUE,
                    silent = TRUE, plot.l = FALSE, SE = FALSE, return.Hinv = TRUE){
  n <- length(y)

  if(n >= n.thres){
    res <- EMM1.cpp(y = y, X = X, ZETA = ZETA, eigen.G = eigen.G, lam.len = lam.len,
                    init.range = init.range, init.one = init.one, conv.param = conv.param,
                    count.max = count.max, bounds = bounds, tol = tol, REML = REML,
                    silent = silent, plot.l = plot.l, SE = SE, return.Hinv = return.Hinv)

    if((res$reest == 1) & reestimation){
      res <- EMM2.cpp(y = y, X = X, ZETA = ZETA, eigen.G = eigen.G, eigen.SGS = eigen.SGS, REML = REML,
                      optimizer = optimizer, traceInside = traceInside, tol = tol, bounds = bounds, SE = SE, return.Hinv = return.Hinv)
    }

  }else{
    res <- EMM2.cpp(y = y, X = X, ZETA = ZETA, eigen.G = eigen.G, eigen.SGS = eigen.SGS, REML = REML,
                    optimizer = optimizer, traceInside = traceInside, bounds = bounds, SE = SE, return.Hinv = return.Hinv, tol = tol)
  }

  return(res)
}





#' Equation of mixed model for multi-kernel (slow, general version)
#'
#' @description This function solves the following multi-kernel linear mixed effects model.
#'
#' \eqn{y = X \beta + \sum _{l=1} ^ {L} Z _ {l} u _ {l} + \epsilon}
#'
#' where \eqn{Var[y] = \sum _{l=1} ^ {L} Z _ {l} K _ {l} Z _ {l}' \sigma _ {l} ^ 2 + I \sigma _ {e} ^ {2}}.
#'
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X0 A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param ZETA A list of variance matrices and its design matrices of random effects. You can use more than one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z.A, K = K.A), D = list(Z = Z.D, K = K.D)) (A for additive, D for dominance)
#' Please set names of lists "Z" and "K"!
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param traceInside Perform trace for the optimzation if traceInside >= 1, and this argument shows the frequency of reports.
#' @param n.thres If \eqn{n >= n.thres}, perform EMM1.cpp. Else perform EMM2.cpp.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param REML You can choose which method you will use, "REML" or "ML".
#' If REML = TRUE, you will perform "REML", and if REML = FALSE, you will perform "ML".
#' @param pred If TRUE, the fitting values of y is returned.
#'
#' @return
#' \describe{
#' \item{$y.pred}{The fitting values of y \eqn{y = X\beta + Zu}}
#' \item{$Vu}{Estimator for \eqn{\sigma^2_u}, all of the genetic variance}
#' \item{$Ve}{Estimator for \eqn{\sigma^2_e}}
#' \item{$beta}{BLUE(\eqn{\beta})}
#' \item{$u}{BLUP(\eqn{u})}
#' \item{$weights}{The proportion of each genetic variance (corresponding to each kernel of ZETA) to Vu}
#' \item{$LL}{Maximized log-likelihood (full or restricted, depending on method)}
#' \item{$Vinv}{The inverse of \eqn{V = Vu \times ZKZ' + Ve \times I}}
#' \item{$Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I}}
#' }
#'
#' @references Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis
#'  for association studies. Nat Genet. 44(7): 821-824.
#'
#'
#' @example R/examples/EM3.cpp_example.R
#'
#'
#'
#'
EM3.cpp <- function (y, X0 = NULL, ZETA, eigen.G = NULL, eigen.SGS = NULL, tol = NULL,
                     optimizer = "nlminb", traceInside = 0, n.thres = 450, REML = TRUE, pred = TRUE){
  n <- length(as.matrix(y))
  y <- matrix(y, n, 1)

  not.NA <- which(!is.na(y))

  if (is.null(X0)) {
    p <- 1
    X0 <- matrix(rep(1, n), n, 1)
  }
  p <- ncol(X0)

  lz <- length(ZETA)
  weights <- rep(1/lz, lz)

  if(lz >= 2){
    if(is.null(eigen.G)){
      Z <- c()
      ms <- rep(NA, lz)
      for (i in 1:lz) {
        Z.now <- ZETA[[i]]$Z

        if(is.null(Z.now)){
          Z.now <- diag(n)
        }
        m <- ncol(Z.now)
        if (is.null(m)) {
          m <- 1
          Z.now <- matrix(Z.now, length(Z.now), 1)
        }

        K.now <- ZETA[[i]]$K

        if (!is.null(K.now)) {
          stopifnot(nrow(K.now) == m)
          stopifnot(ncol(K.now) == m)
        }

        Z <- cbind(Z, Z.now)
        ms[i] <- m
      }

      stopifnot(nrow(Z) == n)
      stopifnot(nrow(X0) == n)

      Z <- Z[not.NA, , drop = FALSE]
      X <- X0[not.NA, , drop = FALSE]
      n <- length(not.NA)
      y <- matrix(y[not.NA], n, 1)

      spI <- diag(n)
      S <- spI - tcrossprod(X %*% solve(crossprod(X)), X)


      minimfunctionouter <- function(weights = rep(1 / lz, lz)) {
        weights <- weights / sum(weights)
        ZKZt <- matrix(0, nrow = n, ncol = n)
        for (i in 1:lz) {
          ZKZt <- ZKZt + weights[i] * tcrossprod(ZETA[[i]]$Z[not.NA, ] %*%
                                                   ZETA[[i]]$K, ZETA[[i]]$Z[not.NA, ])
        }


        res <- EM3_kernel(y, X, ZKZt, S, spI, n, p)
        lambda <- res$lambda
        eta <- res$eta
        phi <- res$phi

        if(REML){
          minimfunc <- function(delta) {
            (n - p) * log(sum(eta ^ 2/{
              lambda + delta
            })) + sum(log(lambda + delta))
          }
        }else{
          minimfunc <- function(delta) {
            n * log(sum(eta ^ 2/{
              lambda + delta
            })) + sum(log(phi + delta))
          }
        }

        optimout <- optimize(minimfunc, lower = 0, upper = 10000)
        return(optimout$objective)
      }


      parInit <- rep(1 / lz, lz)
      parLower <- rep(0, lz)
      parUpper <- rep(1, lz)
      traceNo <- ifelse(traceInside > 0, 3, 0)
      traceREPORT <- ifelse(traceInside > 0, traceInside, 1)

      if (optimizer == "optim") {
        soln <- optim(par = parInit, fn = minimfunctionouter, control = list(trace = traceNo, REPORT = traceREPORT),
                      method = "L-BFGS-B", lower = parLower, upper = parUpper)
        weights <- soln$par
      } else if (optimizer == "optimx") {
        soln <- optimx::optimx(par = parInit, fn = minimfunctionouter, gr = NULL, hess = NULL,
                               lower = parLower, upper = parUpper, method = "L-BFGS-B", hessian = FALSE,
                               control = list(trace = traceNo, starttests = FALSE, maximize = F, REPORT = traceREPORT, kkt = FALSE))
        weights <- as.matrix(soln)[1, 1:lz]
      } else if (optimizer == "nlminb") {
        soln <- nlminb(start = parInit, objective = minimfunctionouter, gradient = NULL, hessian = NULL,
                       lower = parLower, upper = parUpper, control = list(trace = traceInside))
        weights <- soln$par
      } else {
        warning("We offer 'optim', 'optimx', and 'nlminb' as optimzers. Here we use 'optim' instead.")
        soln <- optim(par = parInit, fn = minimfunctionouter, control = list(trace = traceNo, REPORT = traceREPORT),
                      method = "L-BFGS-B", lower = parLower, upper = parUpper)
        weights <- soln$par
      }

    }
  }


  Z <- c()
  ms <- rep(NA, lz)
  for (i in 1:lz) {
    Z.now <- ZETA[[i]]$Z

    if(is.null(Z.now)){
      Z.now <- diag(n)
    }
    m <- ncol(Z.now)
    if (is.null(m)) {
      m <- 1
      Z.now <- matrix(Z.now, length(Z.now), 1)
    }

    K.now <- ZETA[[i]]$K

    if (!is.null(K.now)) {
      stopifnot(nrow(K.now) == m)
      stopifnot(ncol(K.now) == m)
    }

    Z <- cbind(Z, Z.now)
    ms[i] <- m
  }

  Z <- Z[not.NA, , drop = FALSE]
  if((lz == 1) | (!is.null(eigen.G))){
    X <- X0[not.NA, , drop = FALSE]
    n <- length(not.NA)
    y <- matrix(y[not.NA], n, 1)
    spI <- diag(n)
  }

  weights <- weights / sum(weights)

  ZKZt <- matrix(0, nrow = n, ncol = n)
  Klist <- NULL
  for(i in 1:lz){
    K.list <- list(ZETA[[i]]$K)
    Klist <- c(Klist, K.list)
  }
  Klistweighted <- Klist
  for (i in 1:lz) {
    Klistweighted[[i]] <- weights[i] * ZETA[[i]]$K
    ZKZt <- ZKZt + weights[i] *
      tcrossprod(as.matrix(ZETA[[i]]$Z[not.NA, ]) %*% ZETA[[i]]$K, ZETA[[i]]$Z[not.NA, ])
  }
  K <- Matrix::.bdiag(Klistweighted)
  ZK <- as.matrix(Z %*% K)



  if(is.null(eigen.G)){
    if(nrow(Z) <= n.thres){
      # return.SGS <- TRUE
      return.SGS <-  FALSE
    }else{
      return.SGS <- FALSE
    }
    spectralG.res <- spectralG.cpp(ZETA = lapply(ZETA, function(x) list(Z = x$Z[not.NA, , drop = FALSE], K = x$K)),
                                   X = X, weights = weights, return.G = TRUE, return.SGS = return.SGS,
                                   tol = tol, df.H = NULL)

    eigen.G <- spectralG.res[[1]]
    eigen.SGS <- spectralG.res[[2]]
  }

  if(lz >= 2){
    ZETA.list <- list(A = list(Z = diag(n), K = ZKZt))
    EMM.cpp.res <- EMM.cpp(y, X = X, ZETA = ZETA.list, eigen.G = eigen.G, eigen.SGS = eigen.SGS, traceInside = traceInside,
                           optimizer = optimizer, tol = tol, n.thres = n.thres, return.Hinv = TRUE, REML = REML)
  }else{
    EMM.cpp.res <- EMM.cpp(y, X = X, ZETA = lapply(ZETA, function(x) list(Z = x$Z[not.NA, , drop = FALSE], K = x$K)), traceInside = traceInside,
                           optimizer = optimizer, eigen.G = eigen.G, eigen.SGS = eigen.SGS, tol = tol, n.thres = n.thres, return.Hinv = TRUE, REML = REML)
  }


  Vu <- EMM.cpp.res$Vu
  Ve <- EMM.cpp.res$Ve
  beta <- EMM.cpp.res$beta
  LL <- EMM.cpp.res$LL
  Hinv <- EMM.cpp.res$Hinv


  e <- (y - {X %*% beta})
  Hinve <- Hinv %*% e

  u <- crossprod(ZK, Hinve)
  Vinv <- (1 / Ve) * Hinv
  namesu <- c()
  for (i in 1:length(Klist)) {
    namesu <- c(namesu, paste("K", i, colnames(ZETA[[i]]$K), sep = "_"))
  }

  rownames(u) <- namesu
  if(!pred){
    return(list(y.pred = NULL, Vu = Vu, Ve = Ve, beta = beta,
                u = u, weights = weights, LL = LL, Vinv = Vinv, Hinv = Hinv))
  }else{
    n.all <- nrow(ZETA[[1]]$Z)

    Z.all <- c()
    for (i in 1:lz) {
      if(is.null(Z.now)){
        Z.now <- diag(n.all)
      }else{
        Z.now <- ZETA[[i]]$Z
      }
      Z.all <- cbind(Z.all, Z.now)
    }


    if(ncol(X) != 1){
      y.pred <- c(X0 %*% beta) + c(Z.all %*% u)
    }else{
      y.pred <- rep(beta, n.all) + c(Z.all %*% u)
    }

    return(list(y.pred = y.pred, Vu = Vu, Ve = Ve, beta = beta,
                u = u, weights = weights, LL = LL, Vinv = Vinv, Hinv = Hinv))
  }
}





#' Equation of mixed model for multi-kernel (fast, for limited cases)
#'
#' @description This function solves multi-kernel mixed model using fastlmm.snpset approach (Lippert et al., 2014).
#' This function can be used only when the kernels other than genomic relationship matrix are linear kernels.
#'
#'
#' @param y0 A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X0 A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param ZETA A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param Zs0 A list of design matrices (Z; \eqn{n \times m} matrix) for Ws.
#' For example, Zs0 = list(A.part = Z.A.part, D.part = Z.D.part)
#' @param Ws0 A list of low rank matrices (W; \eqn{m \times k} matrix). This forms linear kernel \eqn{K = W \Gamma W'}.
#' For example, Ws0 = list(A.part = W.A, D.part = W.D)
#' @param Gammas0 A list of matrices for weighting SNPs (Gamma; \eqn{k \times k} matrix). This forms linear kernel \eqn{K = W \Gamma W'}.
#' For example, if there is no weighting, Gammas0 = lapply(Ws0, function(x) diag(ncol(x)))
#' @param gammas.diag If each Gamma is the diagonal matrix, please set this argument TRUE. The calculationtime can be saved.
#' @param X.fix If you repeat this function and when X0 is fixed during iterations, please set this argument TRUE.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param tol The tolerance for detecting linear dependencies in the columns of G = ZKZ'.
#' Eigen vectors whose eigen values are less than "tol" argument will be omitted from results.
#' If tol is NULL, top 'n' eigen values will be effective.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param traceInside Perform trace for the optimzation if traceInside >= 1, and this argument shows the frequency of reports.
#' @param n.thres If \eqn{n >= n.thres}, perform EMM1.cpp. Else perform EMM2.cpp.
#' @param bounds Lower and upper bounds for weights.
#' @param spectral.method The method of spectral decomposition.
#' In this function, "eigen" : eigen decomposition and "cholesky" : cholesky and singular value decomposition are offered.
#' If this argument is NULL, either method will be chosen accorsing to the dimension of Z and X.
#' @param REML You can choose which method you will use, "REML" or "ML".
#' If REML = TRUE, you will perform "REML", and if REML = FALSE, you will perform "ML".
#' @param pred If TRUE, the fitting values of y is returned.
#'
#' @return
#' \describe{
#' \item{$y.pred}{The fitting values of y \eqn{y = X\beta + Zu}}
#' \item{$Vu}{Estimator for \eqn{\sigma^2_u}, all of the genetic variance}
#' \item{$Ve}{Estimator for \eqn{\sigma^2_e}}
#' \item{$beta}{BLUE(\eqn{\beta})}
#' \item{$u}{BLUP(\eqn{u})}
#' \item{$weights}{The proportion of each genetic variance (corresponding to each kernel of ZETA) to Vu}
#' \item{$LL}{Maximized log-likelihood (full or restricted, depending on method)}
#' \item{$Vinv}{The inverse of \eqn{V = Vu \times ZKZ' + Ve \times I}}
#' \item{$Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I}}
#' }
#'
#' @references Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis
#'  for association studies. Nat Genet. 44(7): 821-824.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#' @example R/examples/EM3.linker.cpp_example.R
#'
#'
#'
EM3.linker.cpp <- function (y0, X0 = NULL, ZETA = NULL, Zs0 = NULL, Ws0,
                            Gammas0 = lapply(Ws0, function(x) diag(ncol(x))), gammas.diag = TRUE,
                            X.fix = TRUE, eigen.SGS = NULL, eigen.G = NULL, tol = NULL,
                            bounds = c(1e-06, 1e06), optimizer = "nlminb", traceInside = 0,
                            n.thres = 450, spectral.method = NULL, REML = TRUE, pred = TRUE){
  n0 <- length(as.matrix(y0))
  y0 <- matrix(y0, n0, 1)

  if(is.null(ZETA)){
    Z0 <- diag(n0)
    K <- diag(n0)
  }else{
    if(length(ZETA) != 1){
      stop("The kernel for cofounders should be one!")
    }
    Z0 <- ZETA[[1]]$Z
    K <- ZETA[[1]]$K
  }


  if (is.null(X0)) {
    p <- 1
    X0 <- matrix(rep(1, n0), n0, 1)
  }
  p <- ncol(X0)


  not.NA <- which(!is.na(y0))
  n <- length(not.NA)
  y <- y0[not.NA]
  X <- X0[not.NA, , drop = FALSE]
  Z <- Z0[not.NA, , drop = FALSE]

  ZETA <- list(A = list(Z = Z, K = K))

  if(is.null(Zs0)){
    Zs0 <- lapply(Ws0, function(x) diag(nrow(x)))
  }
  Zs <- lapply(Zs0, function(x) x[not.NA, ])

  Ws <- NULL
  for(i in 1:length(Ws0)){
    Ws.part <- list(Zs[[i]] %*% Ws0[[i]])
    Ws <- c(Ws, Ws.part)
  }
  lw <- length(Ws)

  offset <- sqrt(n)
  spI <- diag(n)


  if(!REML){
    if(is.null(eigen.G)){
      eigen.G <- spectralG.cpp(ZETA = ZETA, X = X, return.G = TRUE,
                               return.SGS = FALSE, tol = tol, df.H = NULL)[[1]]
    }
    U <- eigen.G$U
    delta <- eigen.G$delta
  }


  if((!X.fix) | (is.null(eigen.SGS))){
    eigen.SGS <- spectralG.cpp(ZETA = ZETA, X = X, return.G = FALSE,
                               return.SGS = TRUE, tol = tol, df.H = NULL)[[2]]
  }
  U2 <- eigen.SGS$Q
  delta2 <- eigen.SGS$theta


  n.param <- lw + 1
  weights <- rep(1 / n.param, n.param + 1)
  weights[2:(n.param + 1)] <- weights[2:(n.param + 1)] / sum(weights[2:(n.param + 1)])

  minimumfunctionouter <- function(weights){
    weights[2:(n.param + 1)] <- weights[2:(n.param + 1)] / sum(weights[2:(n.param + 1)])

    lambda <- weights[1] / weights[2]
    gammas <- weights[-(1:2)] / weights[2]

    Gammas <- NULL
    for(i in 1:lw){
      Gammas[[i]] <- Gammas0[[i]] * gammas[i]
    }
    Gammas <- lapply(Gammas, as.matrix)
    P.res <- P_calc(lambda, Ws, Gammas, U2, as.matrix(delta2), lw, TRUE, gammas.diag)
    P <- P.res$P
    lnP <- P.res$lnP

    yPy <- c(crossprod(y, P %*% y))

    if(REML){
      LL <- llik_REML(n, p, yPy, lnP)
    }else{
      Hinv.res <- P_calc(lambda, Ws, Gammas, U, as.matrix(delta), lw, FALSE, gammas.diag)
      Hinv <- Hinv.res$P
      lnHinv <- Hinv.res$lnP

      LL <- llik_ML(n, yPy, lnHinv)
    }

    obj <- -LL
    return(obj)
  }


  parInit <- rep(1 / n.param, n.param + 1)
  parLower <- rep(bounds[1], n.param + 1)
  parUpper <- c(bounds[2], rep(1, n.param))
  traceNo <- ifelse(traceInside > 0, 3, 0)
  traceREPORT <- ifelse(traceInside > 0, traceInside, 1)

  if (optimizer == "optim") {
    soln <- optim(par = parInit, fn = minimumfunctionouter, control = list(trace = traceNo, REPORT = traceREPORT),
                  method = "L-BFGS-B", lower = parLower, upper = parUpper)
    weights <- soln$par[-1]
  } else if (optimizer == "optimx") {
    soln <- optimx::optimx(par = parInit, fn = minimumfunctionouter, gr = NULL, hess = NULL,
                           lower = parLower, upper = parUpper, method = "L-BFGS-B", hessian = FALSE,
                           control = list(trace = traceNo, starttests = FALSE, maximize = F, REPORT = traceREPORT, kkt = FALSE))
    weights <- as.matrix(soln)[1, 2:(n.param + 1)]
  } else if (optimizer == "nlminb") {
    soln <- nlminb(start = parInit, objective = minimumfunctionouter, gradient = NULL, hessian = NULL,
                   lower = parLower, upper = parUpper, control = list(trace = traceInside))
    weights <- soln$par[-1]
  } else {
    warning("We offer 'optim', 'optimx', and 'nlminb' as optimzers. Here we use 'optim' instead.")
    soln <- optim(par = parInit, fn = minimumfunctionouter, control = list(trace = traceNo, REPORT = traceREPORT),
                  method = "L-BFGS-B", lower = parLower, upper = parUpper)
    weights <- soln$par[-1]
  }


  Zs.all.list <- c(list(K.A = Z), Zs)
  Zs0.all.list <- c(list(K.A = Z0), Zs0)
  ZWs <- NULL

  Zs.all <- Z
  Zs0.all <- Z0
  weights <- weights / sum(weights)
  ZKZt <- weights[1] * tcrossprod(Z %*% K, Z)
  Klist <-  list(K.A = K)
  Klistweighted <- list(K.A = weights[1] * K)
  for(i in 1:lw){
    Z.now <- Zs[[i]]
    Z0.now <- Zs0[[i]]
    K.now <- tcrossprod(Ws0[[i]] %*% Gammas0[[i]], Ws0[[i]])
    Kweighted.now <- weights[i + 1] * tcrossprod(Ws0[[i]] %*% Gammas0[[i]], Ws0[[i]])
    ZKZt.now <-  weights[i + 1] * tcrossprod(Ws[[i]] %*% Gammas0[[i]], Ws[[i]])

    Zs.all <- cbind(Zs.all, Z.now)
    Zs0.all <- cbind(Zs0.all, Z0.now)

    ZWs.now <- list(list(Z = Z.now, W = Ws0[[i]], Gamma = Gammas0[[i]]))
    names(ZWs.now) <- names(Zs)[i]
    ZWs <- c(ZWs, ZWs.now)

    Klist <- c(Klist, list(K.now))
    Klistweighted <- c(Klistweighted, list(Kweighted.now))
    ZKZt <- ZKZt + ZKZt.now
  }

  K.all <- Matrix::.bdiag(Klistweighted)
  ZK <- as.matrix(Zs.all %*% K.all)


  ZETA.list <- list(A = list(Z = diag(n), K = ZKZt))
  spectralG.all <- spectralG.cpp(ZETA = ZETA, ZWs = ZWs, weights = weights, X = X, return.G = TRUE,
                                 return.SGS = TRUE, tol = tol, df.H = NULL, spectral.method = "eigen")
  eigen.G.all <- spectralG.all[[1]]
  eigen.SGS.all <- spectralG.all[[2]]

  EMM.cpp.res <- EMM.cpp(y, X = X, ZETA = ZETA.list, eigen.G = eigen.G.all, eigen.SGS = eigen.SGS.all,
                         optimizer = optimizer, traceInside = traceInside, n.thres = n.thres, return.Hinv = TRUE, REML = REML, tol = tol)

  Vu <- EMM.cpp.res$Vu
  Ve <- EMM.cpp.res$Ve
  beta <- EMM.cpp.res$beta
  LL <- EMM.cpp.res$LL
  Hinv <- EMM.cpp.res$Hinv

  e <- (y - {X %*% as.matrix(beta)})
  Hinve <- Hinv %*% e

  u <- crossprod(ZK, Hinve)
  Vinv <- (1 / Ve) * Hinv
  namesu <- c()
  for (i in 1:length(Klist)) {
    namesu <- c(namesu, paste("K", i, colnames(Klist[[i]]), sep = "_"))
  }

  rownames(u) <- namesu
  if(!pred){
    return(list(y.pred = NULL, Vu = Vu, Ve = Ve, beta = beta,
                u = u, weights = weights, LL = LL, Vinv = Vinv, Hinv = Hinv))
  }else{
    if(ncol(X) != 1){
      y.pred <- c(X0 %*% beta) + c(Zs0.all %*% u)
    }else{
      y.pred <- rep(beta, n0) + c(Zs0.all %*% u)
    }

    return(list(y.pred = y.pred, Vu = Vu, Ve = Ve, beta = beta,
                u = u, weights = weights, LL = LL, Vinv = Vinv, Hinv = Hinv))
  }
}



#' Calculte -log10(p) by score test (slow, for general cases)
#'
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param Gs A list of kernel matrices you want to test. For example, Gs = list(A.part = K.A.part, D.part = K.D.part)
#' @param Gu A \eqn{n \times n} matrix. You should assign \eqn{ZKZ'}, where K is covariance (relationship) matrix and Z is its design matrix.
#' @param Ge A \eqn{n \times n} matrix. You should assign identity matrix I (diag(n)).
#' @param P0 A \eqn{n \times n} matrix. The Moore-Penrose generalized inverse of \eqn{SV0S}, where \eqn{S = X(X'X)^{-1}X'} and
#' \eqn{V0 = \sigma^2_u Gu + \sigma^2_e Ge}. \eqn{\sigma^2_u} and \eqn{\sigma^2_e} are estimators of the null model.
#' @param chi0.mixture RAINBOW assumes the test statistic \eqn{l1' F l1} is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where l1 is the first derivative of the log-likelihood and F is the Fisher information. And r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#'
#' @return -log10(p) calculated by score test
#'
#'
#'
#'
score.cpp <- function(y, Gs, Gu, Ge, P0, chi0.mixture = 0.5){
  nuisance.no <- 2
  Gs.all <- c(Gs, list(Gu), list(Ge))
  
  l1 <- score_l1(y = as.matrix(y), p0 = P0, Gs = Gs, lg = length(Gs))
  F.info <- score_fisher(p0 = P0, Gs_all = Gs.all, nuisance_no = nuisance.no,
                         lg_all = length(Gs.all))
  
  
  score.stat <- c(crossprod(l1, F.info %*% l1))
  logp <- ifelse(score.stat <= 0, 0, -log10((1 - chi0.mixture) *
                                              pchisq(deviance, df = length(Gs), lower.tail = FALSE)))
  return(logp)
}



#' Calculte -log10(p) by score test (fast, for limited cases)
#'
#'
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param Ws A list of low rank matrices (ZW; \eqn{n \times k} matrix). This forms linear kernel \eqn{ZKZ' = ZW \Gamma (ZW)'}.
#' For example, Ws = list(A.part = ZW.A, D.part = ZW.D)
#' @param Gammas A list of matrices for weighting SNPs (Gamma; \eqn{k \times k} matrix). This forms linear kernel \eqn{ZKZ' = ZW \Gamma (ZW)'}.
#' For example, if there is no weighting, Gammas = lapply(Ws, function(x) diag(ncol(x)))
#' @param gammas.diag If each Gamma is the diagonal matrix, please set this argument TRUE. The calculation time can be saved.
#' @param Gu A \eqn{n \times n} matrix. You should assign \eqn{ZKZ'}, where K is covariance (relationship) matrix and Z is its design matrix.
#' @param Ge A \eqn{n \times n} matrix. You should assign identity matrix I (diag(n)).
#' @param P0 A \eqn{n \times n} matrix. The Moore-Penrose generalized inverse of \eqn{SV0S}, where \eqn{S = X(X'X)^{-1}X'} and
#' \eqn{V0 = \sigma^2_u Gu + \sigma^2_e Ge}. \eqn{\sigma^2_u} and \eqn{\sigma^2_e} are estimators of the null model.
#' @param chi0.mixture RAINBOW assumes the statistic \eqn{l1' F l1} follows the mixture of \eqn{\chi^2_0} and \eqn{\chi^2_r},
#' where l1 is the first derivative of the log-likelihood and F is the Fisher information. And r is the degree of freedom.
#' chi0.mixture determins the proportion of \eqn{\chi^2_0}
#' @return -log10(p) calculated by score test
#'
#'
#'
#'
score.linker.cpp <- function(y, Ws, Gammas, gammas.diag = TRUE, Gu, Ge, P0, chi0.mixture = 0.5){
  nuisance.no <- 2
  if(gammas.diag){
    W2s <- NULL
    for(i in 1:length(Ws)){
      W2s.now <- t(t(Ws[[i]]) * sqrt(diag(Gammas[[i]])))
      W2s <- c(W2s, list(W2s.now))
    }
    names(W2s) <- names(Ws)
    Gs.all <- c(W2s, list(Gu), list(Ge))
    
    l1 <- score_l1_linker_diag(y = as.matrix(y), p0 = P0, W2s = W2s, lw = length(W2s))
    F.info <- score_fisher_linker_diag(p0 = P0, Gs_all = Gs.all,
                                       nuisance_no = nuisance.no, lw_all = length(Gs.all))
  }else{
    Gs.all <- c(Ws, list(Gu), list(Ge))
    
    l1 <- score_l1_linker(y = as.matrix(y), p0 = P0, Ws = Ws, Gammas = Gammas, lw = length(Ws))
    F.info <- score_fisher_linker(p0 = P0, Gs_all = Gs.all, Gammas = Gammas,
                                  nuisance_no = nuisance.no, lw_all = length(Gs.all))
  }
  score.stat <- c(crossprod(l1, F.info %*% l1))
  logp <- ifelse(score.stat <= 0, 0, -log10((1 - chi0.mixture) *
                                              pchisq(score.stat, df = length(Ws), lower.tail = FALSE)))
  return(logp)
}
