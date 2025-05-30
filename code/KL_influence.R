## ID: KL_influence.R, last updated 2025-03-31, F.Osorio

KL.divergence.ridge <- function(model) 
{ ## Kullback-Liebler divergence for ridge regression
  obj <- model
  y <- model.response(obj$model, "numeric")
  x <- model.matrix(obj$terms, obj$model, obj$contrast)
  n <- nrow(x)
  p <- ncol(x)

  # SVD of model matrix
  rs <- svd(x)

  # extracting estimates
  cf <- obj$coef
  s2 <- obj$scale
  lambda <- obj$lambda

  div <- rs$d^2 + lambda
  a <- crossprod(rs$v, cf)
  z <- (rs$d * a) / div

  # leverages (code re-use)
  Delta <- rs$d^2 / div 
  u <- rs$u %*% diag(sqrt(Delta))
  levs <- rowSums(u^2)
  hats <- rowSums(rs$u^2) # leverages form LS estimation

  # computing KL-divergence
  u <- rs$u %*% diag(Delta)
  D2  <- rowSums(u^2)
  qi  <- c(rs$u %*% z)
  rel <- hats / (1 - hats)
  KL  <- .5 * (lambda * qi^2 / s2 + 1 + D2) * rel
  rel <- levs / (1 - hats)
  KL  <- KL + rel - log(1 - levs) + .5 * log(1 - hats)
  KL
}

transition <- function(n = 1)
{ # transition matrix of order 'n'
  mat <- matrix(0, nrow = n * n, ncol = n)
  for (j in 1:n) {
    row <- (j - 1) * n + j
    mat[row,j] <- 1 
  }
  mat
}

curvature.KL.ridge <- function(object)
{ 
  y <- model.response(object$model, "numeric")
  x <- model.matrix(object$terms, object$model, object$contrast)
  n <- nrow(x)
  p <- ncol(x)

  # SVD of model matrix
  rs <- svd(x, nv = 0)

  # extracting ridge elements
  lambda <- object$lambda
  res <- object$residuals
  yfit <- object$fitted
  s2 <- object$scale

  # computing H and H(lambda)
  div <- rs$d^2 + lambda
  Delta <- rs$d^2 / div
  u <- rs$u %*% diag(sqrt(Delta))
  H  <- tcrossprod(u)
  H0 <- tcrossprod(rs$u)

  # computing Q matrix
  Q <- outer(y, y) - outer(yfit, yfit)

  # curvature under scale perturbation
  term1 <- kronecker.prod(H, H)
  term2 <- kronecker.prod(H0, H %*% Q %*% H) / s2
  term3 <- .5 * kronecker.prod(H0, H0) 
  term4 <- 2. * kronecker.prod(H0, H)
  term5 <- .5 * kronecker.prod(H0, H %*% H)

  B <- transition(n = n)
  curv <- (term1 - term2 + term3 - term4 + term5) %*% B
  curv <- crossprod(B, curv)
  curv <- asSymmetric(curv) # just a safeguard
  scaling <- sqrt(sum(diag(crossprod(curv))))
  curv <- curv / scaling
  
  # compute largest eigenvectors and create the output object
  rs <- svd(curv, nv = 0)
  which <- abs(rs$d) < .Machine$double.eps
  hmax <- rs$u[,1]

  out <- list(hmax = abs(hmax), vectors = rs$u[,!which])
  out
}

explanatory.KL.ridge <- function(object, which = 2)
{ 
  y <- model.response(object$model, "numeric")
  x <- model.matrix(object$terms, object$model, object$contrast)
  n <- nrow(x)
  p <- ncol(x)
  k <- which

  # scales 
  scales <- equilibrate(x)
  scales <- attr(scales, "scales")
  a2 <- scales[k]^2

  # SVD of model matrix
  rs <- svd(x, nv = 0)

  # extracting ridge elements
  lambda <- object$lambda
  res <- object$residuals
  yfit <- object$fitted
  s2 <- object$scale

  # computing H and H(lambda)
  div <- rs$d^2 + lambda
  Delta <- rs$d^2 / div
  u <- rs$u %*% diag(sqrt(Delta))
  H  <- tcrossprod(u)
  H0 <- tcrossprod(rs$u)

  Id <- diag(n)
  xx <- crossprod(x)
  S <- xx + diag(lambda, p)
  R <- solve(xx)
  L <- solve(S, xx) %*% solve(S)
  m <- solve(S, crossprod(x, yfit))
  P <- outer(yfit, yfit)
  Q <- Id - .5 * H
  U <- x %*% R
  V <- x %*% solve(S)
  E0 <- matrix(0, nrow = p, ncol = p)
  E0[k,k] <- 1

  # curvature under explanatory variable perturbation
  term1 <- R[k,k] * P - 2 * (R[k,k] * H + U %*% E0 %*% t(V)) %*% P %*% Q
  term1 <- (a2 / s2) * term1
  term2 <- (L[k,k] + m[k]^2 / s2 + R[k,k]^2) * H0
  term2 <- term2 - 2 * R[k,k] * H %*% Q  + 2 * V %*% E0 %*% t(V)
  term2 <- term2 + U %*% E0 %*% (t(U) - 4 * t(V) %*% Q)
  term2 <- a2 * term2
  curv <- term1 + term2
  curv <- asSymmetric(curv) # just a safeguard
  scaling <- sqrt(sum(diag(crossprod(curv))))
  curv <- curv / scaling
  
  # compute largest eigenvectors and create the output object
  rs <- svd(curv, nv = 0)
  which <- abs(rs$d) < .Machine$double.eps
  hmax <- rs$u[,1]

  out <- list(hmax = abs(hmax), vectors = rs$u[,!which])
  out
}
