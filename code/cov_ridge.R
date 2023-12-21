vcov.ridge <- function(model)
{ # return the covariance matrix (and its determinant) for the ridge estimator
  obj <- model
  x <- model.matrix(obj$terms, obj$model, obj$contrast)
  s2 <- obj$scale
  lambda <- obj$lambda

  rs <- svd(x, nu = 0)
  v  <- rs$v
  d2 <- rs$d^2
  p  <- length(d2)

  # computing covariance and its determinant
  mid <- d2 / (d2 + lambda)^2
  COV <- s2 * (v %*% diag(mid) %*% t(v))
  COV <- sqrt(diag(COV))
  DET <- (s2^p) * prod(mid)
  
  attr(COV, 'determinant') <- DET
  COV
}

