## LD_influence.R

curvature.LD.ridge <- function(object, scheme = "scale")
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
  s2 <- object$scale

  # computing H(lambda)
  div <- rs$d^2 + lambda
  Delta <- rs$d^2 / div
  u <- rs$u %*% diag(sqrt(Delta))
  H <- tcrossprod(u)

  switch(scheme,
          "scale" = { # curvature under scale perturbation
          curv <- H + outer(res, res) / (2 * n * s2)
          curv <- diag(res) %*% curv %*% diag(res)
          curv <- curv / s2
          scaling <- sqrt(sum(diag(crossprod(curv))))
          curv <- curv / scaling
         },
          "response" = { # curvature under response perturbation
          curv <- H + 2 * outer(res, res) / (n * s2)
          curv <- curv / s2
          scaling <- sqrt(sum(diag(crossprod(curv))))
          curv <- curv / scaling
         },
         stop("scheme =", scheme, " in not implemented.")
  )
  
  # compute largest eigenvectors and create the output object
  rs <- svd(curv, nv = 0)
  which <- abs(rs$d) < .Machine$double.eps
  hmax <- rs$u[,1]

  out <- list(hmax = abs(hmax), vectors = rs$u[,!which])
  out
}
