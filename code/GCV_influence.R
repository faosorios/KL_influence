## ID: lambda_influence.R, last updated 2023-12-21, F.Osorio

GCV.influence <- function(object, scheme = "response")
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
  edf <- object$edf
  RSS <- object$RSS

  # computing H(lambda)
  div <- rs$d^2 + lambda
  Delta1 <- rs$d^2 / div
  u <- rs$u %*% diag(sqrt(Delta1))
  H <- tcrossprod(u)
  
  # computing M(lambda)
  Delta2 <- (rs$d / div)^2
  u <- rs$u %*% diag(sqrt(Delta2))
  M <- tcrossprod(u)
  tr2 <- sum(Delta2) # tr(M)

  switch(scheme,
          "scale" = { # hmax under perturbation of variances
          id <- diag(n)
          v1 <- c((id - 2 * H) %*% res) * c(H %*% res)
          v2 <- c(y + 4 * res) * c(H %*% (id - H) %*% res)
          v3 <- c((3 * id - 2 * H) %*% y) * c(H %*% res)
          v4 <- c((2 * id - H) %*% y) * res
          v5 <- diag(H %*% (id - H) %*% (id - 2 * H))
          v6 <- diag(H %*% (id - H))
          d  <- n - edf
          a1 <- sum(v6)
          a2 <- sum(res^2)
          a3 <- c(crossprod(res, H %*% res))
          t1 <- v1 + v2 + 2 * (a1 / d) * (v3 - v4)
          t2 <- 2 * (a2 / d) * v5
          t3 <- 4 * (a3 / d) * v6
          t4 <- 6 * (a1 * a2 / d^2) * v6
          hmax <- t1 - t2 + t3 - t4
          hmax <- as.vector(hmax)
         },
          "response" = { # hmax under response perturbation
          id <- diag(n)
          a1 <- sum(diag(H %*% (id - H))) / (n - edf)
          res <- (id - H) %*% y
          M <- (a1 * id - H) %*% (id - H)
          hmax <- M %*% res
          hmax <- as.vector(hmax)
         },
         stop("scheme =", scheme, " in not implemented.")
  )
  
  # compute largest eigenvectors and create the output object
  hmax <- abs(hmax)
  hmax <- hmax / minkowski(hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)

  out <- list(hmax = hmax, cutoff = cutoff)
  out
}
