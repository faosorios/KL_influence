## ID: ridge_par.R, last updated 2023-12-21, F.Osorio

ridge.par <- function(object)
{ 
  nobs <- object$dims[1]
  form <- object$call$formula
  data <- object$model
  lambda0 <- object$lambda

  lambda  <- rep(0, nobs)
  for (i in 1:nobs) {
    rmdata <- data[-i,]
    lambda[i] <- ridge(form, data = rmdata, lambda = lambda0, method = "GCV")$lambda
  }
  
  lambda
}
