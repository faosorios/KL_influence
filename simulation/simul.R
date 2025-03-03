## Id: simul.R, last updated 2025-01-27
## Author: Felipe Osorio

summary.KL <- function(Nsize = 1000, nobs = 20, k = 1) {
  out <- matrix(0, nrow = 9, ncol = 6)
  now <- proc.time()

  # Scenario I:
  alpha <- c(1,1,1,1,1)
  cat(" 1/9:\n")
  set.seed(5)
  out[1,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 1, alpha = alpha)$percentage
  cat(" 2/9:\n")
  set.seed(23)
  out[2,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 5, alpha = alpha)$percentage
  cat(" 3/9:\n")
  set.seed(53)
  out[3,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 10, alpha = alpha)$percentage

  # Scenario II:
  alpha <- c(0,0,0,1,1)
  cat(" 4/9:\n")
  set.seed(167)
  out[4,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 1, alpha = alpha)$percentage
  cat(" 5/9:\n")
  set.seed(239)
  out[5,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 5, alpha = alpha)$percentage
  cat(" 6/9:\n")
  set.seed(347)
  out[6,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 10, alpha = alpha)$percentage

  # Scenario III:
  alpha <- c(1,0,0,0,1)
  cat(" 7/9:\n")
  set.seed(433)
  out[7,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 1, alpha = alpha)$percentage
  cat(" 8/9:\n")
  set.seed(577)
  out[8,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 5, alpha = alpha)$percentage
  cat(" 9/9:\n")
  set.seed(863)
  out[9,] <- simul.KL(Nsize = Nsize, nobs = nobs, k = k, d = 10, alpha = alpha)$percentage

  colnames(out) <- c("cooks","KL","LD.var","LD.res","KL.var","KL.res")
  rownames(out) <- c("I: 1","I: 5","I:10","II: 1","II: 5","II:10","III: 1","III: 5","III:10")
  speed <- proc.time() - now

  list(out = out, speed  = speed) 
}

simul.KL <- function(Nsize = 500, nobs = 20, p = 3, k = 2, d = 5, alpha = rep(1, p + 2))
{ ## function to perform the simulation experiment (Section 4 of the manuscript)
  res <- matrix(0, nrow = Nsize, ncol = 13) # results container
  ok  <- matrix(0, nrow = Nsize, ncol = 6)
  now <- proc.time()

  cf <- rep(1, p +2)

  pb <- txtProgressBar(min = 0, max = Nsize, style = 3)
  # Monte Carlo iterations
  for (i in 1:Nsize) {
    z <- model.sim(n = nobs, p, k, d, a = c(1,1,0), b = c(0,0,1), alpha = alpha, cf = cf)
    res[i,1] <- scaled.condition(z$x)
    fm <- ridge(z$y ~ -1 + z$x)
    res[i,2] <- fm$lambda
    res[i,3] <- fm$edf
    lev <- leverages(fm)
    largest <- order(lev)[nobs]
    cutoff <- 2 * fm$edf / nobs
    ok[i,1] <- lev[largest] > cutoff
    res[i,4] <- largest
    res[i,5] <- cutoff
    CD <- cooks.distance(fm)
    largest <- order(CD)[nobs]
    res[i,6] <- largest
    res[i,7] <- CD[largest]
    KL <- KL.divergence.ridge(fm)
    largest <- order(KL)[nobs]
    res[i,8] <- largest
    res[i,9] <- KL[largest]
    o <- curvature.LD.ridge(fm, scheme = "scale")
    cutoff <- mean(o$hmax) + 2 * sd(o$hmax)
    res[i,10] <- o$hmax[nobs] > cutoff
    o <- curvature.LD.ridge(fm, scheme = "response")
    cutoff <- mean(o$hmax) + 2 * sd(o$hmax)
    res[i,11] <- o$hmax[nobs] > cutoff
    o <- curvature.KL.ridge(fm)
    cutoff <- mean(o$hmax) + 2 * sd(o$hmax)
    res[i,12] <- o$hmax[nobs] > cutoff
    u <- svd(z$x)$u
    u <- abs(u)
    cutoff <- apply(u, 2, mean) + 2 * apply(u, 2, sd)
    res[i,13] <- apply(u[,] > cutoff, 1, sum)[nobs] > 0
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)

  mnames <- c("condition","lambda","edf","leverage","cutoff","ID","cooks","ID","KL","LD.var","LD.res","KL.var","KL.res")
  colnames(res) <- mnames

  ok <- res[,c(6,8,10:13)]
  ok[,1] <- ifelse(ok[,1] == nobs, 1, 0)
  ok[,2] <- ifelse(ok[,2] == nobs, 1, 0)
  percentage <- apply(ok, 2, sum) / Nsize
  names(percentage) <- c("cooks","KL","LD.var","LD.res","KL.var","KL.res")

  speed <- proc.time() - now
  out <- list(results = res, ok = ok, percentage = 100 * percentage, speed = speed)
  out
}

model.sim <- function(n = 20, p = 3, k = 2, d = 5, a = c(1,1,0), b = c(0,0,1), alpha = rep(1, p + 2), cf = rep(1, p + 2))
{
  ret.n <- n
  m <- n - 1
  if (length(a) != length(b))
    stop("'a' and 'b' must have same length.")
  if (p != length(a))
    stop("length of vectors 'a' and 'b' must be equal to 'p'.")

  z <- matrix(rnorm(m * p), nrow = m)
  sigma2 <- 1. / 10^k
  e1 <- rnorm(m, sd = sqrt(sigma2))
  e2 <- rnorm(m, sd = sqrt(sigma2))

  w1 <- c(z %*% a) + e1
  w2 <- c(z %*% b) + e2

  x <- cbind(z, w1, w2)
  s <- crossprod(x)

  rs <- eigen(s, symm = TRUE)
  alpha <- alpha / minkowski(alpha)
  last <- d * c(rs$vectors %*% alpha)

  x <- rbind(x, last)
  rownames(x) <- NULL

  eps <- rnorm(n)
  y <- c(x %*% cf) + eps

  o <- list(x = x, y = y)
  o
}
