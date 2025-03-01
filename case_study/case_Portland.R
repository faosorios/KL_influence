## ID: case_Portland.R, last updated 2025-03-01, F.Osorio

## loading dataset and reading R sources
library(india) # load required package 'fastmatrix'
data(portland)
source("../code/cov_ridge.R")
source("../code/KL_influence.R")
source("../code/LD_influence.R")
source("../code/GCV_influence.R")
source("../code/ridge_par.R")

## fitted model
fm <- ridge(y ~ ., data = portland, method = "GCV", x = TRUE)

## removing observations 3, 8 and 10
nobs <- fm$dims[1]
rm03 <- rm06 <- rm08 <- rm10 <- rep(TRUE, nobs)
rm03[3] <- rm06[6] <- rm08[8] <- rm10[10] <- FALSE
f03 <- ridge(y ~ ., data = portland, subset = rm03, method = "GCV")
f06 <- ridge(y ~ ., data = portland, subset = rm06, method = "GCV")
f08 <- ridge(y ~ ., data = portland, subset = rm08, method = "GCV")
f10 <- ridge(y ~ ., data = portland, subset = rm10, method = "GCV")

## computing standard errors
SE00 <- vcov.ridge(fm)
SE03 <- vcov.ridge(f03)
SE06 <- vcov.ridge(f06)
SE08 <- vcov.ridge(f08)
SE10 <- vcov.ridge(f10)

## Table 1 (slightly recrafted)
tab4 <- matrix(0, nrow = 5, ncol = 14)
tab4[1,1:5] <- fm$coef
tab4[2,1:5] <- f03$coef
tab4[3,1:5] <- f06$coef
tab4[4,1:5] <- f08$coef
tab4[5,1:5] <- f10$coef
tab4[1,6:10] <- SE00
tab4[2,6:10] <- SE03
tab4[3,6:10] <- SE06
tab4[4,6:10] <- SE08
tab4[5,6:10] <- SE10
tab4[1,11] <- fm$scale
tab4[2,11] <- f03$scale
tab4[3,11] <- f06$scale
tab4[4,11] <- f08$scale
tab4[5,11] <- f10$scale
tab4[1,12] <- fm$lambda
tab4[2,12] <- f03$lambda
tab4[3,12] <- f06$lambda
tab4[4,12] <- f08$lambda
tab4[5,12] <- f10$lambda
tab4[1,13] <- fm$edf
tab4[2,13] <- f03$edf
tab4[3,13] <- f06$edf
tab4[4,13] <- f08$edf
tab4[5,13] <- f10$edf
tab4[1,14] <- attr(SE00, "determinant") * 10^14
tab4[2,14] <- attr(SE03, "determinant") * 10^14
tab4[3,14] <- attr(SE06, "determinant") * 10^14
tab4[4,14] <- attr(SE08, "determinant") * 10^14
tab4[5,14] <- attr(SE10, "determinant") * 10^14

rownames(tab4) <- c("all", "3", "6", "8", "10")
colnames(tab4) <- c("Int", "x1", "x2", "x3", "x4", "SE0", "SE1", "SE2", "SE3", "SE4", "sigma2", "lambda", "edf", "det")

tab4 # output below (SE: means standard error)
#            Int       x1       x2       x3       x4      SE0       SE1
# all 0.08544558 2.165338 1.158642 0.738343 0.489504 0.039976 0.1699393
# 3   0.08383780 2.163506 1.159707 0.735742 0.489820 0.026409 0.1789039
# 6   0.08876270 2.132513 1.153963 0.737848 0.493333 0.045248 0.1445421
# 8   0.08397974 2.246701 1.119455 0.903097 0.482145 0.050720 0.1443982
# 10  0.08884449 2.134533 1.163714 0.725962 0.492265 0.042267 0.2794168
#            SE2      SE3      SE4   sigma2   lambda      edf       det
# all 0.04406509 0.146486 0.038286 5.090239 1.971569 3.979456  4.503600
# 3   0.04651311 0.153995 0.040270 5.634344 2.188263 3.976734  3.332202
# 6   0.03723397 0.123718 0.032347 3.599968 1.456899 3.985001  1.653894
# 8   0.03975024 0.138487 0.031847 3.477114 1.254747 3.985361  2.586757
# 10  0.05838536 0.175789 0.044434 5.445638 1.867301 3.959996 19.630096

## Table 1 (percentage changes)
chg <- matrix(0, nrow = 4, ncol = 9)
chg[1,] <- 100 * (tab4[2,-(6:10)] - tab4[1,-(6:10)]) / tab4[1,-(6:10)]
chg[2,] <- 100 * (tab4[3,-(6:10)] - tab4[1,-(6:10)]) / tab4[1,-(6:10)]
chg[3,] <- 100 * (tab4[4,-(6:10)] - tab4[1,-(6:10)]) / tab4[1,-(6:10)]
chg[4,] <- 100 * (tab4[5,-(6:10)] - tab4[1,-(6:10)]) / tab4[1,-(6:10)]

rownames(chg) <- c("3", "6", "8", "10")
colnames(chg) <- c("Int", "x1", "x2", "x3", "x4", "sigma2", "lambda", "edf", "det")

chg # output
#          Int          x1        x2          x3          x4     sigma2
# 3  -1.881645 -0.08460559   0.09187 -0.35225381  0.06468667  10.689165
# 6   3.882143 -1.51595793  -0.40387 -0.06706213  0.78232190 -29.277041
# 8  -1.715531  3.75750277  -3.38218 22.31397079 -1.50321824 -31.690567
# 10  3.977857 -1.42264465   0.43776 -1.67691483  0.56418080   6.981952
#       lambda         edf       det
# 3   10.99092 -0.06840341 -26.01025
# 6  -26.10460  0.13935290 -63.27617
# 8  -36.35794  0.14839300 -42.56245
# 10  -5.28860 -0.48899823 335.87566

## extracting the model matrix
x <- fm$x

## 'classical' influence measures
CD <- cooks.distance(fm, type = "cov")
LD <- logLik.displacement(fm, pars = "coef")
lev <- leverages(fm)
rel <- relative.condition(x)

## local influence based on penalized-likelihood displacement
z1 <- curvature.LD.ridge(fm, scheme = "scale")
z2 <- curvature.LD.ridge(fm, scheme = "response")

## Figure 1(a)
par(pty = "s")
plot(CD, ylab = "Cook's distances", ylim = c(0,.44), lwd = 2, cex.lab = 1.3)
text(8, CD[8], label = as.character(8), pos = 3)

## Figure 1(b)
par(pty = "s")
plot(LD, ylab = "Penalized likelihood displacement", ylim = c(0,2.5), lwd = 2, cex.lab = 1.3)
text(8, LD[8], label = as.character(8), pos = 3)

## Figure 1(c)
cutoff <- 2 * mean(lev) # 2 * edf / n
par(pty = "s")
plot(lev, ylab = "Leverages", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(10, lev[10], label = as.character(10), pos = 3)

## Figure 1(d)
par(pty = "s")
plot(rel, ylab = "Relative condition number", ylim = c(-0.04,0.4), lwd = 2, cex.lab = 1.3)
abline(h = 0, lwd = 2, col = "gray55")
text(3, rel[3], label = as.character(3), pos = 3)

## Figure 1(e)
hmax <- z1$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(8, hmax[8], label = as.character(8), pos = 3)

## Figure 1(f)
hmax <- z2$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(6, hmax[6], label = as.character(6), pos = 3)

## selected value for the shrinkage parameter
lambda <- ridge.par(fm)

## Figure 2(a)
obs <- c(6,8)
par(pty = "s")
plot(lambda, ylab = "GCV selection of lambda", ylim = c(1.25,2.41), lwd = 2, cex.lab = 1.3)
abline(h = fm$lambda, lwd = 2, lty = 2, col = "red")
text(obs, lambda[obs], label = as.character(obs), pos = 3)

## local influence based on generalized cross-validation criterion
z1 <- GCV.influence(fm, scheme = "scale")
z2 <- GCV.influence(fm, scheme = "response")

## Figure 2(b)
hmax <- z1$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(6, hmax[6], label = as.character(6), pos = 3)

## Figure 2(c)
hmax <- z2$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(8, hmax[8], label = as.character(8), pos = 3)

## case deletion Kullback-Leibler influence measure
KL <- KL.divergence.ridge(fm)

## Figure 3
par(pty = "s")
plot(KL, ylab = "KL divergence", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
text(10, KL[10], label = as.character(10), pos = 3)

## local influence based on the Kullback-Leibler divergence
z <- curvature.KL.ridge(fm)
u <- svd(x)$u

## Figure 4(a)
hmax <- z$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(10, hmax[10], label = as.character(10), pos = 1)

## Figure 4(b)
u <- abs(u) # magnitude of the columns of matrix 'U'
cutoffs <- apply(u, 2, mean) + 2 * apply(u, 2, sd)
ok <- u[,] > cutoffs
matplot(u, ylab = "eigenvectors", ylim = c(0,1), pch = c(0:4), col = 1, lwd = 2, cex.lab = 1.3)
abline(h = cutoffs[4], lty = 2, col = "red", lwd = 2)
abline(h = cutoffs[5], lty = 3, col = "red", lwd = 3)
text(3, u[3,5], label = as.character(3), pos = 3)
text(10, u[10,4], label = as.character(10), pos = 3)
