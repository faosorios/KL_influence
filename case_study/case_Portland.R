## ID: case_Portland.R, last updated 2023-12-21, F.Osorio

## loading dataset and reading R sources
library(india) # load required package 'fastmatrix'
data(portland)
source("../code/cov_ridge.R")
source("../code/KL_influence.R")
source("../code/LD_influence.R")
source("../code/ridge_par.R")

## fitted model
fm <- ridge(y ~ ., data = portland, method = "GCV", x = TRUE)

## removing observations 3, 8 and 10
nobs <- fm$dims[1]
rm03 <- rm08 <- rm10 <- rep(TRUE, nobs)
rm03[3] <- rm08[8] <- rm10[10] <- FALSE
f03 <- ridge(y ~ ., data = portland, subset = rm03, method = "GCV")
f08 <- ridge(y ~ ., data = portland, subset = rm08, method = "GCV")
f10 <- ridge(y ~ ., data = portland, subset = rm10, method = "GCV")

## computing standard errors
SE00 <- vcov.ridge(fm)
SE03 <- vcov.ridge(f03)
SE08 <- vcov.ridge(f08)
SE10 <- vcov.ridge(f10)

## Table 1 (slightly recrafted)
tab1 <- matrix(0, nrow = 4, ncol = 14)
tab1[1,1:5] <- fm$coef
tab1[2,1:5] <- f03$coef
tab1[3,1:5] <- f08$coef
tab1[4,1:5] <- f10$coef
tab1[1,6:10] <- SE00
tab1[2,6:10] <- SE03
tab1[3,6:10] <- SE08
tab1[4,6:10] <- SE10
tab1[1,11] <- fm$scale
tab1[2,11] <- f03$scale
tab1[3,11] <- f08$scale
tab1[4,11] <- f10$scale
tab1[1,12] <- fm$lambda
tab1[2,12] <- f03$lambda
tab1[3,12] <- f08$lambda
tab1[4,12] <- f10$lambda
tab1[1,13] <- fm$edf
tab1[2,13] <- f03$edf
tab1[3,13] <- f08$edf
tab1[4,13] <- f10$edf
tab1[1,14] <- attr(SE00, "determinant") * 10^14
tab1[2,14] <- attr(SE03, "determinant") * 10^14
tab1[3,14] <- attr(SE08, "determinant") * 10^14
tab1[4,14] <- attr(SE10, "determinant") * 10^14

rownames(tab1) <- c("all", "3", "8", "10")
colnames(tab1) <- c("Int", "x1", "x2", "x3", "x4", "SE0", "SE1", "SE2", "SE3", "SE4", "sigma2", "lambda", "edf", "det")

tab1 # output below (SE: means standard error)
#            Int       x1       x2        x3        x4        SE0       SE1
# all 0.08544556 2.165338 1.158642 0.7383429 0.4895036 0.03997627 0.1699393
# 3   0.08383785 2.163506 1.159707 0.7357421 0.4898202 0.02640894 0.1789039
# 8   0.08397970 2.246701 1.119455 0.9030965 0.4821453 0.05072011 0.1443982
# 10  0.08884445 2.134533 1.163714 0.7259615 0.4922653 0.04226731 0.2794168
#            SE2      SE3      SE4    sigma2    lambda        edf       det
# all 0.04406509 0.146486 0.038286 5.0902402 1.9715705 3.97945569  4.503596
# 3   0.04651311 0.153995 0.040270 5.6343420 2.1882597 3.97673365  3.332206
# 8   0.03975024 0.138487 0.031847 3.4771145 1.2547484 3.98536092  2.586755
# 10  0.05838536 0.175789 0.044434 5.4456384 1.8673023 3.95999620 19.630077

## Table 1 (percentage changes)
chg <- matrix(0, nrow = 3, ncol = 9)
chg[1,] <- 100 * (tab1[2,-(6:10)] - tab1[1,-(6:10)]) / tab1[1,-(6:10)]
chg[2,] <- 100 * (tab1[3,-(6:10)] - tab1[1,-(6:10)]) / tab1[1,-(6:10)]
chg[3,] <- 100 * (tab1[4,-(6:10)] - tab1[1,-(6:10)]) / tab1[1,-(6:10)]

rownames(chg) <- c("3", "8", "10")
colnames(chg) <- c("Int", "x1", "x2", "x3", "x4", "sigma2", "lambda", "edf", "det")

## extracting the model matrix
x <- fm$x

## 'classical' influence measures
CD <- cooks.distance(fm, type = "cov")
LD <- logLik.displacement(fm, pars = "coef")
lev <- leverages(fm)
rel <- relative.condition(x)

rel # relative condition (reported in Section 4.1)
#           1            2            3            4            5            6 
# 0.036498473 -0.001126192  0.376952965 -0.002284163 -0.029407540 -0.039380700 
#           7            8            9           10           11           12 
#-0.035649700 -0.012596178  0.029996258  0.001835695  0.008286249  0.001986670 
#          13 
# 0.003488182 
#attr(,"scaled condition")
#[1] 249.5783

## selected value for the shrinkage parameter
lambda <- ridge.par(fm)

## local influence based on penalized-likelihood displacement
z <- curvature.LD.ridge(fm, scheme = "scale")

## Figure 1(a)
par(pty = "s")
plot(CD, ylab = "Cook's distances", ylim = c(0,.44), lwd = 2, cex.lab = 1.3)
text(8, CD[8], label = as.character(8), pos = 3)

## Figure 1(b)
par(pty = "s")
plot(LD, ylab = "Penalized likelihood displacement", ylim = c(0,2.5), lwd = 2, cex.lab = 1.3)
text(8, LD[8], label = as.character(8), pos = 3)

## Figure 1(c)
cutoff <- 2 * mean(lev)
par(pty = "s")
plot(lev, ylab = "Leverages", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(10, lev[10], label = as.character(10), pos = 3)

## Figure 1(d)
par(pty = "s")
plot(z$hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(8, z$hmax[8], label = as.character(8), pos = 3)

## entropy-based influence measures for ridge regression
KL <- KL.divergence.ridge(fm)
z <- curvature.KL.ridge(fm)

## Figure 1(e)
par(pty = "s")
plot(KL, ylab = "KL divergence", ylim = c(0,5), lwd = 2, cex.lab = 1.3)
text(10, KL[10], label = as.character(10), pos = 3)

## Figure 1(f)
par(pty = "s")
plot(z$hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(3, z$hmax[3], label = as.character(3), pos = 3)
text(10, z$hmax[10], label = as.character(10), pos = 1)
