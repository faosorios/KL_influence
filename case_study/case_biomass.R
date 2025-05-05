## ID: case_biomass.R, last updated 2025-03-02, F.Osorio

## loading dataset
biomass <- read.csv("../data/biomass.csv") 
## or loading a RDA file 
#biomass <- load("../data/biomass.rda")

## reading R sources
library(india) # load required package 'fastmatrix'
source("../code/cov_ridge.R")
source("../code/KL_influence.R")
source("../code/LD_influence.R")
source("../code/GCV_influence.R")
source("../code/ridge_par.R")

## fitted model
fm <- ridge(y ~ ., data = biomass, x = TRUE)
x <- fm$x

## computing scaled condition
scaled.condition(x)

## 'classical' influence measures
CD <- cooks.distance(fm, type = "cov")
LD <- logLik.displacement(fm, pars = "coef")
levs <- leverages(fm)

## Fig. 5(a)
obs <- c(12,14,29,34)
par(pty = "s")
plot(CD, ylab = "Cook's distances", ylim = c(0,.17), lwd = 2, cex.lab = 1.3)
text(obs, CD[obs], label = as.character(obs), pos = 3)

## Fig. 5(b)
pdf()
par(pty = "s")
plot(LD, ylab = "Penalized likelihood displacement", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(obs, LD[obs], label = as.character(obs), pos = 3)
dev.off()

## Fig. 5(c)
cutoff <- 2 * mean(levs)
par(pty = "s")
plot(levs, ylab = "Leverages", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(5, lev[5], label = as.character(5), pos = 3)

## local influence based on penalized-likelihood displacement
z1 <- curvature.LD.ridge(fm, scheme = "scale")
z2 <- curvature.LD.ridge(fm, scheme = "response")

## Figure 6(a)
hmax <- z1$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- c(33,34)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## Figure 6(b)
hmax <- z2$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- c(12,14,33,34)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## local influence based on generalized cross-validation criterion
z1 <- GCV.influence(fm, scheme = "scale")
z2 <- GCV.influence(fm, scheme = "response")

## Figure 6(c)
hmax <- z1$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- c(12,34)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## Figure 6(d)
hmax <- z2$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- c(7,11,15)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## case deletion Kullback-Leibler influence measure
KL <- KL.divergence.ridge(fm)

## Figure 7
par(pty = "s")
plot(KL, ylab = "KL divergence", ylim = c(0,4), lwd = 2, cex.lab = 1.3)
text(5, KL[5], label = as.character(5), pos = 3)

## local influence based on the Kullback-Leibler divergence
z <- curvature.KL.ridge(fm)
u <- svd(x)$u

## Fig. 8(a)
hmax <- z$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- 29:30
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(obs, z$hmax[obs], label = as.character(obs), pos = 3)

## Fig. 8(b)
u <- abs(u) # magnitude of the columns of matrix 'U'
cutoffs <- apply(u, 2, mean) + 2 * apply(u, 2, sd)
matplot(u, ylab = "eigenvectors", ylim = c(0,1), pch = c(0:4,8), col = 1, lwd = 2, cex.lab = 1.3)
abline(h = cutoffs[2], lty = 2, col = "red", lwd = 2)
abline(h = cutoffs[4], lty = 3, col = "red", lwd = 3)
text(5, u[5,2], label = as.character(5), pos = 3)
text(15, u[15,6], label = as.character(15), pos = 3)
text(c(28,29), u[28:29,4], label = as.character(28:29), pos = 3)
text(27, u[27,4], label = as.character(27), pos = 1)

## Fig. 9(a)
z2 <- explanatory.KL.ridge(fm, which = 2)
hmax <- z2$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(28:29, hmax[28:29], label = as.character(28:29), pos = 3)
text(27, hmax[27], label = as.character(27), pos = 1)

## Fig. 9(b)
z4 <- explanatory.KL.ridge(fm, which = 4)
hmax <- z4$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- (1:45)[hmax > cutoff]
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## Fig. 9(c)
z5 <- explanatory.KL.ridge(fm, which = 5)
hmax <- z5$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- (1:45)[hmax > cutoff]
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## Fig. 9(d)
z6 <- explanatory.KL.ridge(fm, which = 6)
hmax <- z6$hmax
cutoff <- mean(hmax) + 2 * sd(hmax)
obs <- (1:45)[hmax > cutoff]
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lty = 2, lwd = 2, col = "red")
text(obs, hmax[obs], label = as.character(obs), pos = 3)

## removing individual observations 5, 7, 11, 12, 14, 15, 29, 30, 33 and 34
nobs <- nrow(x)
rm05 <- rm07 <- rm11 <- rm12 <- rm14 <- rm15 <- rm29 <- rm30 <- rm33 <- rm34 <- rep(TRUE, nobs)
rm05[5] <- rm07[7] <- rm11[11] <- rm12[12] <- rm14[14] <- rm15[15] <- rm29[29] <- rm30[30] <- rm33[33] <- rm34[34] <- FALSE

f05 <- ridge(y ~ ., data = biomass, subset = rm05, method = "GCV")
f07 <- ridge(y ~ ., data = biomass, subset = rm07, method = "GCV")
f11 <- ridge(y ~ ., data = biomass, subset = rm11, method = "GCV")
f12 <- ridge(y ~ ., data = biomass, subset = rm12, method = "GCV")
f14 <- ridge(y ~ ., data = biomass, subset = rm14, method = "GCV")
f15 <- ridge(y ~ ., data = biomass, subset = rm15, method = "GCV")
f29 <- ridge(y ~ ., data = biomass, subset = rm29, method = "GCV")
f30 <- ridge(y ~ ., data = biomass, subset = rm30, method = "GCV")
f33 <- ridge(y ~ ., data = biomass, subset = rm33, method = "GCV")
f34 <- ridge(y ~ ., data = biomass, subset = rm34, method = "GCV")

SE00 <- vcov.ridge(fm)
SE05 <- vcov.ridge(f05)
SE07 <- vcov.ridge(f07)
SE11 <- vcov.ridge(f11)
SE12 <- vcov.ridge(f12)
SE14 <- vcov.ridge(f14)
SE15 <- vcov.ridge(f15)
SE29 <- vcov.ridge(f29)
SE30 <- vcov.ridge(f30)
SE33 <- vcov.ridge(f33)
SE34 <- vcov.ridge(f34)

## Table 5 and 6 (slightly recrafted)
tab5 <- matrix(0, nrow = 11, ncol = 16)
tab5[1,1:6] <- fm$coef
tab5[2,1:6] <- f05$coef
tab5[3,1:6] <- f07$coef
tab5[4,1:6] <- f11$coef
tab5[5,1:6] <- f12$coef
tab5[6,1:6] <- f14$coef
tab5[7,1:6] <- f15$coef
tab5[8,1:6] <- f29$coef
tab5[9,1:6] <- f30$coef
tab5[10,1:6] <- f33$coef
tab5[11,1:6] <- f34$coef
tab5[1,7:12] <- SE00
tab5[2,7:12] <- SE05
tab5[3,7:12] <- SE07
tab5[4,7:12] <- SE11
tab5[5,7:12] <- SE12
tab5[6,7:12] <- SE14
tab5[7,7:12] <- SE15
tab5[8,7:12] <- SE29
tab5[9,7:12] <- SE30
tab5[10,7:12] <- SE33
tab5[11,7:12] <- SE34
tab5[1,13] <- fm$scale
tab5[2,13] <- f05$scale
tab5[3,13] <- f07$scale
tab5[4,13] <- f11$scale
tab5[5,13] <- f12$scale
tab5[6,13] <- f14$scale
tab5[7,13] <- f15$scale
tab5[8,13] <- f29$scale
tab5[9,13] <- f30$scale
tab5[10,13] <- f33$scale
tab5[11,13] <- f34$scale
tab5[1,14] <- fm$lambda
tab5[2,14] <- f05$lambda
tab5[3,14] <- f07$lambda
tab5[4,14] <- f11$lambda
tab5[5,14] <- f12$lambda
tab5[6,14] <- f14$lambda
tab5[7,14] <- f15$lambda
tab5[8,14] <- f29$lambda
tab5[9,14] <- f30$lambda
tab5[10,14] <- f33$lambda
tab5[11,14] <- f34$lambda
tab5[1,15] <- fm$edf
tab5[2,15] <- f05$edf
tab5[3,15] <- f07$edf
tab5[4,15] <- f11$edf
tab5[5,15] <- f12$edf
tab5[6,15] <- f14$edf
tab5[7,15] <- f15$edf
tab5[8,15] <- f29$edf
tab5[9,15] <- f30$edf
tab5[10,15] <- f33$edf
tab5[11,15] <- f34$edf
tab5[1,16] <- attr(SE00, "determinant")
tab5[2,16] <- attr(SE05, "determinant")
tab5[3,16] <- attr(SE07, "determinant")
tab5[4,16] <- attr(SE11, "determinant")
tab5[5,16] <- attr(SE12, "determinant")
tab5[6,16] <- attr(SE14, "determinant")
tab5[7,16] <- attr(SE15, "determinant")
tab5[8,16] <- attr(SE29, "determinant")
tab5[9,16] <- attr(SE30, "determinant")
tab5[10,16] <- attr(SE33, "determinant")
tab5[11,16] <- attr(SE34, "determinant")

rownames(tab5) <- c("all", "5", "7", "11", "12", "14", "15", "29", "30", "33", "34")
colnames(tab5) <- c("Int", "x1", "x2", "x3", "x4", "x5", "SE0", "SE1", "SE2", "SE3", "SE4", "SE5", "sigma2", "lambda", "edf", "det")

tab5 # output below (SE: standard error)
#          Int         x1       x2         x3           x4         x5       SE0       SE1
#all 121.90659  -8.799582 364.8848 -0.1774392 -0.014186431  -9.096526  92.00414  9.147152
#5    91.52244  -9.138297 372.4452 -0.3285094 -0.007133257  -8.241800  70.03870  9.565807
#7    94.17728  -7.481725 360.9065 -0.1293103 -0.015701311  -9.302207  72.05527  9.338031
#11  337.61667 -12.920229 355.5875 -0.1975490 -0.013536498 -11.340874 261.61926 10.301489
#12  528.83381 -23.618894 368.3152 -0.1141567 -0.014155193 -11.604447 383.27492 11.281467
#14   99.27776 -14.803613 384.5190 -0.1729750 -0.012223916  -5.841072  76.65146  9.144028
#15   77.66601  -8.563058 366.7485 -0.1489074 -0.014576404  -8.702108  61.22241  9.285688
#29  468.51604 -16.558863 378.5370 -0.2802240 -0.012841168 -14.370798 346.40227 10.584212
#30  125.93099  -8.957841 366.7685 -0.1828996 -0.014131983  -9.269722  95.23977  9.262139
#33  293.64801  -7.521296 337.6847 -0.2055611 -0.014130768 -13.729926 226.97408  9.682270
#34   49.17317  -1.985675 350.1561 -0.1627616 -0.019399493 -10.066801  49.22825  8.233237
#         SE2       SE3        SE4      SE5   sigma2    lambda      edf        det
#all 44.69188 0.3146619 0.01404216 6.820932 143560.8 1.2199199 5.060138  14381.509
#5   48.19625 0.4852119 0.02212836 7.027682 147465.7 1.6005198 5.029286  24800.036
#7   46.43163 0.3416534 0.01471752 6.977788 147373.5 1.5063771 5.032542  12125.636
#11  47.79576 0.3167852 0.01428605 7.374957 143274.5 0.3464895 5.205237 127790.097
#12  48.69116 0.2995123 0.01328196 7.565089 125431.3 0.1930693 5.345292 151019.675
#14  43.33393 0.3009836 0.01346354 6.637481 131424.4 1.4141168 5.046202   7202.690
#15  44.87199 0.3202119 0.01422623 6.864399 147483.5 1.8001710 5.021405   7604.825
#29  48.92690 0.3003570 0.01337077 7.449089 128112.4 0.2281310 5.307858 145181.796
#30  47.03015 0.3196542 0.01419325 6.967548 146633.7 1.1847782 5.060835  20652.098
#33  45.29159 0.2972436 0.01327628 6.938311 127385.2 0.4015173 5.198760  54269.847
#34  39.90737 0.2841534 0.01275698 6.126590 117334.7 2.1000826 5.014303   1563.688

## joint effect of obs. 29 and 30 
rmtwo <- rep(TRUE, nobs)
rmtwo[29:30] <- FALSE
two <- ridge(y ~ ., data = biomass, subset = rmtwo, method = "GCV")
SEtwo <- vcov.ridge(two)
