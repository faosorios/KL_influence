## ID: case_biomass.R, last updated 2024-03-27, F.Osorio

## loading dataset
biomass <- read.csv("../data/biomass.csv") 
## or loading a RDA file 
#biomass <- load("../data/biomass.rda")

## reading R sources
library(india) # load required package 'fastmatrix'
source("../code/cov_ridge.R")
source("../code/KL_influence.R")
source("../code/LD_influence.R")

## fitted model
biomass <- read.csv("biomass.csv")
fm <- ridge(y ~ ., data = biomass, x = TRUE)
x <- fm$x

## computing scaled condition
scaled.condition(x)

## 'classical' influence measures
levs <- leverages(fm)
CD <- cooks.distance(fm, type = "cov")
LD <- logLik.displacement(fm, pars = "coef")
z <- curvature.LD.ridge(fm, scheme = "scale")
hmax <- z$hmax

## KL-based influence diagnostics
KL <- KL.divergence.ridge(fm)
z <- curvature.KL.ridge(fm)

## Fig. 2.a
obs <- c(12,14,29,34)
pdf()
par(pty = "s")
plot(CD, ylab = "Cook's distances", ylim = c(0,.17), lwd = 2, cex.lab = 1.3)
text(obs, CD[obs], label = as.character(obs), pos = 3)
dev.off()

## Fig. 2.b
pdf()
par(pty = "s")
plot(LD, ylab = "Penalized likelihood displacement", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(obs, LD[obs], label = as.character(obs), pos = 3)
dev.off()

## Fig. 2.c
cutoff <- 2 * mean(levs)
pdf()
par(pty = "s")
plot(levs, ylab = "Leverages", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(5, lev[5], label = as.character(5), pos = 3)
dev.off()

## Fig. 2.d
obs <- c(33,34)
pdf()
par(pty = "s")
plot(hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(obs, hmax[obs], label = as.character(obs), pos = 3)
dev.off()

## Fig. 2.e
pdf()
par(pty = "s")
plot(KL, ylab = "KL divergence", ylim = c(0,4), lwd = 2, cex.lab = 1.3)
text(5, KL[5], label = as.character(5), pos = 3)
dev.off()

## Fig. 2.f
obs <- 29:30
pdf()
par(pty = "s")
plot(z$hmax, ylab = "hmax", ylim = c(0,1), lwd = 2, cex.lab = 1.3)
text(obs, z$hmax[obs], label = as.character(obs), pos = 3)
dev.off()

## removing individual observations 5, 12, 14, 29, 30, 33 and 34
nobs <- nrow(x)
rm05 <- rm12 <- rm14 <- rm29 <- rm30 <- rm33 <- rm34 <- rep(TRUE, nobs)
rm05[5] <- rm12[12] <- rm14[14] <- rm29[29] <- rm30[30] <- rm33[33] <- rm34[34] <- FALSE

f05 <- ridge(y ~ ., data = biomass, subset = rm05, method = "GCV")
f12 <- ridge(y ~ ., data = biomass, subset = rm12, method = "GCV")
f14 <- ridge(y ~ ., data = biomass, subset = rm14, method = "GCV")
f29 <- ridge(y ~ ., data = biomass, subset = rm29, method = "GCV")
f30 <- ridge(y ~ ., data = biomass, subset = rm30, method = "GCV")
f33 <- ridge(y ~ ., data = biomass, subset = rm33, method = "GCV")
f34 <- ridge(y ~ ., data = biomass, subset = rm34, method = "GCV")

SE00 <- vcov.ridge(fm)
SE05 <- vcov.ridge(f05)
SE12 <- vcov.ridge(f12)
SE14 <- vcov.ridge(f14)
SE29 <- vcov.ridge(f29)
SE30 <- vcov.ridge(f30)
SE33 <- vcov.ridge(f33)
SE34 <- vcov.ridge(f34)

## Table 2 (slightly recrafted)
tab2 <- matrix(0, nrow = 8, ncol = 16)
tab2[1,1:6] <- fm$coef
tab2[2,1:6] <- f05$coef
tab2[3,1:6] <- f12$coef
tab2[4,1:6] <- f14$coef
tab2[5,1:6] <- f29$coef
tab2[6,1:6] <- f30$coef
tab2[7,1:6] <- f33$coef
tab2[8,1:6] <- f34$coef
tab2[1,7:12] <- SE00
tab2[2,7:12] <- SE05
tab2[3,7:12] <- SE12
tab2[4,7:12] <- SE14
tab2[5,7:12] <- SE29
tab2[6,7:12] <- SE30
tab2[7,7:12] <- SE33
tab2[8,7:12] <- SE34
tab2[1,13] <- fm$scale
tab2[2,13] <- f05$scale
tab2[3,13] <- f12$scale
tab2[4,13] <- f14$scale
tab2[5,13] <- f29$scale
tab2[6,13] <- f30$scale
tab2[7,13] <- f33$scale
tab2[8,13] <- f34$scale
tab2[1,14] <- fm$lambda
tab2[2,14] <- f05$lambda
tab2[3,14] <- f12$lambda
tab2[4,14] <- f14$lambda
tab2[5,14] <- f29$lambda
tab2[6,14] <- f30$lambda
tab2[7,14] <- f33$lambda
tab2[8,14] <- f34$lambda
tab2[1,15] <- fm$edf
tab2[2,15] <- f05$edf
tab2[3,15] <- f12$edf
tab2[4,15] <- f14$edf
tab2[5,15] <- f29$edf
tab2[6,15] <- f30$edf
tab2[7,15] <- f33$edf
tab2[8,15] <- f34$edf
tab2[1,16] <- attr(SE00, "determinant")
tab2[2,16] <- attr(SE05, "determinant")
tab2[3,16] <- attr(SE12, "determinant")
tab2[4,16] <- attr(SE14, "determinant")
tab2[5,16] <- attr(SE29, "determinant")
tab2[6,16] <- attr(SE30, "determinant")
tab2[7,16] <- attr(SE33, "determinant")
tab2[8,16] <- attr(SE34, "determinant")

rownames(tab2) <- c("all", "5", "12", "14", "29", "30", "33", "34")
colnames(tab2) <- c("Int", "x1", "x2", "x3", "x4", "x5", "SE0", "SE1", "SE2", "SE3", "SE4", "SE5", "sigma2", "lambda", "edf", "det")

tab2 # output below (SE: standard error)
#          Int         x1       x2         x3           x4         x5  
#all 121.90659  -8.799582 364.8848 -0.1774392 -0.014186431  -9.096526
#5    91.52244  -9.138297 372.4452 -0.3285094 -0.007133257  -8.241800
#12  528.83381 -23.618894 368.3152 -0.1141567 -0.014155193 -11.604447
#14   99.27776 -14.803613 384.5190 -0.1729750 -0.012223916  -5.841072
#29  468.51604 -16.558863 378.5370 -0.2802240 -0.012841168 -14.370798
#30  125.93099  -8.957841 366.7685 -0.1828996 -0.014131983  -9.269722
#33  293.64801  -7.521296 337.6847 -0.2055611 -0.014130768 -13.729926
#34   49.17317  -1.985675 350.1561 -0.1627616 -0.019399493 -10.066801
#          SE0        SE1      SE2        SE3         SE4       SE5
#all  92.00414   9.147152 44.69188  0.3146619  0.01404216  6.820932
#5    70.03870   9.565807 48.19625  0.4852119  0.02212836  7.027682
#12  383.27492  11.281467 48.69116  0.2995123  0.01328196  7.565089
#14   76.65146   9.144028 43.33393  0.3009836  0.01346354  6.637481
#29  346.40227  10.584212 48.92690  0.3003570  0.01337077  7.449089
#30   95.23977   9.262139 47.03015  0.3196542  0.01419325  6.967548
#33  226.97408   9.682270 45.29159  0.2972436  0.01327628  6.938311
#34   49.22825   8.233237 39.90737  0.2841534  0.01275698  6.126590
#       sigma2     lambda      edf        det
#all  143560.8  1.2199199 5.060138  14381.509
#5    147465.7  1.6005198 5.029286  24800.036
#12   125431.3  0.1930693 5.345292 151019.675
#14   131424.4  1.4141168 5.046202   7202.690
#29   128112.4  0.2281310 5.307858 145181.796
#30   146633.7  1.1847782 5.060835  20652.098
#33   127385.2  0.4015173 5.198760  54269.847
#34   117334.7  2.1000826 5.014303   1563.688

## Table 2 (percentage changes)
all <- tab2[1,-(7:12)]
chg <- tab2[-1,-(7:12)]
chg[1,] <- 100 * (chg[1,] - all) / all
chg[2,] <- 100 * (chg[2,] - all) / all
chg[3,] <- 100 * (chg[3,] - all) / all
chg[4,] <- 100 * (chg[4,] - all) / all
chg[5,] <- 100 * (chg[5,] - all) / all
chg[6,] <- 100 * (chg[6,] - all) / all
chg[7,] <- 100 * (chg[7,] - all) / all

rownames(chg) <- c("5", "12", "14", "29", "30", "33", "34")
colnames(chg) <- c("Int", "x1", "x2", "x3", "x4", "x5", "sigma2", "lambda", "edf", "det")

chg
#         Int         x1         x2         x3          x4         x5     sigma2     lambda         edf       det
#5  -24.92412   3.849213  2.0719947  85.139150 -49.7177463  -9.396177   2.720003  31.198764 -0.60972106  72.44390
#12 333.80250 168.409277  0.9401202 -35.664327  -0.2201992  27.570100 -12.628447 -84.173607  5.63529032 950.09614
#14 -18.56243  68.230863  5.3809399  -2.515878 -13.8337485 -35.787883  -8.453848  15.918830 -0.27541002 -49.91701
#29 284.32381  88.177836  3.7415124  57.926790  -9.4827459  57.981174 -10.760900 -81.299506  4.89551623 909.50319
#30   3.30122   1.798484  0.5162340   3.077362  -0.3838075   1.903985   2.140499  -2.880651  0.01376437  43.60174
#33 140.87952 -14.526665 -7.4544445  15.848745  -0.3923720  50.935933 -11.267428 -67.086582  2.73947166 277.35849
#34 -59.66324 -77.434437 -4.0365332  -8.271891  36.7468144  10.666436 -18.268315  72.149228 -0.90582197 -89.12710

## joint effect of obs. 29 and 30 
rmtwo <- rep(TRUE, nobs)
rmtwo[29:30] <- FALSE
two <- ridge(y ~ ., data = biomass, subset = rmtwo, method = "GCV")
SEtwo <- vcov.ridge(two)
