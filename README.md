# Influence diagnostics for ridge regression

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![fastmatrix](https://img.shields.io/badge/fastmatrix-0.5--77-orange)](https://cran.r-project.org/package=fastmatrix)
[![india](https://img.shields.io/badge/india-0.1-orange)](https://cran.r-project.org/package=india)

Supplementary material to **Entropy-based influence diagnostics for ridge regression** by Alonso Ogueda and Felipe Osorio

Code tested on R under development (2022-01-27 r81578), running Linux Zorin 16 (64 bits)

Attached packages: fastmatrix 0.5-77, india 0.1

CONTENTS:
- case_study/case_Portland.R: R commands for the analysis of Portland cement dataset (described/analyzed at Section 4.1 from manuscript).
- case_study/case_biomass.R: R commands for the analysis of Aerial biomass dataset (described/analyzed at Section 4.2 from manuscript)
- code/KL_influence.R: R functions to compute diagnostic measures for ridge regression based on Kullback-Liebler divergence.
- code/LD_influence.R: R function to assess the local influence based on penalized likelihood displacement.
- code/cov_ridge.R: R function to obtain the covariance matrix (and its determinant) for the ridge estimator.
- code/ridge_par.R: R function to carry out the shrinkage parameter selection when the i-th observation has been deleted from the dataset.
- data/biomass.csv: Aerial biomass production of the marsh grass Spartina alterniflora dataset, in CSV format
- data/biomass.rda: Aerial biomass production of the marsh grass Spartina alterniflora dataset, in RDA format.
- README.md: this file.
