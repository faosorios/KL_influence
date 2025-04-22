# Influence diagnostics for ridge regression based on the Kullback-Leibler divergence

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![fastmatrix](https://img.shields.io/badge/fastmatrix-0.5--7721-orange)](https://cran.r-project.org/package=fastmatrix)
[![india](https://img.shields.io/badge/india-0.1-orange)](https://cran.r-project.org/package=india)

Supplementary material to **Influence diagnostics for ridge regression based on the Kullback-Leibler divergence** by Alonso Ogueda and Felipe Osorio (Statistical Papers, Accepted)

Code tested on R under development (2025-02-20 r87772), running Linux Mint 22.1 (64 bits)

Attached packages: fastmatrix 0.5-7721, india 0.1

CONTENTS:
- case_study/case_Portland.R: R commands for the analysis of Portland cement dataset (described/analyzed at Section 4.2.1 from manuscript).
- case_study/case_biomass.R: R commands for the analysis of Aerial biomass dataset (described/analyzed at Section 4.2.2 from manuscript)
- code/GCV_influence.R: R functions to compute diagnostic measures on the ridge (shrinkage) parameter based on the GCV criterion.
- code/KL_influence.R: R functions to compute diagnostic measures for ridge regression based on Kullback-Liebler divergence.
- code/LD_influence.R: R function to assess the local influence based on penalized likelihood displacement.
- code/cov_ridge.R: R function to obtain the covariance matrix (and its determinant) for the ridge estimator.
- code/ridge_par.R: R function to carry out the shrinkage parameter selection when the i-th observation has been deleted from the dataset.
- simulation/simul.R: R function to compute the outlier detection percentage using different influence measures.
- simulation/oiutput.R: R commands to perform the simulation study described at Section 4.1 from manuscript.
- data/biomass.csv: Aerial biomass production of the marsh grass Spartina alterniflora dataset, in CSV format
- data/biomass.rda: Aerial biomass production of the marsh grass Spartina alterniflora dataset, in RDA format.
- README.md: this file.
