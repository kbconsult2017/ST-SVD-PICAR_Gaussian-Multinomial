# Code

R scripts for the ST SVD-PICAR paper. Each script is self-contained: it loads
packages, builds the spatio-temporal basis, sweeps the candidate basis grid
*K* ∈ {200, 300, 400, 500, 600}, selects the working rank by out-of-sample MSPE,
runs MCMC (NIMBLE), and writes figures/tables to a user-set `out_dir`.

## Contents

| Path | Setting | Method |
|---|---|---|
| `simulation/gaussian/sim_SVD-PICAR_Gaussian.R` | Gaussian (simulation) | ST SVD-PICAR |
| `simulation/gaussian/sim_PICAR_Gaussian.R` | Gaussian (simulation) | ST PICAR (no SVD) |
| `simulation/multinomial/sim_SVD_PICAR_mutltinomial.R` | Multinomial (simulation) | ST SVD-PICAR |
| `simulation/multinomial/sim_PICAR__multinomial__code.R` | Multinomial (simulation) | ST PICAR (no SVD) |
| `data-application/gaussian/data_svd-picar_Gaussian.R` | Gaussian (real data) | ST SVD-PICAR |
| `data-application/gaussian/data_picar_Gaussian.R` | Gaussian (real data) | ST PICAR (no SVD) |
| `data-application/multinomial/data_svd_picar_multinomial.R` | Multinomial (real data) | ST SVD-PICAR |
| `data-application/multinomial/data_picar_multinomial.R` | Multinomial (real data) | ST PICAR (no SVD) |

## Before running

1. **Simulation scripts** generate their own data; set `out_dir` (output folder) near the top.
2. **Data-application scripts** additionally need `data_path` set to the dataset in `../../data/real-data/`:
   - Gaussian → `NASA_GRACE-FO_DATA.xlsx`
   - Multinomial → `FINAL__1.CSV`

## R package dependencies

Extracted from the scripts:

`INLA`, `nimble`, `rsvd`, `irlba`, `Matrix`, `mvtnorm`, `VGAM`, `coda`, `sf`,
`fields`, `ggplot2`, `gridExtra`, `reshape2`, `dplyr`, `tidyr`, plus base
packages `grid` and `splines`.

```r
install.packages(c("nimble","rsvd","irlba","Matrix","mvtnorm","VGAM","coda",
                   "sf","fields","ggplot2","gridExtra","reshape2","dplyr","tidyr","readxl"))
# INLA is not on CRAN:
install.packages("INLA", repos = c(getOption("repos"),
                 INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```
