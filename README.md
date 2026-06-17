# ST SVD-PICAR

Reproducibility materials for the paper

> **Compressed spatio-temporal basis functions: An SVD-PICAR framework for scalable Gaussian and multinomial inference**
> Kyei Baffour Afari and Yeongjin Gwon, Department of Biostatistics, University of Nebraska Medical Center.

ST SVD-PICAR integrates randomized singular value decomposition (rSVD) into spatio-temporal
basis construction (mesh-based Moran eigenvectors for space, B-splines for time) to select an
optimal basis size *K* in a data-driven way. This repository holds the code, data, and figures
needed to reproduce the simulation studies and the real-data applications under both Gaussian
and baseline-category multinomial likelihoods.

## Repository layout

```
ST-SVD-PICAR/
├── README.md
├── LICENSE
├── CITATION.cff
├── code/
│   ├── README.md
│   ├── simulation/
│   │   ├── gaussian/        # sim_SVD-PICAR_Gaussian.R, sim_PICAR_Gaussian.R
│   │   └── multinomial/     # sim_SVD_PICAR_mutltinomial.R, sim_PICAR__multinomial__code.R
│   └── data-application/
│       ├── gaussian/        # data_svd-picar_Gaussian.R, data_picar_Gaussian.R
│       └── multinomial/     # data_svd_picar_multinomial.R, data_picar_multinomial.R
├── data/
│   ├── README.md
│   └── real-data/           # NASA_GRACE-FO_DATA.xlsx, FINAL__1.CSV
└── figures/
    ├── README.md
    ├── main/                # Figures 1–8 (main paper)
    └── supplement/          # Figures A.1–A.8 (online appendix)
```

## Contents at a glance

- **8 R scripts** — SVD-PICAR and PICAR (no SVD), for Gaussian and multinomial,
  in both the simulation and real-data settings. See [`code/README.md`](code/README.md).
- **2 datasets** — GRACE-FO (Gaussian) and lung-cancer county-month (multinomial).
  See [`data/README.md`](data/README.md).
- **16 figures** — main Figures 1–8 and appendix Figures A.1–A.8, with a
  panel-by-panel index in [`figures/README.md`](figures/README.md).

All experiments use the candidate basis-size grid *K* ∈ {200, 300, 400, 500, 600}.

## Requirements

R (≥ 4.1) with:
`INLA`, `nimble`, `rsvd`, `irlba`, `Matrix`, `mvtnorm`, `VGAM`, `coda`, `sf`,
`fields`, `ggplot2`, `gridExtra`, `reshape2`, `dplyr`, `tidyr`, `readxl`
(and base `grid`, `splines`).

```r
install.packages(c("nimble","rsvd","irlba","Matrix","mvtnorm","VGAM","coda",
                   "sf","fields","ggplot2","gridExtra","reshape2","dplyr","tidyr","readxl"))
# INLA is not on CRAN:
install.packages("INLA", repos = c(getOption("repos"),
                 INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

## How to run

1. Open a script in `code/`.
2. Set `out_dir` (output folder) near the top.
3. For data-application scripts, also set `data_path` to the matching file in
   `data/real-data/` (Gaussian → `NASA_GRACE-FO_DATA.xlsx`; multinomial → `FINAL__1.CSV`).
4. Run the script in R; figures and tables are written to `out_dir`.

## Citation

See [`CITATION.cff`](CITATION.cff).

Afari, K. B., and Gwon, Y. Compressed spatio-temporal basis functions: An SVD-PICAR framework
for scalable Gaussian and multinomial inference.

## License

Code is released under the MIT License ([`LICENSE`](LICENSE)). Datasets retain the terms of
their original sources — see [`data/README.md`](data/README.md).

## Contact

Kyei Baffour Afari — kafari@unmc.edu — https://kbconsult2017.github.io
