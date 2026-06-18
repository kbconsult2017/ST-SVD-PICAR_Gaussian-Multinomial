# Data

Inputs for the real-data applications. (The simulation scripts generate their
own data and do not read any file here.)

| File | Application | Description |
|---|---|---|
| `real-data/NASA_GRACE-FO_DATA.xlsx` | Gaussian | Terrestrial water storage / mass-change inputs derived from NASA's GRACE-FO mission; covariates include groundwater storage and root-zone soil moisture. |
| `real-data/FINAL__1.CSV` | Multinomial | County-month records (50,136 rows). Lung-cancer death counts discretized into four ordered categories; covariates include precipitation and wind. |

## Provenance / terms

- GRACE-FO products originate from NASA and are publicly available.
- The lung-cancer mortality response was derived via the George Mason University
  Air Quality lab; meteorological variables were derived from NASA NLDAS-2.

Please confirm and cite the original data sources and any redistribution terms
before sharing. If a dataset cannot be redistributed, replace the file here
with a small example or a download script and update `data_path` accordingly.
