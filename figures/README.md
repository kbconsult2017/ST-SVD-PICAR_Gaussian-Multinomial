# Figures

Candidate basis sizes: *K* ∈ {200, 300, 400, 500, 600}. "ESS/sec" is mean effective sample
size per second across all monitored parameters; "best rank" is the working rank *r* ≤ *K* that
minimizes out-of-sample MSPE.

## Main paper (`main/`)

| File | Paper | Description |
|---|---|---|
| `Fig1_Gaussian_SVDPICAR_runtime_vs_K.png` | Figure 1 | Gaussian, ST SVD-PICAR: best-rank MCMC sampling runtime (sec) vs *K* — runtime falls then stabilizes near *K* = 500–600. |
| `Fig2_Gaussian_SVDPICAR_ESSsec_vs_runtime.png` | Figure 2 | Gaussian, ST SVD-PICAR: mean ESS/sec vs best-rank runtime, points labeled by *K* — efficiency peaks at *K* = 500. |
| `Fig3_Gaussian_PICAR_runtime_vs_K.png` | Figure 3 | Gaussian, ST PICAR (no SVD): best-rank runtime vs *K* — peaks near *K* = 300, no stable floor. |
| `Fig4_Gaussian_PICAR_ESSsec_vs_runtime.png` | Figure 4 | Gaussian, ST PICAR (no SVD): mean ESS/sec vs runtime, labeled by *K* — efficiency highest at *K* = 200. |
| `Fig5_Multinomial_SVDPICAR_runtime_vs_K.png` | Figure 5 | Multinomial, ST SVD-PICAR: best-rank total MCMC runtime (NIMBLE) vs *K* — peak at *K* = 300, lowest at *K* = 400. |
| `Fig6_Multinomial_SVDPICAR_ESSsec_vs_runtime.png` | Figure 6 | Multinomial, ST SVD-PICAR: mean ESS/sec vs runtime, labeled by *K* — best efficiency near *K* = 400. |
| `Fig7_Multinomial_PICAR_runtime_vs_K.png` | Figure 7 | Multinomial, ST PICAR (no SVD): best-rank total runtime vs *K* — sharp increase at *K* = 600. |
| `Fig8_Multinomial_PICAR_ESSsec_vs_runtime.png` | Figure 8 | Multinomial, ST PICAR (no SVD): mean ESS/sec vs runtime, labeled by *K* — *K* = 600 collapses to lowest ESS/sec. |

## Online appendix (`supplement/`)

Latent-field maps share a common color scale within each figure; *t\** is the median training
time. The roughness ratio is √R(û)/sd(û); a value near the dashed line at 1 indicates the fitted
field matches the truth's local variation, while values well below 1 indicate oversmoothing.

| File | Paper | Description |
|---|---|---|
| `FigA1_Gaussian_SVDPICAR_latent_field.png` | Figure A.1 | Gaussian ST SVD-PICAR posterior-mean latent field û(s, t\*) per *K* (t\* = 2017.4167); clear local structure retained. |
| `FigA2_Gaussian_SVDPICAR_roughness_vs_rank.png` | Figure A.2 | Gaussian ST SVD-PICAR scale-free roughness ratio vs rank *r* per *K*; near 1 at the selected rank. |
| `FigA3_Gaussian_PICAR_latent_field.png` | Figure A.3 | Gaussian ST PICAR (no SVD) posterior-mean latent field; nearly flat (oversmoothed). |
| `FigA4_Gaussian_PICAR_roughness_vs_rank.png` | Figure A.4 | Gaussian ST PICAR (no SVD) roughness ratio vs rank *r*; stays well below 1. |
| `FigA5_Multinomial_SVDPICAR_latent_field.png` | Figure A.5 | Multinomial ST SVD-PICAR posterior-mean shared latent field (t\* = 2017.5833); visible local structure. |
| `FigA6_Multinomial_SVDPICAR_roughness_vs_rank.png` | Figure A.6 | Multinomial ST SVD-PICAR scale-free roughness (RMS/sd) vs rank *r* per *K*. |
| `FigA7_Multinomial_PICAR_latent_field.png` | Figure A.7 | Multinomial ST PICAR (no SVD) posterior-mean shared latent field; nearly uniform (oversmoothed). |
| `FigA8_Multinomial_PICAR_roughness_vs_rank.png` | Figure A.8 | Multinomial ST PICAR (no SVD) roughness (RMS/sd) vs rank *r* per *K*. |
