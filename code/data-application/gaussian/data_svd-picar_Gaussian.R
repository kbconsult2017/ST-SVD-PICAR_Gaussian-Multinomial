# Gaussian data application: ST SVD-PICAR

rm(list = ls())
options(stringsAsFactors = FALSE)

# 0) USER CONTROLS

CENTER_ST_BASIS <- FALSE

data_path <- "C:/Users/Admin/OneDrive - University of Nebraska Medical Center/Desktop/NASA GRACE-FO DATA.xlsx"

out_dir <- "C:/Users/Admin/OneDrive - University of Nebraska Medical Center/Desktop/d_svd_continuous"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

constructed_K_values <- c(200, 300, 400, 500, 600)

train_frac <- 0.70

beta0_true <- NA_real_
beta1_true <- NA_real_
beta2_true <- NA_real_

tau_true <- 2.0

num_eigen_full_req <- 30
num_t_basis_req    <- 20

SAFE_RANK_FRAC   <- 0.95
SAFE_RANK_BUFFER <- 20
pX <- 2

safe_max_rank_gaussian <- function(n_train, pX, frac = 0.95, buffer = 20) {
  rhs  <- floor(frac * (n_train - buffer))
  rmax <- rhs - (1 + pX)
  max(5, rmax)
}

USE_RANK_GRID <- TRUE
make_rank_grid <- function(Kmax) {
  if (!USE_RANK_GRID) return(5:Kmax)
  v1 <- 5:25
  v2 <- if (Kmax >= 30)  seq(30, min(200, Kmax), by=10) else numeric(0)
  v3 <- if (Kmax >= 225) seq(225, Kmax, by=25)          else numeric(0)
  v  <- sort(unique(c(v1, v2, v3, Kmax)))
  v[v >= 5 & v <= Kmax]
}

RUN_CV_BLOCK <- TRUE
RUN_CV_BEST_RANK_ONLY <- TRUE

RUN_MCMC_BLOCK <- TRUE
RUN_MCMC_ALL_RANKS <- FALSE

niter_mcmc_rank <- 60000
burn_prop_rank  <- 0.50
MONITOR_DELTA <- FALSE

# 1) PACKAGES
required_packages <- c(
  "readxl","fields","mvtnorm","INLA","ggplot2","sf","Matrix","irlba","splines",
  "reshape2","dplyr","tidyr","nimble","coda","rsvd","gridExtra","grid"
)

installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    if (pkg == "INLA") {
      install.packages("INLA",
                       repos = c(getOption("repos"),
                                 INLA = "https://inla.r-inla-download.org/R/stable"),
                       dep = TRUE)
    } else {
      install.packages(pkg)
    }
  }
}

library(readxl)
library(fields)
library(mvtnorm)
library(INLA)
library(ggplot2)
library(sf)
library(Matrix)
library(irlba)
library(splines)
library(reshape2)
library(dplyr)
library(tidyr)
library(nimble)
library(coda)
library(rsvd)
library(gridExtra)
library(grid)

# 2) HELPERS

pdf_text_page <- function(title, lines, cex_title=1.1, cex_body=0.9) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow=20, ncol=1)))
  grid.text(title, x=0.02, y=unit(19, "lines"), just="left",
            gp=gpar(fontsize=16*cex_title, fontface="bold"))
  yy <- 17
  for (ln in lines) {
    grid.text(paste0("• ", ln), x=0.03, y=unit(yy, "lines"), just="left",
              gp=gpar(fontsize=12*cex_body))
    yy <- yy - 1.2
    if (yy < 1) break
  }
  popViewport()
}

bold_table_cell2 <- function(tg, row, col) {
  nm <- paste0("core-", row, "-", col)
  idx <- which(tg$layout$name == nm)
  if (length(idx) == 1) tg$grobs[[idx]]$gp <- gpar(fontface="bold")
  tg
}

mesh_to_laplacian <- function(mesh) {
  Nnew <- nrow(mesh$loc)
  triangles <- mesh$graph$tv
  edges <- rbind(triangles[, c(1,2)], triangles[, c(1,3)], triangles[, c(2,3)])
  edges <- t(apply(edges, 1, function(x) sort(x)))
  edges <- unique(edges)

  i_idx <- c(edges[,1], edges[,2])
  j_idx <- c(edges[,2], edges[,1])
  W <- sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(Nnew, Nnew))
  D <- Diagonal(x = rowSums(W))
  Q <- D - W
  list(W = W, Q = Q, Nnew = Nnew)
}

build_moran_basis <- function(W_adj, num_eigen) {
  N <- nrow(W_adj)
  one <- Matrix(1, N, 1)
  OrthSpace <- Diagonal(N) - (1/N) * (one %*% t(one))
  MoransOperator <- OrthSpace %*% (W_adj %*% OrthSpace)
  eig <- irlba(MoransOperator, nv = num_eigen, nu = num_eigen)
  eig$v[, 1:num_eigen, drop = FALSE]
}

make_st_scores <- function(spatial_scores, temporal_scores) {
  n <- nrow(spatial_scores)
  p <- ncol(spatial_scores)
  q <- ncol(temporal_scores)
  ST <- matrix(NA_real_, nrow=n, ncol=p*q)
  for (i in 1:n) {
    ST[i, ] <- as.vector(tcrossprod(temporal_scores[i, ], spatial_scores[i, ]))
  }
  ST
}

lmfit_with_ci <- function(Z, y) {
  fit <- try(lm.fit(x = Z, y = y), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  if (any(!is.finite(fit$coefficients))) return(NULL)

  n <- nrow(Z)
  df <- n - fit$rank
  if (!is.finite(df) || df <= 0) return(NULL)

  rss <- sum(fit$residuals^2)
  sig2_hat <- rss / df

  R <- try(qr.R(fit$qr), silent = TRUE)
  if (inherits(R, "try-error")) return(NULL)

  XtX_inv <- try(chol2inv(R), silent = TRUE)
  if (inherits(XtX_inv, "try-error")) return(NULL)

  se <- sqrt(pmax(0, diag(XtX_inv) * sig2_hat))
  tcrit <- qt(0.975, df = df)

  beta <- fit$coefficients
  lcl  <- beta - tcrit * se
  ucl  <- beta + tcrit * se

  list(beta = beta, lcl = lcl, ucl = ucl, fit = fit)
}

# 3) NIMBLE PICAR MCMC

run_picar_mcmc_one_rank <- function(y_train, X_matrix, M_matrix, MQM_reduced,
                                    niter = 20000, burn_prop = 0.70,
                                    monitor_delta = FALSE) {

  n_train <- nrow(X_matrix)
  p_rank  <- ncol(M_matrix)

  data_list_nimble <- list(Z = as.numeric(y_train))

  consts_nimble <- list(
    n    = n_train,
    p    = p_rank,
    X    = as.matrix(X_matrix),
    M    = as.matrix(M_matrix),
    MQM  = as.matrix(MQM_reduced),
    zero = rep(0, p_rank)
  )

  inits_nimble <- list(beta0=0, beta1=0, beta2=0, delta=rep(0,p_rank), sigma2=1, tau=1)

  linear_PICAR_code <- nimbleCode({
    for (j in 1:n) {
      mu[j] <- beta0 + beta1*X[j,1] + beta2*X[j,2] + inprod(M[j,1:p], delta[1:p])
      Z[j] ~ dnorm(mean=mu[j], var=sigma2)
    }
    precMat[1:p,1:p] <- tau * MQM[1:p,1:p]
    delta[1:p] ~ dmnorm(mean=zero[1:p], prec=precMat[1:p,1:p])

    beta0 ~ dnorm(0, var=100)
    beta1 ~ dnorm(0, var=100)
    beta2 ~ dnorm(0, var=100)
    sigma2 ~ dinvgamma(shape=2, scale=4)
    tau ~ dgamma(shape=2, rate=1)
  })

  Rmodel <- nimbleModel(code=linear_PICAR_code, data=data_list_nimble,
                        constants=consts_nimble, inits=inits_nimble)

  nimbleOptions(MCMCprogressBar = FALSE)

  conf <- configureMCMC(Rmodel, print=FALSE)
  monitors <- c("beta0","beta1","beta2","sigma2","tau")
  if (isTRUE(monitor_delta)) monitors <- c(monitors, "delta")
  conf$setMonitors(monitors)

  mcmc  <- buildMCMC(conf)

  t_compile0 <- proc.time()
  try(nimble::clearAllCompiledNimbleFunctions(), silent = TRUE)
  try(nimble:::clearAllCompiledNimbleFunctions(), silent = TRUE)
  gc()

  Cmodel <- compileNimble(Rmodel)
  Cmcmc  <- compileNimble(mcmc, project = Cmodel)

  compile_time_sec <- as.numeric((proc.time() - t_compile0)[3])

  t_run0 <- proc.time()
  Cmcmc$run(niter)
  runtime_sec <- as.numeric((proc.time() - t_run0)[3])

  samples <- as.matrix(Cmcmc$mvSamples)

  burnin <- floor(burn_prop * nrow(samples))
  post   <- samples[(burnin+1):nrow(samples), , drop=FALSE]

  ess_vec <- coda::effectiveSize(coda::as.mcmc(post))
  ess_vec <- as.numeric(ess_vec)
  names(ess_vec) <- names(coda::effectiveSize(coda::as.mcmc(post)))

  essps_vec <- ess_vec / max(runtime_sec, 1e-12)

  ess_long <- data.frame(
    Parameter   = names(ess_vec),
    ESS         = as.numeric(ess_vec),
    ESS_per_sec = as.numeric(essps_vec),
    stringsAsFactors = FALSE
  )

  ess_long$Group <- ifelse(grepl("^delta\\[", ess_long$Parameter), "delta",
                           ifelse(grepl("^beta", ess_long$Parameter), "beta",
                                  ifelse(ess_long$Parameter == "sigma2", "sigma2",
                                         ifelse(ess_long$Parameter == "tau", "tau", "other"))))

  get_or_na <- function(nm) if (nm %in% names(ess_vec)) as.numeric(ess_vec[nm]) else NA_real_
  getps_or_na <- function(nm) if (nm %in% names(essps_vec)) as.numeric(essps_vec[nm]) else NA_real_

  essps_mean_beta12 <- mean(c(getps_or_na("beta1"), getps_or_na("beta2")), na.rm=TRUE)

  list(
    samples = samples,
    post    = post,

    compile_time_sec = compile_time_sec,
    runtime_sec = runtime_sec,

    beta0_ESS = get_or_na("beta0"),
    beta1_ESS = get_or_na("beta1"),
    beta2_ESS = get_or_na("beta2"),

    beta0_ESS_per_sec = getps_or_na("beta0"),
    beta1_ESS_per_sec = getps_or_na("beta1"),
    beta2_ESS_per_sec = getps_or_na("beta2"),

    ESS_per_sec_mean_beta   = essps_mean_beta12,
    ESS_per_sec_mean_beta12 = essps_mean_beta12,

    ess_long = ess_long
  )
}

# 4) READ REAL DATA + TRAIN/TEST SPLIT + BUILD MESH + BASIS INGREDIENTS

overall_start <- Sys.time()
set.seed(123)

lat_range <- c(32.5343, 42.0095)
lon_range <- c(-124.4096, -113.0460)

dat0 <- readxl::read_excel(data_path)
dat0 <- as.data.frame(dat0)
names(dat0) <- trimws(names(dat0))

need_cols <- c("PM25_MEAN", "gws_MEAN", "rtz_MEAN")
stopifnot(all(need_cols %in% names(dat0)))

time_col <- "year_mo"
if (!time_col %in% names(dat0)) {
  stop(paste0(
    'Expected time column "', time_col, '" was not found in the dataset.\n',
    "Column names in your file are:\n  ",
    paste(names(dat0), collapse = ", ")
  ))
}
has_time <- TRUE

id_col <- if ("GEOID" %in% names(dat0)) {
  "GEOID"
} else if ("County" %in% names(dat0)) {
  "County"
} else if ("Location" %in% names(dat0)) {
  "Location"
} else {
  stop("No GEOID/County/Location column found to assign synthetic coordinates.")
}

ids <- sort(unique(dat0[[id_col]]))
set.seed(123)
coord_tbl <- data.frame(
  id  = ids,
  lon = runif(length(ids), min = lon_range[1], max = lon_range[2]),
  lat = runif(length(ids), min = lat_range[1], max = lat_range[2]),
  stringsAsFactors = FALSE
)
names(coord_tbl)[1] <- id_col

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
dat0 <- dplyr::left_join(dat0, coord_tbl, by = id_col)

stopifnot(all(c("lon","lat") %in% names(dat0)))

keep_cols <- unique(c(need_cols, id_col, time_col, "lon", "lat"))
dat <- dat0[, keep_cols, drop = FALSE]
names(dat)[names(dat) == time_col] <- "year_mo"

dat <- dat[
  is.finite(dat$PM25_MEAN) & is.finite(dat$gws_MEAN) & is.finite(dat$rtz_MEAN) &
    is.finite(dat$lon) & is.finite(dat$lat) & !is.na(dat$year_mo),
]

dat <- dat[
  dat$lat >= lat_range[1] & dat$lat <= lat_range[2] &
    dat$lon >= lon_range[1] & dat$lon <= lon_range[2],
]

if (nrow(dat) < 50) {
  stop(paste0("After cleaning, only n=", nrow(dat), " rows remain."))
}

if (inherits(dat$year_mo, c("Date","POSIXct","POSIXt"))) {
  dat$year_mo <- as.numeric(dat$year_mo) / 365.25
} else {
  dat$year_mo <- as.numeric(dat$year_mo)
}

n_total <- nrow(dat)
n_label <- paste0("n", n_total)

scale2 <- function(x) as.numeric(scale(x))
dat$gws_MEAN <- scale2(dat$gws_MEAN)
dat$rtz_MEAN <- scale2(dat$rtz_MEAN)

set.seed(123)
train_idx <- sample(seq_len(nrow(dat)), size = floor(train_frac * nrow(dat)))
train_data <- dat[train_idx, ]
test_data  <- dat[-train_idx, ]

response_train   <- train_data$PM25_MEAN
response_test    <- test_data$PM25_MEAN
covariates_train <- train_data[, c("gws_MEAN","rtz_MEAN")]
covariates_test  <- test_data[,  c("gws_MEAN","rtz_MEAN")]

train_coords <- as.matrix(train_data[, c("lon","lat")])
max.edge0 <- 1.5

mesh <- INLA::inla.mesh.2d(
  loc      = train_coords,
  max.edge = c(1,2) * max.edge0,
  cutoff   = max.edge0/5,
  offset   = c(max.edge0, 6.0)
)

pdf(file.path(out_dir, paste0("01_mesh_and_locations_", n_label, "_", timestamp, ".pdf")), width=8, height=6)
plot(mesh, main="INLA Mesh + Training Locations (synthetic coords)")
points(train_coords, col="blue", pch=16, cex=0.5)
dev.off()

AMat_train <- INLA::inla.spde.make.A(mesh, loc = as.matrix(train_data[, c("lon","lat")]))
AMat_test  <- INLA::inla.spde.make.A(mesh, loc = as.matrix(test_data[,  c("lon","lat")]))

lap <- mesh_to_laplacian(mesh)
W_adj <- lap$W
Q_lap <- as.matrix(lap$Q) + diag(1e-6, lap$Nnew)

num_eigen_full <- min(num_eigen_full_req, max(5L, lap$Nnew - 2L))
mBase <- build_moran_basis(W_adj, num_eigen = num_eigen_full)

n_time_unique <- length(unique(train_data$year_mo))
if (n_time_unique <= 1L) {
  num_t_basis <- 1L
  tBase_train <- matrix(1, nrow=nrow(train_data), ncol=1)
  tBase_test  <- matrix(1, nrow=nrow(test_data),  ncol=1)
} else {
  num_t_basis <- min(num_t_basis_req, max(2L, n_time_unique - 1L))
  tBase_train <- splines::bs(train_data$year_mo, df=num_t_basis, intercept=FALSE)
  tBase_test  <- splines::bs(test_data$year_mo,  df=num_t_basis, intercept=FALSE)
}

mBase_full <- mBase[, 1:num_eigen_full, drop=FALSE]
tBase_full_train <- tBase_train[, 1:num_t_basis, drop=FALSE]
tBase_full_test  <- tBase_test[,  1:num_t_basis, drop=FALSE]

spatial_full_train <- as.matrix(AMat_train %*% mBase_full)
spatial_full_test  <- as.matrix(AMat_test  %*% mBase_full)

ST_train_full_GLOBAL <- make_st_scores(spatial_full_train, tBase_full_train)
ST_test_full_GLOBAL  <- make_st_scores(spatial_full_test,  tBase_full_test)

P_full <- ncol(ST_train_full_GLOBAL)

MQM_space_full <- as.matrix(t(mBase_full) %*% Q_lap %*% mBase_full)
MQM_space_full <- 0.5*(MQM_space_full + t(MQM_space_full)) + diag(1e-6, nrow(MQM_space_full))

MQM_full_GLOBAL <- as.matrix(kronecker(MQM_space_full, Diagonal(num_t_basis)))
MQM_full_GLOBAL <- 0.5*(MQM_full_GLOBAL + t(MQM_full_GLOBAL)) + diag(1e-6, nrow(MQM_full_GLOBAL))

constructed_K_values <- constructed_K_values[constructed_K_values <= P_full]
if (length(constructed_K_values) == 0) {
  stop(paste0("All requested constructed_K_values exceed P_full=", P_full,
              ". Reduce constructed_K_values or reduce num_eigen_full/num_t_basis."))
}

# 5) OUTPUT NAMING HELPERS

pdf_file_for_K <- function(K) {
  file.path(out_dir, paste0("SVD_PICAR_REALDATA_", n_label,
                            "_constructed_bases_", K,
                            "_ALL_OUTPUTS_", timestamp, ".pdf"))
}

pdf_file_mcmc <- file.path(out_dir, paste0("SVD_PICAR_REALDATA_", n_label,
                                           "_constructed_bases_200_600_ALL_OUTPUTS_MCMC.pdf"))

csv_ess_allparams_best     <- file.path(out_dir, "ESS_AllParams_AllKbasis.csv")
csv_ess_allparams_allranks <- file.path(out_dir, "ESS_AllParams_AllKbasis_ALL_RANKS.csv")
csv_ess_summary_byK        <- file.path(out_dir, "ESS_Summary_ByKbasis.csv")
csv_ess_summary_byGroup    <- file.path(out_dir, "ESS_Summary_ByGroup_AllKbasis.csv")

pdf_ess_all_Kbasis         <- file.path(out_dir, "06_ess_all_Kbasis.pdf")
pdf_runtime_vs_essps       <- file.path(out_dir, "07_bestRank_runtime_vs_ESSperSec.pdf")

pdf_rank_selection_all_Kbasis <- file.path(out_dir, "02_rank_selection_all_Kbasis.pdf")
pdf_coeff_vs_rank_all_Kbasis  <- file.path(out_dir, "03_coeff_vs_rank_all_Kbasis.pdf")
pdf_best_metrics_vs_Kbasis    <- file.path(out_dir, "04_best_metrics_vs_Kbasis.pdf")

# 6) GLOBAL COLLECTORS

ESS_allparams_ALLRANKS <- data.frame()
ESS_allparams_BEST     <- data.frame()
ESS_summary_byK        <- data.frame()
ESS_summary_byGroup    <- data.frame()

mcmc_perf_by_K <- list()
best_summary_by_K <- data.frame()

all_results <- list()
basis_store <- list()

ranksel_allK_long <- data.frame()
best_ranksel_byK  <- data.frame()

# 7) MAIN LOOP OVER K

for (K in constructed_K_values) {

  pdf_k <- pdf_file_for_K(K)
  pdf(pdf_k, width=11, height=8.5, onefile=TRUE)

  pdf_text_page(paste0("SVD–PICAR (continuous, REAL DATA) — K = ", K),
                c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                  paste0("n_total = ", n_total, ", train_frac = ", train_frac),
                  paste0("num_eigen_full = ", num_eigen_full, ", num_t_basis = ", num_t_basis),
                  "Per-K contents: rank selection (FAST LM), CV, MCMC (best-rank only default).",
                  "TABLE NOTE: beta columns shown are beta1 and beta2 only (exclude beta0)."),
                cex_title=1.05, cex_body=0.9)

  max_combined_rank <- min(K, P_full)

  svd_train <- rsvd::rsvd(ST_train_full_GLOBAL, k=max_combined_rank)
  available_rank <- min(max_combined_rank, length(svd_train$d))
  if (available_rank < max_combined_rank) max_combined_rank <- available_rank

  U_reduced <- svd_train$u[, 1:max_combined_rank, drop=FALSE]
  D_reduced <- diag(svd_train$d[1:max_combined_rank])
  V_reduced <- svd_train$v[, 1:max_combined_rank, drop=FALSE]

  ST_train_reduced <- U_reduced %*% D_reduced
  ST_test_reduced  <- ST_test_full_GLOBAL %*% V_reduced

  colnames(ST_train_reduced) <- paste0("STBase_", seq_len(ncol(ST_train_reduced)))
  colnames(ST_test_reduced)  <- paste0("STBase_", seq_len(ncol(ST_test_reduced)))

  MQM_full <- MQM_full_GLOBAL

  n_train <- length(response_train)
  safe_cap <- safe_max_rank_gaussian(n_train = n_train, pX = pX,
                                     frac = SAFE_RANK_FRAC, buffer = SAFE_RANK_BUFFER)
  r_cap <- min(max_combined_rank, safe_cap)
  rank_grid <- make_rank_grid(r_cap)

  model_results <- data.frame(
    Constructed_ST_Basis=integer(), ST_Rank=integer(), MSPE=numeric(),
    beta0=numeric(), beta0_Lower=numeric(), beta0_Upper=numeric(),
    beta1=numeric(), beta1_Lower=numeric(), beta1_Upper=numeric(),
    beta2=numeric(), beta2_Lower=numeric(), beta2_Upper=numeric(),
    Time_sec=numeric(), stringsAsFactors=FALSE
  )

  Xtr_mat <- as.matrix(covariates_train)
  Xte_mat <- as.matrix(covariates_test)

  for (rnk in rank_grid) {
    start_time <- Sys.time()

    Btr <- as.matrix(ST_train_reduced[, 1:rnk, drop=FALSE])
    Bte <- as.matrix(ST_test_reduced[,  1:rnk, drop=FALSE])

    Ztr <- cbind(Intercept = 1, gws_MEAN = Xtr_mat[,1], rtz_MEAN = Xtr_mat[,2], Btr)
    Zte <- cbind(Intercept = 1, gws_MEAN = Xte_mat[,1], rtz_MEAN = Xte_mat[,2], Bte)

    out_fit <- lmfit_with_ci(Ztr, response_train)
    if (is.null(out_fit)) next

    beta <- out_fit$beta
    lcl  <- out_fit$lcl
    ucl  <- out_fit$ucl

    preds <- as.vector(Zte %*% beta)
    mspe  <- mean((response_test - preds)^2)

    time_sec <- as.numeric(difftime(Sys.time(), start_time, units="secs"))

    model_results <- rbind(model_results, data.frame(
      Constructed_ST_Basis=max_combined_rank, ST_Rank=rnk, MSPE=mspe,
      beta0=unname(beta["Intercept"]), beta0_Lower=lcl["Intercept"], beta0_Upper=ucl["Intercept"],
      beta1=unname(beta["gws_MEAN"]),  beta1_Lower=lcl["gws_MEAN"],  beta1_Upper=ucl["gws_MEAN"],
      beta2=unname(beta["rtz_MEAN"]),  beta2_Lower=lcl["rtz_MEAN"],  beta2_Upper=ucl["rtz_MEAN"],
      Time_sec=time_sec, stringsAsFactors=FALSE
    ))
  }

  if (nrow(model_results) == 0) {
    pdf_text_page(paste0("K = ", K, " FAILED"),
                  c("No successful LM fits over the rank grid.",
                    "Try reducing K, increasing SAFE_RANK_BUFFER, or lowering SAFE_RANK_FRAC."),
                  cex_title=1.05, cex_body=0.9)
    dev.off()
    next
  }

  best_idx  <- which.min(model_results$MSPE)
  best_rank <- model_results$ST_Rank[best_idx]
  min_mspe  <- model_results$MSPE[best_idx]
  model_results$BestRank <- (model_results$ST_Rank == best_rank)

  tmp_rs <- model_results
  tmp_rs$K_requested <- K
  tmp_rs$K_used      <- max_combined_rank
  ranksel_allK_long  <- rbind(ranksel_allK_long, tmp_rs)

  best_ranksel_byK <- rbind(best_ranksel_byK, data.frame(
    Constructed_ST_Basis = K,
    K_used               = max_combined_rank,
    Best_Rank            = best_rank,
    Min_MSPE             = min_mspe,
    beta1                = tmp_rs$beta1[best_idx],
    beta1_Lower          = tmp_rs$beta1_Lower[best_idx],
    beta1_Upper          = tmp_rs$beta1_Upper[best_idx],
    beta2                = tmp_rs$beta2[best_idx],
    beta2_Lower          = tmp_rs$beta2_Lower[best_idx],
    beta2_Upper          = tmp_rs$beta2_Upper[best_idx],
    Time_sec_LM          = tmp_rs$Time_sec[best_idx],
    safe_cap             = safe_cap,
    stringsAsFactors = FALSE
  ))

  write.csv(model_results,
            file.path(out_dir, paste0("RankSelection_REALDATA_", n_label, "_K",K,"_",timestamp,".csv")),
            row.names=FALSE)

  all_results[[as.character(K)]] <- model_results
  basis_store[[as.character(K)]] <- list(
    K=K, K_used=max_combined_rank,
    ST_train_reduced=ST_train_reduced, ST_test_reduced=ST_test_reduced,
    V_reduced=V_reduced, MQM_full=MQM_full, model_results=model_results,
    rank_grid_used=rank_grid, safe_cap=safe_cap
  )

  df_show <- model_results[, c("ST_Rank","MSPE",
                               "beta1","beta1_Lower","beta1_Upper",
                               "beta2","beta2_Lower","beta2_Upper",
                               "Time_sec")]
  df_disp <- df_show
  df_disp$MSPE <- round(df_disp$MSPE,4)
  df_disp$beta1 <- round(df_disp$beta1,3); df_disp$beta1_Lower <- round(df_disp$beta1_Lower,3); df_disp$beta1_Upper <- round(df_disp$beta1_Upper,3)
  df_disp$beta2 <- round(df_disp$beta2,3); df_disp$beta2_Lower <- round(df_disp$beta2_Lower,3); df_disp$beta2_Upper <- round(df_disp$beta2_Upper,3)
  df_disp$Time_sec <- round(df_disp$Time_sec,3)

  tg <- tableGrob(df_disp, rows=NULL)
  st_rank_col <- which(colnames(df_disp)=="ST_Rank")
  tg <- bold_table_cell2(tg, row=best_idx, col=st_rank_col)
  grid.newpage(); grid.draw(tg)

  print(
    ggplot(model_results, aes(x=ST_Rank, y=MSPE)) +
      geom_line() + geom_point(size=2) +
      geom_point(data=subset(model_results, BestRank), size=4) +
      labs(title=paste0("K=",K,": MSPE vs ST Rank (best highlighted)"),
           x="ST Rank", y="MSPE") +
      theme_minimal()
  )

  print(
    ggplot(model_results, aes(x=ST_Rank, y=beta1)) +
      geom_line() + geom_point(size=2) +
      geom_point(data=subset(model_results, BestRank), size=4) +
      labs(title=paste0("K=",K,": beta1 vs ST Rank"),
           x="ST Rank", y="beta1") +
      theme_minimal()
  )

  print(
    ggplot(model_results, aes(x=ST_Rank, y=beta2)) +
      geom_line() + geom_point(size=2) +
      geom_point(data=subset(model_results, BestRank), size=4) +
      labs(title=paste0("K=",K,": beta2 vs ST Rank"),
           x="ST Rank", y="beta2") +
      theme_minimal()
  )

  if (RUN_CV_BLOCK) {

    pdf_text_page(paste0("CROSS-VALIDATION (CV) — K = ", K),
                  c("Refit FAST lm.fit() and evaluate MSPE on held-out TEST set.",
                    paste0("CV setting: RUN_CV_BEST_RANK_ONLY = ", RUN_CV_BEST_RANK_ONLY),
                    "CV table includes beta1 and beta2 only (exclude beta0)."),
                  cex_title=1.05, cex_body=0.9)

    cv_rank_grid <- if (isTRUE(RUN_CV_BEST_RANK_ONLY)) best_rank else rank_grid

    cv_results <- data.frame(
      Constructed_ST_Basis=integer(), ST_Rank=integer(), MSPE=numeric(),
      beta0=numeric(), beta0_Lower=numeric(), beta0_Upper=numeric(),
      beta1=numeric(), beta1_Lower=numeric(), beta1_Upper=numeric(),
      beta2=numeric(), beta2_Lower=numeric(), beta2_Upper=numeric(),
      Time_sec=numeric(), stringsAsFactors=FALSE
    )

    for (current_rank in cv_rank_grid) {
      start_time <- Sys.time()

      Btr <- as.matrix(ST_train_reduced[, 1:current_rank, drop=FALSE])
      Bte <- as.matrix(ST_test_reduced[,  1:current_rank, drop=FALSE])

      Ztr <- cbind(Intercept = 1, gws_MEAN = Xtr_mat[,1], rtz_MEAN = Xtr_mat[,2], Btr)
      Zte <- cbind(Intercept = 1, gws_MEAN = Xte_mat[,1], rtz_MEAN = Xte_mat[,2], Bte)

      out_fit <- lmfit_with_ci(Ztr, response_train)
      if (is.null(out_fit)) next

      beta <- out_fit$beta
      lcl  <- out_fit$lcl
      ucl  <- out_fit$ucl

      preds <- as.vector(Zte %*% beta)
      MSPE  <- mean((response_test - preds)^2)

      time_sec <- as.numeric(difftime(Sys.time(), start_time, units="secs"))

      cv_results <- rbind(cv_results, data.frame(
        Constructed_ST_Basis=max_combined_rank, ST_Rank=current_rank, MSPE=MSPE,
        beta0=unname(beta["Intercept"]), beta0_Lower=lcl["Intercept"], beta0_Upper=ucl["Intercept"],
        beta1=unname(beta["gws_MEAN"]),  beta1_Lower=lcl["gws_MEAN"],  beta1_Upper=ucl["gws_MEAN"],
        beta2=unname(beta["rtz_MEAN"]),  beta2_Lower=lcl["rtz_MEAN"],  beta2_Upper=ucl["rtz_MEAN"],
        Time_sec=time_sec, stringsAsFactors=FALSE
      ))
    }

    write.csv(cv_results,
              file.path(out_dir, paste0("CrossValidation_REALDATA_", n_label, "_K", K, "_", timestamp, ".csv")),
              row.names = FALSE)

    if (nrow(cv_results) > 0) {
      best_cv_idx  <- which.min(cv_results$MSPE)
      best_cv_rank <- cv_results$ST_Rank[best_cv_idx]
      cv_results$BestRank <- (cv_results$ST_Rank == best_cv_rank)

      df_cv_show <- cv_results[, c("ST_Rank","MSPE",
                                   "beta1","beta1_Lower","beta1_Upper",
                                   "beta2","beta2_Lower","beta2_Upper",
                                   "Time_sec")]
      df_cv_disp <- df_cv_show
      df_cv_disp$MSPE <- round(df_cv_disp$MSPE,4)
      df_cv_disp$beta1 <- round(df_cv_disp$beta1,3); df_cv_disp$beta1_Lower <- round(df_cv_disp$beta1_Lower,3); df_cv_disp$beta1_Upper <- round(df_cv_disp$beta1_Upper,3)
      df_cv_disp$beta2 <- round(df_cv_disp$beta2,3); df_cv_disp$beta2_Lower <- round(df_cv_disp$beta2_Lower,3); df_cv_disp$beta2_Upper <- round(df_cv_disp$beta2_Upper,3)
      df_cv_disp$Time_sec <- round(df_cv_disp$Time_sec,3)

      tg2 <- tableGrob(df_cv_disp, rows=NULL)
      st_rank_col2 <- which(colnames(df_cv_disp)=="ST_Rank")
      tg2 <- bold_table_cell2(tg2, row=best_cv_idx, col=st_rank_col2)
      grid.newpage(); grid.draw(tg2)

      print(
        ggplot(cv_results, aes(x=ST_Rank, y=MSPE)) +
          geom_line() + geom_point(size=2) +
          geom_point(data=subset(cv_results, BestRank), size=4) +
          labs(title=paste0("CV (K=",K,"): MSPE vs ST Rank (best highlighted)"),
               x="ST Rank", y="MSPE") +
          theme_minimal()
      )
    } else {
      pdf_text_page("CV BLOCK NOTE",
                    c("No successful CV fits (cv_results empty).",
                      "This can happen if the best rank is near the safe cap or LM becomes ill-conditioned."),
                    cex_title=1.05, cex_body=0.9)
    }
  }

  if (RUN_MCMC_BLOCK) {

    pdf_text_page(paste0("MCMC PERFORMANCE (K = ", K, ")"),
                  c(paste0("MCMC settings: niter = ", niter_mcmc_rank,
                           ", burn-in = ", round(100*burn_prop_rank), "%"),
                    paste0("RUN_MCMC_ALL_RANKS = ", RUN_MCMC_ALL_RANKS, " (FALSE = best-rank only)"),
                    "Runtime = sampling-only; compile time tracked separately.",
                    "Tables focus on beta1 and beta2 only (exclude beta0)."),
                  cex_title=1.05, cex_body=0.9)

    X_matrix <- as.matrix(covariates_train)
    y_train  <- as.numeric(response_train)

    mcmc_rank_grid <- if (isTRUE(RUN_MCMC_ALL_RANKS)) rank_grid else best_rank

    mcmc_rank_perf <- data.frame(
      Constructed_ST_Basis = integer(),
      ST_Rank              = integer(),
      compile_time_sec     = numeric(),
      runtime_sec          = numeric(),

      beta0_ESS            = numeric(), beta1_ESS = numeric(), beta2_ESS = numeric(),
      beta0_ESS_per_sec    = numeric(), beta1_ESS_per_sec = numeric(), beta2_ESS_per_sec = numeric(),

      ESS_per_sec_mean_beta   = numeric(),
      ESS_per_sec_mean_beta12 = numeric(),
      stringsAsFactors = FALSE
    )

    for (rnk in mcmc_rank_grid) {
      cat("MCMC: K =", K, "| Rank =", rnk, "\n")

      M_matrix <- as.matrix(ST_train_reduced[, 1:rnk, drop=FALSE])

      V_opt <- V_reduced[, 1:rnk, drop=FALSE]
      MQM_reduced <- as.matrix(t(V_opt) %*% MQM_full %*% V_opt)
      MQM_reduced <- 0.5*(MQM_reduced + t(MQM_reduced)) + diag(1e-6, nrow(MQM_reduced))

      out_mcmc <- run_picar_mcmc_one_rank(
        y_train=y_train, X_matrix=X_matrix,
        M_matrix=M_matrix, MQM_reduced=MQM_reduced,
        niter=niter_mcmc_rank, burn_prop=burn_prop_rank,
        monitor_delta = MONITOR_DELTA
      )

      mcmc_rank_perf <- rbind(mcmc_rank_perf, data.frame(
        Constructed_ST_Basis = K,
        ST_Rank              = rnk,
        compile_time_sec     = out_mcmc$compile_time_sec,
        runtime_sec          = out_mcmc$runtime_sec,

        beta0_ESS            = out_mcmc$beta0_ESS,
        beta1_ESS            = out_mcmc$beta1_ESS,
        beta2_ESS            = out_mcmc$beta2_ESS,

        beta0_ESS_per_sec    = out_mcmc$beta0_ESS_per_sec,
        beta1_ESS_per_sec    = out_mcmc$beta1_ESS_per_sec,
        beta2_ESS_per_sec    = out_mcmc$beta2_ESS_per_sec,

        ESS_per_sec_mean_beta   = out_mcmc$ESS_per_sec_mean_beta,
        ESS_per_sec_mean_beta12 = out_mcmc$ESS_per_sec_mean_beta12,
        stringsAsFactors = FALSE
      ))

      tmp_long <- out_mcmc$ess_long
      tmp_long$Constructed_ST_Basis <- K
      tmp_long$ST_Rank <- rnk
      tmp_long$compile_time_sec <- out_mcmc$compile_time_sec
      tmp_long$runtime_sec <- out_mcmc$runtime_sec
      ESS_allparams_ALLRANKS <- rbind(ESS_allparams_ALLRANKS, tmp_long)
    }

    write.csv(mcmc_rank_perf,
              file.path(out_dir, paste0("MCMC_Perf_REALDATA_", n_label, "_K", K, "_", timestamp, ".csv")),
              row.names = FALSE)

    mcmc_best_row <- mcmc_rank_perf[mcmc_rank_perf$ST_Rank == best_rank, , drop=FALSE]
    if (nrow(mcmc_best_row) == 0) mcmc_best_row <- mcmc_rank_perf[1, , drop=FALSE]

    mcmc_perf_by_K[[as.character(K)]] <- mcmc_rank_perf

    best_summary_by_K <- rbind(best_summary_by_K, data.frame(
      Constructed_ST_Basis = K,
      Best_Rank            = best_rank,
      Min_MSPE             = min_mspe,
      BestRank_Time_sec    = mcmc_best_row$runtime_sec[1],
      BestRank_Compile_sec = mcmc_best_row$compile_time_sec[1],

      beta0_ESS            = mcmc_best_row$beta0_ESS[1],
      beta1_ESS            = mcmc_best_row$beta1_ESS[1],
      beta2_ESS            = mcmc_best_row$beta2_ESS[1],

      beta0_ESS_per_sec    = mcmc_best_row$beta0_ESS_per_sec[1],
      beta1_ESS_per_sec    = mcmc_best_row$beta1_ESS_per_sec[1],
      beta2_ESS_per_sec    = mcmc_best_row$beta2_ESS_per_sec[1],

      ESS_per_sec_mean_beta   = mcmc_best_row$ESS_per_sec_mean_beta[1],
      ESS_per_sec_mean_beta12 = mcmc_best_row$ESS_per_sec_mean_beta12[1],
      stringsAsFactors = FALSE
    ))

    essK_all  <- ESS_allparams_ALLRANKS[ESS_allparams_ALLRANKS$Constructed_ST_Basis == K, , drop=FALSE]
    essK_best <- essK_all[essK_all$ST_Rank == best_rank, , drop=FALSE]

    if (nrow(essK_best) > 0) {
      ESS_allparams_BEST <- rbind(ESS_allparams_BEST, essK_best)

      rt_best <- unique(essK_best$runtime_sec)[1]

      essps_all    <- essK_best$ESS_per_sec
      essps_beta   <- essK_best$ESS_per_sec[essK_best$Parameter %in% c("beta0","beta1","beta2")]
      essps_beta12 <- essK_best$ESS_per_sec[essK_best$Parameter %in% c("beta1","beta2")]

      ESS_summary_byK <- rbind(ESS_summary_byK, data.frame(
        Constructed_ST_Basis = K,
        Best_Rank            = best_rank,
        Min_MSPE             = min_mspe,
        BestRank_RuntimeSec  = rt_best,

        ESSps_mean_all       = mean(essps_all,    na.rm=TRUE),
        ESSps_median_all     = median(essps_all,  na.rm=TRUE),
        ESSps_min_all        = min(essps_all,     na.rm=TRUE),

        ESSps_mean_beta012   = mean(essps_beta,   na.rm=TRUE),
        ESSps_min_beta012    = min(essps_beta,    na.rm=TRUE),

        ESSps_mean_beta12    = mean(essps_beta12, na.rm=TRUE),
        ESSps_min_beta12     = min(essps_beta12,  na.rm=TRUE),

        n_params_bestRank    = nrow(essK_best),
        stringsAsFactors = FALSE
      ))

      tmpG <- aggregate(
        cbind(ESS, ESS_per_sec) ~ Group,
        data = essK_best,
        FUN = function(x) c(mean=mean(x,na.rm=TRUE), median=median(x,na.rm=TRUE), min=min(x,na.rm=TRUE))
      )

      tmp_out <- data.frame()
      for (i in 1:nrow(tmpG)) {
        g <- tmpG$Group[i]
        ESS_stats   <- tmpG$ESS[i, ]
        ESSps_stats <- tmpG$ESS_per_sec[i, ]
        tmp_out <- rbind(tmp_out, data.frame(
          Constructed_ST_Basis = K,
          Best_Rank            = best_rank,
          Group                = g,
          RuntimeSec           = rt_best,

          ESS_mean             = ESS_stats["mean"],
          ESS_median           = ESS_stats["median"],
          ESS_min              = ESS_stats["min"],

          ESSps_mean           = ESSps_stats["mean"],
          ESSps_median         = ESSps_stats["median"],
          ESSps_min            = ESSps_stats["min"],

          stringsAsFactors = FALSE
        ))
      }
      ESS_summary_byGroup <- rbind(ESS_summary_byGroup, tmp_out)
    }
  }

  pdf_text_page("DONE (this K only)",
                c(paste0("Per-K PDF saved to: ", pdf_k),
                  paste0("Best rank: ", best_rank, " | Min MSPE: ", signif(min_mspe, 6)),
                  paste0("Rank grid used: ", paste(rank_grid, collapse=", "))),
                cex_title=1.05, cex_body=0.9)

  dev.off()
}

# 8) GLOBAL PDFs

pdf(pdf_rank_selection_all_Kbasis, width=11, height=8.5, onefile=TRUE)
pdf_text_page("02 — RANK SELECTION (ALL Kbasis) — CONTINUOUS SVD–PICAR (REAL DATA)",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                "Aggregated rank-selection results across all K.",
                "TABLE NOTE: beta columns shown are beta1 and beta2 only (exclude beta0)."),
              cex_title=1.05, cex_body=0.9)

if (nrow(best_ranksel_byK) > 0) {
  tab <- best_ranksel_byK
  tab$Min_MSPE <- signif(tab$Min_MSPE, 6)
  tab$beta1 <- round(tab$beta1, 3); tab$beta1_Lower <- round(tab$beta1_Lower, 3); tab$beta1_Upper <- round(tab$beta1_Upper, 3)
  tab$beta2 <- round(tab$beta2, 3); tab$beta2_Lower <- round(tab$beta2_Lower, 3); tab$beta2_Upper <- round(tab$beta2_Upper, 3)
  tab$Time_sec_LM <- round(tab$Time_sec_LM, 3)

  tab_show <- tab[, c("Constructed_ST_Basis","K_used","Best_Rank","Min_MSPE",
                      "beta1","beta1_Lower","beta1_Upper",
                      "beta2","beta2_Lower","beta2_Upper",
                      "Time_sec_LM","safe_cap")]

  tg <- tableGrob(tab_show, rows=NULL)
  grid.newpage(); grid.draw(tg)

  print(
    ggplot(best_ranksel_byK, aes(x=Constructed_ST_Basis, y=Min_MSPE)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(title="Best (min MSPE) vs Kbasis", x="K", y="Min MSPE")
  )

  print(
    ggplot(best_ranksel_byK, aes(x=Constructed_ST_Basis, y=Best_Rank)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(title="Selected best rank vs Kbasis", x="K", y="Best Rank")
  )
} else {
  pdf_text_page("NOTE", c("best_ranksel_byK is empty (no successful rank-selection fits)."),
                cex_title=1.05, cex_body=0.9)
}
dev.off()

pdf(pdf_coeff_vs_rank_all_Kbasis, width=11, height=8.5, onefile=TRUE)
pdf_text_page("03 — COEFFICIENTS vs RANK (ALL Kbasis) — CONTINUOUS SVD–PICAR (REAL DATA)",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                "Shows how beta1 and beta2 vary with ST rank across K."),
              cex_title=1.05, cex_body=0.9)

if (nrow(ranksel_allK_long) > 0) {
  print(
    ggplot(ranksel_allK_long, aes(x=ST_Rank, y=beta1, group=factor(K_requested))) +
      geom_line() + geom_point(size=1.6) +
      facet_wrap(~K_requested, scales="free_x") +
      theme_minimal() +
      labs(title="beta1 vs ST Rank (faceted by K)", x="ST Rank", y="beta1")
  )
  print(
    ggplot(ranksel_allK_long, aes(x=ST_Rank, y=beta2, group=factor(K_requested))) +
      geom_line() + geom_point(size=1.6) +
      facet_wrap(~K_requested, scales="free_x") +
      theme_minimal() +
      labs(title="beta2 vs ST Rank (faceted by K)", x="ST Rank", y="beta2")
  )
} else {
  pdf_text_page("NOTE", c("ranksel_allK_long is empty."),
                cex_title=1.05, cex_body=0.9)
}
dev.off()

pdf(pdf_best_metrics_vs_Kbasis, width=11, height=8.5, onefile=TRUE)
pdf_text_page("04 — BEST METRICS vs Kbasis — CONTINUOUS SVD–PICAR (REAL DATA)",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                "Best-rank selection metrics (FAST LM) vs K.",
                "If MCMC ran, includes best-rank runtime and ESS/sec(beta1,beta2)."),
              cex_title=1.05, cex_body=0.9)

if (nrow(best_ranksel_byK) > 0) {
  bm <- best_ranksel_byK
  if (nrow(best_summary_by_K) > 0) {
    mcmc_best <- best_summary_by_K[, c("Constructed_ST_Basis","Best_Rank",
                                       "BestRank_Time_sec","BestRank_Compile_sec",
                                       "ESS_per_sec_mean_beta12")]
    bm <- merge(bm, mcmc_best, by=c("Constructed_ST_Basis","Best_Rank"), all.x=TRUE)
  } else {
    bm$BestRank_Time_sec <- NA_real_
    bm$BestRank_Compile_sec <- NA_real_
    bm$ESS_per_sec_mean_beta12 <- NA_real_
  }

  tab <- bm
  tab$Min_MSPE <- signif(tab$Min_MSPE, 6)
  tab$Time_sec_LM <- round(tab$Time_sec_LM, 3)
  tab$BestRank_Time_sec <- round(tab$BestRank_Time_sec, 2)
  tab$BestRank_Compile_sec <- round(tab$BestRank_Compile_sec, 2)
  tab$ESS_per_sec_mean_beta12 <- round(tab$ESS_per_sec_mean_beta12, 4)

  tab_show <- tab[, c("Constructed_ST_Basis","K_used","Best_Rank","Min_MSPE",
                      "Time_sec_LM","BestRank_Time_sec","BestRank_Compile_sec",
                      "ESS_per_sec_mean_beta12","safe_cap")]

  tg <- tableGrob(tab_show, rows=NULL)
  grid.newpage(); grid.draw(tg)

  print(
    ggplot(bm, aes(x=Constructed_ST_Basis, y=Min_MSPE)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(title="Best MSPE vs Kbasis", x="K", y="Min MSPE")
  )

  if (any(is.finite(bm$BestRank_Time_sec))) {
    print(
      ggplot(bm, aes(x=Constructed_ST_Basis, y=BestRank_Time_sec)) +
        geom_line() + geom_point(size=2) +
        theme_minimal() +
        labs(title="MCMC runtime (sampling-only) at best rank vs Kbasis", x="K", y="Runtime (sec)")
    )
  }
} else {
  pdf_text_page("NOTE", c("best_ranksel_byK is empty; cannot summarize best metrics."),
                cex_title=1.05, cex_body=0.9)
}
dev.off()

# 9) GLOBAL ESS OUTPUTS + GLOBAL MCMC PDF

write.csv(ESS_allparams_BEST,     csv_ess_allparams_best,     row.names=FALSE)
write.csv(ESS_allparams_ALLRANKS, csv_ess_allparams_allranks, row.names=FALSE)
write.csv(ESS_summary_byK,        csv_ess_summary_byK,        row.names=FALSE)
write.csv(ESS_summary_byGroup,    csv_ess_summary_byGroup,    row.names=FALSE)

pdf(pdf_ess_all_Kbasis, width=11, height=8.5, onefile=TRUE)
pdf_text_page("GLOBAL ESS SUMMARY (BEST RANK PER K) — REAL DATA",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                paste0("RUN_MCMC_ALL_RANKS = ", RUN_MCMC_ALL_RANKS)),
              cex_title=1.05, cex_body=0.9)

if (nrow(ESS_summary_byK) > 0) {
  df_tab <- ESS_summary_byK
  df_tab$Min_MSPE <- signif(df_tab$Min_MSPE, 6)
  df_tab$BestRank_RuntimeSec <- round(df_tab$BestRank_RuntimeSec, 2)
  df_tab$ESSps_mean_all <- round(df_tab$ESSps_mean_all, 4)
  df_tab$ESSps_min_all <- round(df_tab$ESSps_min_all, 4)
  df_tab$ESSps_mean_beta12 <- round(df_tab$ESSps_mean_beta12, 4)
  df_tab$ESSps_min_beta12 <- round(df_tab$ESSps_min_beta12, 4)
  tg <- tableGrob(df_tab, rows=NULL)
  grid.newpage(); grid.draw(tg)

  print(
    ggplot(ESS_summary_byK, aes(x=Constructed_ST_Basis, y=ESSps_mean_beta12)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(title="Best-rank mean ESS/sec (beta1,beta2) vs K",
           x="K", y="Mean ESS/sec (beta1,beta2)")
  )
}
dev.off()

pdf(pdf_runtime_vs_essps, width=10, height=7, onefile=TRUE)
pdf_text_page("BEST-RANK RUNTIME vs ESS/sec (REAL DATA)",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                "Runtime is sampling-only; compile time tracked separately."),
              cex_title=1.05, cex_body=0.9)

if (nrow(ESS_summary_byK) > 0) {
  print(
    ggplot(ESS_summary_byK, aes(x=BestRank_RuntimeSec, y=ESSps_mean_beta12)) +
      geom_point(size=2) +
      geom_text(aes(label=Constructed_ST_Basis), vjust=-0.8, size=3) +
      theme_minimal() +
      labs(title="Best-rank runtime vs mean ESS/sec (beta1,beta2) — labels are K",
           x="Runtime (sec; sampling only)", y="Mean ESS/sec (beta1,beta2)")
  )
}
dev.off()

mcmc_perf_all <- do.call(rbind, lapply(names(mcmc_perf_by_K), function(kstr) {
  df <- mcmc_perf_by_K[[kstr]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$Constructed_ST_Basis <- as.integer(kstr)
  df
}))

pdf(pdf_file_mcmc, width=11, height=8.5, onefile=TRUE)
pdf_text_page("MCMC SUMMARY + COMBINED CURVES (ALL K) — REAL DATA",
              c(paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                paste0("n_total = ", n_total),
                paste0("Per-rank MCMC: niter=", niter_mcmc_rank, ", burn=", round(100*burn_prop_rank), "%"),
                paste0("RUN_MCMC_ALL_RANKS = ", RUN_MCMC_ALL_RANKS)),
              cex_title=1.05, cex_body=0.9)

if (nrow(best_summary_by_K) > 0) {
  best_sum_disp <- best_summary_by_K
  best_sum_disp$Min_MSPE <- signif(best_sum_disp$Min_MSPE, 6)
  best_sum_disp$BestRank_Time_sec <- round(best_sum_disp$BestRank_Time_sec, 2)
  best_sum_disp$BestRank_Compile_sec <- round(best_sum_disp$BestRank_Compile_sec, 2)
  best_sum_disp$beta1_ESS <- round(best_sum_disp$beta1_ESS, 1)
  best_sum_disp$beta2_ESS <- round(best_sum_disp$beta2_ESS, 1)
  best_sum_disp$beta1_ESS_per_sec <- round(best_sum_disp$beta1_ESS_per_sec, 4)
  best_sum_disp$beta2_ESS_per_sec <- round(best_sum_disp$beta2_ESS_per_sec, 4)
  best_sum_disp$ESS_per_sec_mean_beta12 <- round(best_sum_disp$ESS_per_sec_mean_beta12, 4)

  best_sum_tbl <- best_sum_disp[, c("Constructed_ST_Basis","Best_Rank","Min_MSPE",
                                    "BestRank_Time_sec","BestRank_Compile_sec",
                                    "beta1_ESS","beta2_ESS",
                                    "beta1_ESS_per_sec","beta2_ESS_per_sec",
                                    "ESS_per_sec_mean_beta12")]
  tg_best <- tableGrob(best_sum_tbl, rows=NULL)
  grid.newpage(); grid.draw(tg_best)
}

if (!is.null(mcmc_perf_all) && nrow(mcmc_perf_all) > 0) {
  print(
    ggplot(mcmc_perf_all,
           aes(x=ST_Rank, y=runtime_sec,
               group=factor(Constructed_ST_Basis),
               color=factor(Constructed_ST_Basis))) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(title="ST_Rank vs Time (MCMC runtime; sampling only) — All K",
           x="ST Rank", y="Time (sec)", color="K")
  )
}
dev.off()

# 10) SAVE EVERYTHING

save(
  n_total, train_frac,
  beta0_true, beta1_true, beta2_true, tau_true,
  constructed_K_values, num_eigen_full, num_t_basis,
  mesh, AMat_train, AMat_test, W_adj, Q_lap, mBase, tBase_train, tBase_test,
  train_data, test_data,
  all_results, basis_store, mcmc_perf_by_K, best_summary_by_K,
  ESS_allparams_ALLRANKS, ESS_allparams_BEST, ESS_summary_byK, ESS_summary_byGroup,
  file = file.path(out_dir, paste0("All_Objects_Continuous_SVD_PICAR_REALDATA_FAST_", n_label, "_", timestamp, ".RData"))
)

overall_end <- Sys.time()
cat(sprintf("\nDONE. Total runtime: %.2f sec\nOutputs in: %s\n",
            as.numeric(difftime(overall_end, overall_start, units="secs")),
            out_dir))

cat("\nKey global outputs:\n",
    "  ", csv_ess_allparams_best, "\n",
    "  ", csv_ess_allparams_allranks, "\n",
    "  ", csv_ess_summary_byK, "\n",
    "  ", csv_ess_summary_byGroup, "\n",
    "  ", pdf_ess_all_Kbasis, "\n",
    "  ", pdf_runtime_vs_essps, "\n",
    "  ", pdf_file_mcmc, "\n")

# END REAL DATA SCRIPT

# 12) NEW ANALYSIS PDF (REAL DATA VERSION)

pdf_new_analysis <- file.path(out_dir, "08_SVD_PICAR_Oversmoothing_and_MeshDecoupling_REALDATA.pdf")
csv_mesh_sens    <- file.path(out_dir, "MeshDensity_Sensitivity_SVD_PICAR_REALDATA.csv")
csv_rank_sweep   <- file.path(out_dir, "Roughness_ByRank_SVD_PICAR_Gaussian_REALDATA.csv")

# Packages used ONLY in this analysis block
need_pkgs <- c("grid", "gridExtra", "ggplot2")
for (pp in need_pkgs) {
  if (!requireNamespace(pp, quietly = TRUE)) install.packages(pp)
}
library(grid)
library(gridExtra)
library(ggplot2)

# Helper: "title + bullets" PDF page
pdf_text_page_12 <- function(title, bullets, cex_title=1.1, cex_body=0.95) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 9), c("null","null")))))
  grid.text(title, vp = viewport(layout.pos.row=1, layout.pos.col=1),
            gp=gpar(fontsize=16*cex_title, fontface="bold"))
  grid.text(paste0("\u2022 ", bullets, collapse="\n"),
            x=unit(0.02,"npc"), y=unit(0.98,"npc"),
            just=c("left","top"),
            vp = viewport(layout.pos.row=2, layout.pos.col=1),
            gp=gpar(fontsize=12*cex_body))
  popViewport()
}

# Helper: bs() evaluation at new time using train basis attributes
bs_eval_like_train <- function(x_new, bs_train_mat) {
  if (ncol(bs_train_mat) == 1L && isTRUE(all(bs_train_mat == 1))) {
    return(matrix(1, nrow=length(x_new), ncol=1))
  }
  splines::bs(
    x_new,
    df             = ncol(bs_train_mat),
    intercept      = FALSE,
    knots          = attr(bs_train_mat, "knots"),
    Boundary.knots = attr(bs_train_mat, "Boundary.knots")
  )
}

# Helper: unique undirected edge list from adjacency W
get_edge_list_from_W <- function(W_sparse) {
  ij <- which(W_sparse != 0, arr.ind = TRUE)
  ij <- ij[ij[,1] < ij[,2], , drop=FALSE]
  colnames(ij) <- c("i","j")
  ij
}

# Helper: edge roughness index R(u)=mean_{edges}(u_i-u_j)^2
roughness_index_edges <- function(u, edges_ij) {
  dif2 <- (u[edges_ij[,"i"]] - u[edges_ij[,"j"]])^2
  mean(dif2, na.rm = TRUE)
}

# Helper: compute SVD–PICAR latent field u on mesh vertices at fixed time
compute_u_mesh_at_tstar <- function(mesh_obj, mBase_full, tBase_train,
                                    V_reduced, rnk, delta_bar, t_star) {

  t_row <- bs_eval_like_train(t_star, tBase_train)
  if (is.vector(t_row)) t_row <- matrix(t_row, nrow=1)
  if (nrow(t_row) != 1L) t_row <- matrix(t_row[1, ], nrow=1)

  Nmesh <- nrow(mesh_obj$loc)
  temporal_mesh <- matrix(rep(t_row, each=Nmesh), nrow=Nmesh, byrow=FALSE)

  ST_mesh_full <- make_st_scores(mBase_full, temporal_mesh)

  ST_mesh_red  <- as.matrix(ST_mesh_full %*% V_reduced)
  M_mesh       <- ST_mesh_red[, 1:rnk, drop=FALSE]

  as.vector(M_mesh %*% delta_bar)
}

# Choose a representative time slice t*
if (exists("has_time") && isTRUE(has_time)) {
  t_star <- median(train_data$year_mo)
} else {
  t_star <- 1
}

# Edge list on the ORIGINAL mesh
edges_ij <- get_edge_list_from_W(W_adj)

# Rank sweep settings (REAL DATA)
K_vals_sweep <- constructed_K_values
niter_rank_sweep <- 30000
burn_rank_sweep  <- 0.60

make_rank_sweep_grid <- function(K_used) {
  base <- unique(sort(c(5, 10, 15, 20, 30, 40, 50, 75, 100)))
  base <- base[base <= K_used]
  if (length(base) == 0) base <- unique(sort(c(1, min(5, K_used), K_used)))
  if (K_used > max(base)) base <- unique(sort(c(base, min(K_used, 150), K_used)))
  base
}

# Storage for rank-sweep results
rank_sweep_res <- data.frame()

cat("\n[NEW ANALYSIS — REAL DATA] Part 1: fitting MANY ranks for each K (Gaussian SVD-PICAR)...\n")

for (K_map in K_vals_sweep) {

  cat("\n  --- K =", K_map, " ---\n")

  bsK <- basis_store[[as.character(K_map)]]
  if (is.null(bsK)) stop(paste0("basis_store missing K=", K_map))

  ST_train_red <- as.matrix(bsK$ST_train_reduced)
  V_reduced    <- as.matrix(bsK$V_reduced)
  MQM_full     <- as.matrix(bsK$MQM_full)

  K_used <- ncol(ST_train_red)
  if (!is.null(bsK$K_used)) K_used <- as.integer(bsK$K_used)

  r_grid <- make_rank_sweep_grid(K_used)

  r_best <- NA_integer_
  if (exists("best_ranksel_byK") && nrow(best_ranksel_byK) > 0) {
    rr <- best_ranksel_byK[best_ranksel_byK$Constructed_ST_Basis == K_map, , drop=FALSE]
    if (nrow(rr) > 0) r_best <- as.integer(rr$Best_Rank[1])
  }

  for (rnk in r_grid) {

    cat("    fitting rank r =", rnk, "...\n")

    V_opt  <- V_reduced[, 1:rnk, drop=FALSE]
    MQM_red <- as.matrix(t(V_opt) %*% MQM_full %*% V_opt)
    MQM_red <- 0.5*(MQM_red + t(MQM_red)) + diag(1e-6, nrow(MQM_red))

    M_train <- ST_train_red[, 1:rnk, drop=FALSE]

    out_fit <- run_picar_mcmc_one_rank(
      y_train       = as.numeric(response_train),
      X_matrix      = as.matrix(covariates_train),
      M_matrix      = M_train,
      MQM_reduced   = MQM_red,
      niter         = niter_rank_sweep,
      burn_prop     = burn_rank_sweep,
      monitor_delta = TRUE
    )

    delta_cols <- grep("^delta\\[", colnames(out_fit$post), value=TRUE)
    if (length(delta_cols) < 1) stop("No delta columns found in MCMC output (monitor_delta=TRUE expected).")

    delta_bar <- colMeans(out_fit$post[, delta_cols, drop=FALSE])

    u_hat_mesh <- compute_u_mesh_at_tstar(
      mesh_obj    = mesh,
      mBase_full  = mBase_full,
      tBase_train = tBase_train,
      V_reduced   = V_reduced[, 1:K_used, drop=FALSE],
      rnk         = rnk,
      delta_bar   = delta_bar[1:rnk],
      t_star      = t_star
    )

    R_hat  <- roughness_index_edges(u_hat_mesh, edges_ij)
    sd_hat <- sd(u_hat_mesh)

    sf_hat <- sqrt(R_hat) / max(sd_hat, 1e-12)

    rank_sweep_res <- rbind(rank_sweep_res, data.frame(
      K_map         = K_map,
      K_used        = K_used,
      rank          = rnk,
      runtime_sec   = out_fit$runtime_sec,
      time_per_iter = out_fit$runtime_sec / niter_rank_sweep,
      R_hat         = R_hat,
      sd_hat        = sd_hat,
      sf_hat        = sf_hat,
      is_best_rank  = if (!is.na(r_best)) (rnk == r_best) else FALSE,
      stringsAsFactors = FALSE
    ))
  }
}

# Within-K normalization: ratio to the *largest rank fitted* for that K
if (nrow(rank_sweep_res) > 0) {
  ref_tbl <- do.call(rbind, lapply(split(rank_sweep_res, rank_sweep_res$K_map), function(dfK) {
    dfK <- dfK[order(dfK$rank), , drop=FALSE]
    ref <- dfK[nrow(dfK), , drop=FALSE]
    data.frame(
      K_map   = ref$K_map[1],
      R_ref   = ref$R_hat[1],
      sf_ref  = ref$sf_hat[1],
      rank_ref= ref$rank[1],
      stringsAsFactors = FALSE
    )
  }))
  rank_sweep_res <- merge(rank_sweep_res, ref_tbl, by="K_map", all.x=TRUE)

  rank_sweep_res$R_ratio_ref  <- rank_sweep_res$R_hat / pmax(rank_sweep_res$R_ref, 1e-12)
  rank_sweep_res$sf_ratio_ref <- rank_sweep_res$sf_hat / pmax(rank_sweep_res$sf_ref, 1e-12)
}

write.csv(rank_sweep_res, csv_rank_sweep, row.names=FALSE)

# Part 2) Mesh-density decoupling experiment
K_sens <- max(constructed_K_values)

r_sens_target <- NA_integer_
tmpK <- rank_sweep_res[rank_sweep_res$K_map == K_sens, , drop=FALSE]
if (nrow(tmpK) > 0) {
  rr_best <- tmpK[tmpK$is_best_rank, , drop=FALSE]
  if (nrow(rr_best) > 0) r_sens_target <- rr_best$rank[1]
  if (is.na(r_sens_target)) r_sens_target <- tmpK$rank[ceiling(nrow(tmpK)/2)]
} else {
  r_sens_target <- 20
}

niter_sens <- 8000
burn_sens  <- 0.50
max_edge_grid <- c(2.5, 2.0, 1.5, 1.0)

mesh_sens_res <- data.frame()

cat("\n[NEW ANALYSIS — REAL DATA] Part 2: mesh-density sensitivity experiment...\n")

for (me in max_edge_grid) {

  cat("  mesh max.edge0 =", me, "\n")

  t0_build <- proc.time()

  mesh_tmp <- INLA::inla.mesh.2d(
    loc      = train_coords,
    max.edge = c(1,2) * me,
    cutoff   = me/5,
    offset   = c(me, 6.0)
  )

  Nmesh_tmp <- nrow(mesh_tmp$loc)
  A_tr_tmp <- INLA::inla.spde.make.A(mesh_tmp, loc = as.matrix(train_data[,c("lon","lat")]))

  lap_tmp <- mesh_to_laplacian(mesh_tmp)
  W_tmp   <- lap_tmp$W
  Q_tmp   <- as.matrix(lap_tmp$Q) + diag(1e-6, lap_tmp$Nnew)

  num_eigen_tmp <- min(num_eigen_full_req, max(5L, lap_tmp$Nnew - 2L))
  mBase_tmp <- build_moran_basis(W_tmp, num_eigen_tmp)

  spatial_tr_tmp  <- as.matrix(A_tr_tmp %*% mBase_tmp)
  temporal_tr_tmp <- tBase_full_train
  ST_tr_full_tmp  <- make_st_scores(spatial_tr_tmp, temporal_tr_tmp)
  P_full_tmp      <- ncol(ST_tr_full_tmp)

  K_use_tmp <- min(K_sens, P_full_tmp)
  k_svd_tmp <- max(5L, min(K_use_tmp, min(dim(ST_tr_full_tmp)) - 2L))
  if (k_svd_tmp < 5L) {
    warning("Skipping this mesh density: not enough dimension for SVD basis.")
    next
  }

  if (!requireNamespace("rsvd", quietly = TRUE)) install.packages("rsvd")
  svd_tmp <- rsvd::rsvd(ST_tr_full_tmp, k = k_svd_tmp)

  U_tmp <- svd_tmp$u[, 1:k_svd_tmp, drop=FALSE]
  D_tmp <- diag(svd_tmp$d[1:k_svd_tmp])
  V_tmp <- svd_tmp$v[, 1:k_svd_tmp, drop=FALSE]
  ST_tr_red_tmp <- U_tmp %*% D_tmp

  MQM_space_tmp <- as.matrix(t(mBase_tmp) %*% Q_tmp %*% mBase_tmp)
  MQM_space_tmp <- 0.5*(MQM_space_tmp + t(MQM_space_tmp)) + diag(1e-6, nrow(MQM_space_tmp))
  MQM_full_tmp  <- as.matrix(kronecker(MQM_space_tmp, Diagonal(num_t_basis)))
  MQM_full_tmp  <- 0.5*(MQM_full_tmp + t(MQM_full_tmp)) + diag(1e-6, nrow(MQM_full_tmp))

  r_sens <- min(as.integer(r_sens_target), k_svd_tmp)
  r_sens <- max(1L, r_sens)

  V_opt_tmp <- V_tmp[, 1:r_sens, drop=FALSE]
  MQM_red_tmp <- as.matrix(t(V_opt_tmp) %*% MQM_full_tmp %*% V_opt_tmp)
  MQM_red_tmp <- 0.5*(MQM_red_tmp + t(MQM_red_tmp)) + diag(1e-6, nrow(MQM_red_tmp))

  M_tr_tmp <- as.matrix(ST_tr_red_tmp[, 1:r_sens, drop=FALSE])

  build_time_sec <- as.numeric((proc.time() - t0_build)[3])

  out_sens <- run_picar_mcmc_one_rank(
    y_train       = as.numeric(response_train),
    X_matrix      = as.matrix(covariates_train),
    M_matrix      = M_tr_tmp,
    MQM_reduced   = MQM_red_tmp,
    niter         = niter_sens,
    burn_prop     = burn_sens,
    monitor_delta = FALSE
  )

  mesh_sens_res <- rbind(mesh_sens_res, data.frame(
    max_edge0          = me,
    N_mesh_vertices    = Nmesh_tmp,
    K_used             = k_svd_tmp,
    rank_used          = r_sens,
    basis_build_sec    = build_time_sec,
    mcmc_runtime_sec   = out_sens$runtime_sec,
    mcmc_time_per_iter = out_sens$runtime_sec / niter_sens,
    compile_time_sec   = out_sens$compile_time_sec,
    stringsAsFactors = FALSE
  ))
}

write.csv(mesh_sens_res, csv_mesh_sens, row.names=FALSE)

# Build the NEW PDF
pdf(pdf_new_analysis, width=11, height=8.5, onefile=TRUE)

pdf_text_page_12(
  "NEW ANALYSIS — SVD–PICAR (Gaussian, REAL DATA): Oversmoothing (maps+roughness) + Mesh-density decoupling",
  c(
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Time slice: t* = ", round(t_star, 4)),
    "Part 1: Fit many ranks r for each K; compute u_hat(s,t*) and edge-roughness R(u).",
    "Part 1 Maps: One panel per K at that K's best MSPE rank (common color scale).",
    "Part 2: Vary mesh density; measure basis build time and MCMC time/iter (fixed K, fixed r).",
    paste0("Roughness CSV: ", csv_rank_sweep),
    paste0("Mesh sensitivity CSV: ", csv_mesh_sens)
  )
)

# Part 1: Roughness curves vs rank
pdf_text_page_12(
  "PART 1 — Oversmoothing vs fine-scale: rank-sweep roughness curves (REAL DATA)",
  c(
    "For each K and rank r, we refit the Gaussian SVD–PICAR model and reconstruct the posterior-mean u_hat(s,t*).",
    "Roughness is computed on mesh vertices using neighbor differences along mesh edges:",
    "   R(u) = mean_{(i~j) in edges} (u_i - u_j)^2",
    "We report absolute roughness R_hat and scale-free roughness sqrt(R_hat)/sd(u_hat).",
    "We also normalize within each K by the *largest rank fitted* (ratio-to-reference), to show relative smoothing."
  )
)

if (nrow(rank_sweep_res) > 0) {

  p_abs <- ggplot(rank_sweep_res, aes(x=rank, y=R_hat)) +
    geom_line() + geom_point(size=1.7) +
    facet_wrap(~K_map, scales="free_x") +
    theme_minimal() +
    labs(
      title="Absolute edge roughness R_hat vs rank r",
      x="Rank r",
      y="R_hat"
    )
  print(p_abs)

  p_sf <- ggplot(rank_sweep_res, aes(x=rank, y=sf_hat)) +
    geom_line() + geom_point(size=1.7) +
    facet_wrap(~K_map, scales="free_x") +
    theme_minimal() +
    labs(
      title="Scale-free roughness sqrt(R_hat)/sd(u_hat) vs rank r",
      x="Rank r",
      y="sqrt(R_hat)/sd(u_hat)"
    )
  print(p_sf)

  if (all(c("R_ratio_ref","sf_ratio_ref") %in% names(rank_sweep_res))) {

    rank_sweep_res$logRratio_ref <- log10(rank_sweep_res$R_ratio_ref)

    p1 <- ggplot(rank_sweep_res, aes(x=rank, y=logRratio_ref)) +
      geom_hline(yintercept=0, linetype=2) +
      geom_line() + geom_point(size=1.7) +
      facet_wrap(~K_map, scales="free_x") +
      theme_minimal() +
      labs(
        title="log10( R_hat / R_ref ) vs rank r  (0 at the within-K reference rank)",
        subtitle="Reference is the largest rank fitted for that K",
        x="Rank r",
        y="log10(R_hat / R_ref)"
      )
    print(p1)

    p2 <- ggplot(rank_sweep_res, aes(x=rank, y=sf_ratio_ref)) +
      geom_hline(yintercept=1, linetype=2) +
      geom_line() + geom_point(size=1.7) +
      facet_wrap(~K_map, scales="free_x") +
      theme_minimal() +
      labs(
        title="Scale-free roughness ratio (sf_hat / sf_ref) vs rank r  (1 at reference rank)",
        subtitle="Reference is the largest rank fitted for that K",
        x="Rank r",
        y="sf_hat / sf_ref"
      )
    print(p2)
  }
}

# Part 1B: Latent field panels at t* for EACH K at best MSPE rank
pdf_text_page_12(
  "PART 1 — Latent field maps at t* for each K (best MSPE rank per K)",
  c(
    "For each K, we refit the Gaussian SVD–PICAR model at r*(K) = best MSPE rank.",
    "We reconstruct u_hat(s,t*) on mesh vertices and show one panel per K.",
    "All panels share a common color scale for direct visual comparison across K."
  )
)

get_best_rank_for_K <- function(K_map, K_used) {
  r_best <- NA_integer_
  if (exists("best_ranksel_byK") && nrow(best_ranksel_byK) > 0) {
    rr <- best_ranksel_byK[best_ranksel_byK$Constructed_ST_Basis == K_map, , drop=FALSE]
    if (nrow(rr) > 0) r_best <- as.integer(rr$Best_Rank[1])
  }
  if (is.na(r_best) || length(r_best) == 0) r_best <- min(20L, as.integer(K_used))
  r_best <- max(1L, min(as.integer(r_best), as.integer(K_used)))
  r_best
}

latent_byK <- data.frame()

for (K_map in K_vals_sweep) {

  cat("\n[MAP PANELS] Fitting best-rank map for K =", K_map, "...\n")

  bsK <- basis_store[[as.character(K_map)]]
  if (is.null(bsK)) stop(paste0("basis_store missing K=", K_map))

  ST_train_red <- as.matrix(bsK$ST_train_reduced)
  V_reduced    <- as.matrix(bsK$V_reduced)
  MQM_full     <- as.matrix(bsK$MQM_full)

  K_used <- ncol(ST_train_red)
  if (!is.null(bsK$K_used)) K_used <- as.integer(bsK$K_used)

  r_best <- get_best_rank_for_K(K_map, K_used)

  V_opt  <- V_reduced[, 1:r_best, drop=FALSE]
  MQM_red <- as.matrix(t(V_opt) %*% MQM_full %*% V_opt)
  MQM_red <- 0.5 * (MQM_red + t(MQM_red)) + diag(1e-6, nrow(MQM_red))

  M_train <- ST_train_red[, 1:r_best, drop=FALSE]

  out_fit <- run_picar_mcmc_one_rank(
    y_train       = as.numeric(response_train),
    X_matrix      = as.matrix(covariates_train),
    M_matrix      = M_train,
    MQM_reduced   = MQM_red,
    niter         = niter_rank_sweep,
    burn_prop     = burn_rank_sweep,
    monitor_delta = TRUE
  )

  delta_cols <- grep("^delta\\[", colnames(out_fit$post), value=TRUE)
  if (length(delta_cols) < 1) stop("No delta columns found for map fit (monitor_delta=TRUE expected).")

  delta_bar <- colMeans(out_fit$post[, delta_cols, drop=FALSE])
  delta_bar <- delta_bar[1:r_best]

  u_hat_mesh <- compute_u_mesh_at_tstar(
    mesh_obj    = mesh,
    mBase_full  = mBase_full,
    tBase_train = tBase_train,
    V_reduced   = V_reduced[, 1:K_used, drop=FALSE],
    rnk         = r_best,
    delta_bar   = delta_bar,
    t_star      = t_star
  )

  latent_byK <- rbind(
    latent_byK,
    data.frame(
      lon    = mesh$loc[,1],
      lat    = mesh$loc[,2],
      u_hat  = u_hat_mesh,
      K_map  = K_map,
      r_best = r_best,
      stringsAsFactors = FALSE
    )
  )
}

latent_byK$K_label <- paste0("K = ", latent_byK$K_map, "   (r* = ", latent_byK$r_best, ")")
latent_byK$K_label <- factor(latent_byK$K_label,
                             levels = unique(latent_byK$K_label[order(latent_byK$K_map)]))

lims_u <- range(latent_byK$u_hat, finite=TRUE)

p_panels <- ggplot(latent_byK, aes(x=lon, y=lat, color=u_hat)) +
  geom_point(size=0.55) +
  coord_equal() +
  facet_wrap(~K_label, ncol=3) +
  theme_minimal() +
  labs(
    title    = "Gaussian SVD–PICAR (REAL DATA): posterior-mean latent field u_hat(s,t*) at best MSPE rank per K",
    subtitle = paste0("t* = ", sprintf("%.4f", t_star), "   |   one panel per K (common color scale)"),
    x = "Longitude", y = "Latitude", color = "u_hat"
  ) +
  scale_color_continuous(limits = lims_u)

print(p_panels)

# Part 2: Mesh-density decoupling curves
pdf_text_page_12(
  "PART 2 — Mesh-density decoupling (build time vs mesh density; MCMC time/iter vs mesh density)",
  c(
    "We vary mesh density via max.edge0 (larger max.edge0 = coarser mesh).",
    "For each mesh we record:",
    "  (i) one-time basis build time (mesh + A + Moran + ST + SVD + MQM reduction)",
    "  (ii) MCMC sampling time per iteration at fixed K and fixed rank",
    paste0("Sensitivity settings: K_sens=",K_sens,", rank_sens=",r_sens_target,", niter_sens=",niter_sens," (pilot).")
  )
)

if (nrow(mesh_sens_res) > 0) {

  p_build <- ggplot(mesh_sens_res, aes(x=N_mesh_vertices, y=basis_build_sec)) +
    geom_line() + geom_point(size=2) +
    theme_minimal() +
    labs(title="One-time basis build time vs mesh density",
         x="Number of mesh vertices", y="Basis build time (sec)")

  p_mcmc <- ggplot(mesh_sens_res, aes(x=N_mesh_vertices, y=mcmc_time_per_iter)) +
    geom_line() + geom_point(size=2) +
    theme_minimal() +
    labs(title="MCMC time per iteration vs mesh density (fixed K, fixed rank)",
         x="Number of mesh vertices", y="Time per iteration (sec)")

  print(p_build)
  print(p_mcmc)

  tab_show <- mesh_sens_res
  tab_show$basis_build_sec    <- round(tab_show$basis_build_sec, 2)
  tab_show$mcmc_runtime_sec   <- round(tab_show$mcmc_runtime_sec, 2)
  tab_show$mcmc_time_per_iter <- signif(tab_show$mcmc_time_per_iter, 4)
  tab_show$compile_time_sec   <- round(tab_show$compile_time_sec, 2)

  tg <- tableGrob(tab_show, rows=NULL)
  grid.newpage(); grid.draw(tg)
}

pdf_text_page_12(
  "ANALYSIS OUTPUTS SAVED",
  c(
    paste0("PDF: ", pdf_new_analysis),
    paste0("Rank-sweep roughness CSV: ", csv_rank_sweep),
    paste0("Mesh-density sensitivity CSV: ", csv_mesh_sens)
  )
)

dev.off()

cat("\n[NEW ANALYSIS — REAL DATA] Saved:\n",
    "  PDF:", pdf_new_analysis, "\n",
    "  CSV (rank sweep):", csv_rank_sweep, "\n",
    "  CSV (mesh sensitivity):", csv_mesh_sens, "\n")

# END Section 12
