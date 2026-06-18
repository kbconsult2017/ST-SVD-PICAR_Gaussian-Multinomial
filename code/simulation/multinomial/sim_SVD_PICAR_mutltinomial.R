# Multinomial simulation study: ST SVD-PICAR

rm(list = ls())
options(stringsAsFactors = FALSE)

# 0) USER CONTROLS

# Multinomial categories
Kcat  <- 4
Km1   <- Kcat - 1

# Kbasis values
Kbasis_vals <- c(200, 300, 400, 500, 600)

# Full basis sizes (must satisfy p_full*q_full >= max(Kbasis_vals))
p_full <- 30
q_full <- 20
stopifnot(p_full * q_full >= max(Kbasis_vals))

# Rank selection grid
USE_RANK_GRID <- TRUE
make_rank_grid <- function(Kbasis) {
  if (!USE_RANK_GRID) return(5:Kbasis)

  v1 <- 5:25
  v2 <- if (Kbasis >= 30)  seq(30, min(200, Kbasis), by=10) else numeric(0)
  v3 <- if (Kbasis >= 225) seq(225, Kbasis, by=25)          else numeric(0)

  v <- sort(unique(c(v1, v2, v3, Kbasis)))
  v[v >= 5 & v <= Kbasis]
}

# Run 5-fold CV only for the best Kbasis
RUN_CV_BEST_KBASIS <- TRUE
Kfold <- 5

# Run NIMBLE MCMC for each Kbasis
RUN_MCMC_EACH_KBASIS <- TRUE

# MCMC length
NITER_MCMC <- 100000
BURN_FRAC  <- 0.50

# Output folder
out_dir <- "C:/Users/Admin/OneDrive - University of Nebraska Medical Center/Desktop/Svd-Multinomial"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "PerKbasis"), showWarnings = FALSE, recursive = TRUE)

# PDF filenames
pdf_mesh_file         <- file.path(out_dir, "01_mesh_and_locations.pdf")
pdf_rank_all_file     <- file.path(out_dir, "02_rank_selection_all_Kbasis.pdf")
pdf_coef_all_file     <- file.path(out_dir, "03_coeff_vs_rank_all_Kbasis.pdf")
pdf_bestK_file        <- file.path(out_dir, "04_best_metrics_vs_Kbasis.pdf")
pdf_cv_bestK_file     <- file.path(out_dir, "05_cv_best_Kbasis.pdf")

# PDFs for ESS analysis
pdf_ess_all_file       <- file.path(out_dir, "06_ess_all_Kbasis.pdf")
pdf_runtime_ess_file   <- file.path(out_dir, "07_bestRank_runtime_vs_ESSperSec.pdf")

# CSVs for ESS analysis
csv_ess_all_file       <- file.path(out_dir, "ESS_AllParams_AllKbasis.csv")
csv_ess_sum_file       <- file.path(out_dir, "ESS_Summary_ByKbasis.csv")
csv_ess_group_file     <- file.path(out_dir, "ESS_Summary_ByGroup_AllKbasis.csv")
csv_best_metrics_file  <- file.path(out_dir, "BestMetrics_VGAM_plus_MCMC.csv")

# 1) PACKAGES
required_packages <- c(
  "fields","mvtnorm","INLA","ggplot2","sf","Matrix","irlba","splines",
  "reshape2","dplyr","tidyr","nimble","coda","VGAM","rsvd"
)
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) if (!pkg %in% installed_packages) install.packages(pkg)

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
library(VGAM)
library(rsvd)

# 2) MULTI-COVARIATE SPECIFICATION
cov_spec <- list(
  gws_MEAN = list(dist="norm", mean=5,  sd=2, beta=c( 0.40, -0.10,  0.20)),
  rtz_MEAN = list(dist="norm", mean=10, sd=3, beta=c(-0.30,  0.25, -0.15))
)

cov_names <- names(cov_spec)
pX <- length(cov_names)

alpha_true <- c(0.30, -0.20, 0.10)
stopifnot(length(alpha_true) == Km1)

beta_true <- sapply(cov_names, function(nm) cov_spec[[nm]]$beta)
beta_true <- as.matrix(beta_true)
stopifnot(nrow(beta_true) == Km1, ncol(beta_true) == pX)

# 3) HELPER FUNCTIONS

multinom_nll <- function(Y, P, eps = 1e-12) {
  P <- pmax(P, eps)
  -sum(Y * log(P))
}

make_pi_baseline <- function(eta_mat) {
  exp_eta <- exp(eta_mat)
  den <- 1 + rowSums(exp_eta)
  cbind(exp_eta / den, 1 / den)
}

make_vglm_formula <- function(Kcat, cov_names, st_rank) {
  resp <- paste0("cbind(", paste0("Y", 1:Kcat, collapse = ","), ")")
  rhs  <- paste(c(cov_names, paste0("STBase_", 1:st_rank)), collapse = " + ")
  as.formula(paste(resp, "~", rhs))
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

build_moran_basis <- function(W, num_eigen) {
  N <- nrow(W)
  OrthSpace <- Diagonal(N) - (1/N) * (Matrix(1, N, 1) %*% Matrix(1, 1, N))
  MoransOperator <- OrthSpace %*% (W %*% OrthSpace)
  eig <- irlba(MoransOperator, nv = num_eigen, nu = num_eigen)
  eig$v[, 1:num_eigen, drop = FALSE]
}

make_st_scores <- function(spatial_scores, temporal_scores) {
  n <- nrow(spatial_scores); p <- ncol(spatial_scores); q <- ncol(temporal_scores)
  ST <- matrix(NA_real_, nrow = n, ncol = p*q)
  for (i in 1:n) ST[i, ] <- as.vector(tcrossprod(temporal_scores[i, ], spatial_scores[i, ]))
  ST
}

center_basis_train_test <- function(B_train, B_test) {
  mu <- colMeans(B_train)
  list(
    B_train = sweep(B_train, 2, mu, "-"),
    B_test  = sweep(B_test,  2, mu, "-"),
    center  = mu
  )
}

# Coef extraction from VGAM with multiple covariates
extract_vgam_coef_ci <- function(coef_table, term, k) {
  rn <- rownames(coef_table)
  if (term == "(Intercept)") {
    pat <- paste0("^\\(Intercept\\):\\s*", k, "$")
  } else {
    term_esc <- gsub("([\\+\\*\\?\\^\\$\\(\\)\\[\\]\\{\\}\\.\\|\\\\])", "\\\\\\1", term)
    pat <- paste0("^", term_esc, ":\\s*", k, "$")
  }
  hit <- grep(pat, rn)
  if (length(hit) == 0) return(c(est=NA_real_, lcl=NA_real_, ucl=NA_real_))
  est <- coef_table[hit[1], "Estimate"]
  se  <- coef_table[hit[1], "Std. Error"]
  z   <- qnorm(0.975)
  c(est=est, lcl=est - z*se, ucl=est + z*se)
}

# 4) BUILD FULL ST BASIS ONCE (p_full * q_full columns) + MQM_full
build_full_st_once <- function(AMat_train, AMat_test, Q_lap, mBase, tBase_train, tBase_test,
                               p_full, q_full) {

  mBase_sub <- mBase[, 1:p_full, drop = FALSE]
  tTr_sub   <- tBase_train[, 1:q_full, drop = FALSE]
  tTe_sub   <- tBase_test[,  1:q_full, drop = FALSE]

  spatial_train <- as.matrix(AMat_train %*% mBase_sub)
  spatial_test  <- as.matrix(AMat_test  %*% mBase_sub)

  ST_train_full <- make_st_scores(spatial_train, tTr_sub)
  ST_test_full  <- make_st_scores(spatial_test,  tTe_sub)

  MQM_s <- t(mBase_sub) %*% Q_lap %*% mBase_sub
  MQM_s <- 0.5*(MQM_s + t(MQM_s)) + diag(1e-6, nrow(MQM_s))
  MQM_full <- kronecker(diag(1, q_full), MQM_s)
  MQM_full <- 0.5*(MQM_full + t(MQM_full)) + diag(1e-6, nrow(MQM_full))

  list(
    ST_train_full = ST_train_full,
    ST_test_full  = ST_test_full,
    MQM_full      = MQM_full,
    mBase_sub     = mBase_sub,
    tBase_train_sub = tTr_sub,
    tBase_test_sub  = tTe_sub
  )
}

# SAFETY HELPERS

safe_max_rank <- function(n_train, Kcat, pX, frac = 0.90, buffer = 10) {
  Km1 <- Kcat - 1
  rhs  <- floor(frac * (n_train - buffer))
  rmax <- floor(rhs / Km1) - (1 + pX)
  max(5, rmax)
}

clip_eta <- function(eta, lim = 30) pmax(pmin(eta, lim), -lim)

# VGAM RANK SELECTION
run_rank_selection_vgam <- function(Kcat, cov_names,
                                    Y_train, Y_test,
                                    cov_train_df, cov_test_df,
                                    nsize_test,
                                    B_train_full, B_test_full,
                                    rank_sequence) {

  res <- data.frame()

  for (r in rank_sequence) {
    t0 <- Sys.time()

    Btr <- B_train_full[, 1:r, drop=FALSE]
    Bte <- B_test_full[,  1:r, drop=FALSE]
    colnames(Btr) <- paste0("STBase_", 1:r)
    colnames(Bte) <- paste0("STBase_", 1:r)

    df_train <- data.frame(Y_train, cov_train_df, Btr, check.names = FALSE)
    df_test  <- data.frame(Y_test,  cov_test_df,  Bte, check.names = FALSE)

    for (k in 1:Kcat) {
      df_train[[paste0("Y",k)]] <- as.integer(df_train[[paste0("Y",k)]])
      df_test[[paste0("Y",k)]]  <- as.integer(df_test[[paste0("Y",k)]])
    }

    form_r <- make_vglm_formula(Kcat, cov_names, r)

    fit <- try(
      vglm(form_r, family = multinomial(refLevel = Kcat), data = df_train,
           control = vglm.control(maxit = 100)),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      cat(sprintf("rank %d | VGAM failed -> skipped\n", r))
      next
    }

    P_test <- try(predict(fit, newdata = df_test, type = "response"), silent = TRUE)
    if (inherits(P_test, "try-error") || anyNA(P_test) || any(!is.finite(P_test))) {
      cat(sprintf("rank %d | predict produced NA/Inf -> skipped\n", r))
      next
    }

    rs <- rowSums(P_test)
    if (any(!is.finite(rs)) || max(abs(rs - 1)) > 1e-6) {
      cat(sprintf("rank %d | invalid probability rows -> skipped\n", r))
      next
    }

    TestNLL <- multinom_nll(Y_test, P_test)

    pred_counts <- sweep(P_test, 1, nsize_test, `*`)
    MSPE_count  <- mean((as.matrix(Y_test) - pred_counts)^2)

    Y_prop <- sweep(as.matrix(Y_test), 1, rowSums(as.matrix(Y_test)), `/`)
    MSPE_prop <- mean((Y_prop - P_test)^2)

    time_sec <- as.numeric(difftime(Sys.time(), t0, units="secs"))

    ct <- as.data.frame(coef(summary(fit)))
    row <- data.frame(ST_Rank=r, TestNLL=TestNLL, MSPE_count=MSPE_count, MSPE_prop=MSPE_prop, TimeSec=time_sec)

    for (kk in 1:(Kcat-1)) {
      a <- extract_vgam_coef_ci(ct, "(Intercept)", kk)
      row[[paste0("alpha",kk)]]        <- a["est"]
      row[[paste0("alpha",kk,"_LCL")]] <- a["lcl"]
      row[[paste0("alpha",kk,"_UCL")]] <- a["ucl"]

      for (nm in cov_names) {
        b <- extract_vgam_coef_ci(ct, nm, kk)
        cname <- paste0("b_", nm, "_", kk)
        row[[cname]]                <- b["est"]
        row[[paste0(cname,"_LCL")]] <- b["lcl"]
        row[[paste0(cname,"_UCL")]] <- b["ucl"]
      }
    }

    res <- rbind(res, row)
    cat(sprintf("rank %d | TestNLL %.2f | MSPE(count) %.4f | time %.3f s\n",
                r, TestNLL, MSPE_count, time_sec))
  }

  if (nrow(res) == 0) stop("All ranks failed. Reduce rank grid or increase stabilization.")
  best_rank <- res$ST_Rank[which.min(res$MSPE_count)]
  list(results=res, best_rank=best_rank)
}

# K-FOLD CV
run_kfold_cv_vgam <- function(Kcat, cov_names, Y, cov_df, B_full, Kfold=5, rank_sequence, seed=123) {
  set.seed(seed)
  n <- nrow(Y)
  fold_id <- sample(rep(1:Kfold, length.out=n))
  out <- data.frame()

  for (r in rank_sequence) {
    t0 <- Sys.time()
    B_r <- B_full[, 1:r, drop=FALSE]
    colnames(B_r) <- paste0("STBase_", 1:r)

    df_all <- data.frame(Y, cov_df, B_r, check.names = FALSE)
    for (k in 1:Kcat) df_all[[paste0("Y",k)]] <- as.integer(df_all[[paste0("Y",k)]])

    form_r <- make_vglm_formula(Kcat, cov_names, r)

    fold_nll <- rep(NA_real_, Kfold)
    for (f in 1:Kfold) {
      idx_te <- which(fold_id == f)
      idx_tr <- which(fold_id != f)

      fit <- try(
        vglm(form_r, family = multinomial(refLevel = Kcat),
             data = df_all[idx_tr, , drop=FALSE],
             control = vglm.control(maxit = 100)),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) next

      P_te <- try(predict(fit, newdata = df_all[idx_te, , drop=FALSE], type="response"), silent=TRUE)
      if (inherits(P_te, "try-error") || anyNA(P_te) || any(!is.finite(P_te))) next

      Y_te <- as.matrix(df_all[idx_te, paste0("Y",1:Kcat), drop=FALSE])
      fold_nll[f] <- multinom_nll(Y_te, P_te)
    }

    ok <- is.finite(fold_nll)
    if (sum(ok) < 2) next

    dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    out <- rbind(out, data.frame(
      ST_Rank=r,
      CV_NLL=mean(fold_nll[ok]),
      CV_NLL_SD=sd(fold_nll[ok]),
      Folds_OK=sum(ok),
      Folds_Failed=Kfold-sum(ok),
      TimeSec=dt
    ))
    cat(sprintf("CV rank %d | mean NLL %.2f | ok %d/%d\n", r, mean(fold_nll[ok]), sum(ok), Kfold))
  }

  if (nrow(out) == 0) stop("All CV ranks failed. Reduce rank grid.")
  best_rank_cv <- out$ST_Rank[which.min(out$CV_NLL)]
  list(cv_results=out, best_rank_cv=best_rank_cv)
}

# 5) GIVEN FULL ST, BUILD Kbasis REDUCED BASIS VIA rSVD
build_Kbasis <- function(ST_train_full, ST_test_full, MQM_full, Kbasis, seed = 1234) {
  set.seed(seed)
  k_svd <- min(Kbasis, nrow(ST_train_full), ncol(ST_train_full))
  sv <- rsvd::rsvd(ST_train_full, k = k_svd)

  U <- sv$u
  d <- sv$d
  V <- sv$v

  B_train <- U %*% diag(d, nrow=length(d), ncol=length(d))
  B_test  <- ST_test_full %*% V

  colnames(B_train) <- paste0("STBase_", 1:ncol(B_train))
  colnames(B_test)  <- paste0("STBase_", 1:ncol(B_test))

  MQM_reduced <- t(V) %*% MQM_full %*% V
  MQM_reduced <- 0.5*(MQM_reduced + t(MQM_reduced)) + diag(1e-6, nrow(MQM_reduced))

  list(B_train = B_train, B_test = B_test, V = V, MQM_reduced = MQM_reduced)
}

# 8) NIMBLE FIT FUNCTION
make_tau_gibbs_sampler <- function() {
  nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      a0  <- control$a0
      b0  <- control$b0
      MQM <- as.matrix(control$MQM)
      r   <- as.integer(dim(MQM)[1])
      calcNodes <- model$getDependencies(target)
    },
    run = function() {
      deltaVec <- model[["delta"]]
      quad <- 0.0
      for (i in 1:r) {
        tmpi <- 0.0
        for (j in 1:r) tmpi <- tmpi + MQM[i, j] * deltaVec[j]
        quad <- quad + deltaVec[i] * tmpi
      }
      shape_post <- a0 + 0.5 * r
      rate_post  <- b0 + 0.5 * quad
      model[[target]] <<- rgamma(1, shape = shape_post, rate = rate_post)
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, nodes = target, row = 1)
    },
    methods = list(reset = function() {})
  )
}

run_nimble_multinom <- function(Kcat, cov_names,
                                Y_train_mat, X_train_std, B_train_rank, MQM_rank,
                                nsize_train,
                                niter=60000, burn_frac=0.7,
                                seed=1) {

  set.seed(seed)

# TOTAL wall-clock timer
  pt_total <- proc.time()

  Km1 <- Kcat - 1
  n   <- nrow(Y_train_mat)
  pX  <- ncol(X_train_std)
  r   <- ncol(B_train_rank)

  data_list <- list(Y = Y_train_mat)
  consts <- list(
    n=n, K=Kcat, Km1=Km1, pX=pX, r=r,
    X=X_train_std, B=B_train_rank, MQM=MQM_rank,
    zero=rep(0, r),
    nsize=as.integer(nsize_train)
  )
  inits <- list(alpha=rep(0, Km1), beta=matrix(0, nrow=Km1, ncol=pX), delta=rep(0, r), tau=1)

  code <- nimbleCode({
    for (j in 1:n) {
      W[j] <- inprod(B[j, 1:r], delta[1:r])
      for (k in 1:Km1) {
        eta[j, k] <- alpha[k] + inprod(X[j, 1:pX], beta[k, 1:pX]) + W[j]
        exp_eta[j, k] <- exp(eta[j, k])
      }
      sum_exp[j] <- sum(exp_eta[j, 1:Km1])
      for (k in 1:Km1) {
        prob[j, k] <- exp_eta[j, k] / (1 + sum_exp[j])
      }
      prob[j, K] <- 1 / (1 + sum_exp[j])
      Y[j, 1:K] ~ dmulti(prob[j, 1:K], size = nsize[j])
    }

    precMat[1:r, 1:r] <- tau * MQM[1:r, 1:r]
    delta[1:r] ~ dmnorm(mean = zero[1:r], prec = precMat[1:r, 1:r])

    for (k in 1:Km1) {
      alpha[k] ~ dnorm(0, var=100)
      for (l in 1:pX) beta[k, l] ~ dnorm(0, var=100)
    }
    tau ~ dgamma(shape=2, rate=1)
  })

# Build + compile model
  pt0 <- proc.time()
  Rmodel <- nimbleModel(code=code, data=data_list, constants=consts, inits=inits)
  t_build_model_sec <- as.numeric((proc.time() - pt0)[3])

  pt1 <- proc.time()
  Cmodel <- compileNimble(Rmodel)
  compile_model_sec <- as.numeric((proc.time() - pt1)[3])

# MCMC configuration
  tauGibbsSampler <- make_tau_gibbs_sampler()
  conf <- configureMCMC(Rmodel, print=FALSE)

  conf$removeSamplers(paste0("alpha[", 1:Km1, "]"))
  for (k in 1:Km1) conf$removeSamplers(paste0("beta[", k, ",", 1:pX, "]"))
  conf$removeSamplers(paste0("delta[", 1:r, "]"))
  conf$removeSamplers("tau")

  for (k in 1:Km1) {
    conf$addSampler(target=paste0("alpha[",k,"]"), type="RW",
                    control=list(adaptive=TRUE, scale=0.5))
    conf$addSampler(target=paste0("beta[",k,",1:",pX,"]"), type="RW_block",
                    control=list(adaptive=TRUE, scale=0.2))
  }
  conf$addSampler(target=paste0("delta[1:",r,"]"), type="RW_block",
                  control=list(adaptive=TRUE, scale=0.1))

  conf$addSampler(target="tau", type=tauGibbsSampler,
                  control=list(a0=2, b0=1, MQM=as.matrix(MQM_rank)))

  conf$setMonitors(c("alpha","beta","delta","tau"))

  mcmc <- buildMCMC(conf)

# Compile MCMC (timed) — IMPORTANT: project = Cmodel
  pt2 <- proc.time()
  Cmcmc <- compileNimble(mcmc, project=Cmodel)
  compile_mcmc_sec <- as.numeric((proc.time() - pt2)[3])

# Run MCMC
  pt3 <- proc.time()
  Cmcmc$run(niter)

# Sample extraction time (timed) — often nontrivial for large niter
  pt4 <- proc.time()
  samples <- as.matrix(Cmcmc$mvSamples)
  extract_sec <- as.numeric((proc.time() - pt4)[3])

  burn <- floor(burn_frac * nrow(samples))
  post <- samples[(burn+1):nrow(samples), , drop=FALSE]

# TOTAL wall-clock runtime (ONLY): build + compile + run + extract
  runtime_sec <- as.numeric((proc.time() - pt_total)[3])
  time_per_iter_sec <- runtime_sec / max(niter, 1)

  list(
    samples = samples,
    post    = post,

    build_model_sec   = t_build_model_sec,
    compile_model_sec = compile_model_sec,
    compile_mcmc_sec  = compile_mcmc_sec,
    extract_sec       = extract_sec,

    runtime_sec       = runtime_sec,
    time_per_iter_sec = time_per_iter_sec
  )
}

# 9) SIMULATE DATA (GEOID, TIME, COVARS, NSIZE) + BUILD MESH + MORAN +
overall_start <- Sys.time()
set.seed(123)

n_total <- 2500
n_unique_geoids <- 200
unique_geoids <- paste0("GID_", 1:n_unique_geoids)
GEOID_sim <- sample(unique_geoids, size = n_total, replace = TRUE)

time_points <- seq(from = 2015, to = 2020, by = 1/12)
year_mo_sim <- sample(time_points, n_total, replace = TRUE)

cov_mat <- matrix(NA_real_, nrow=n_total, ncol=pX)
colnames(cov_mat) <- cov_names
for (j in 1:pX) {
  sp <- cov_spec[[cov_names[j]]]
  if (sp$dist != "norm") stop("Only normal covariates implemented; extend if needed.")
  cov_mat[, j] <- rnorm(n_total, mean=sp$mean, sd=sp$sd)
}

nsize_sim <- sample(30:80, n_total, replace = TRUE)

lat_range <- c(32.5343, 42.0095)
lon_range <- c(-124.4096, -113.0460)

coords_lookup <- data.frame(
  GEOID = unique_geoids,
  lon   = runif(length(unique_geoids), lon_range[1], lon_range[2]),
  lat   = runif(length(unique_geoids), lat_range[1], lat_range[2])
)

data0 <- data.frame(GEOID=GEOID_sim, year_mo=as.numeric(year_mo_sim), nsize=as.integer(nsize_sim))
for (j in 1:pX) data0[[cov_names[j]]] <- cov_mat[, j]
data0 <- merge(data0, coords_lookup, by="GEOID", all.x=TRUE, sort=FALSE)

set.seed(123)
train_idx <- sample(1:nrow(data0), size = floor(0.7*nrow(data0)))
train_data <- data0[train_idx, ]
test_data  <- data0[-train_idx, ]

train_coords <- as.matrix(train_data[, c("lon","lat")])
max.edge0 <- 1.5
mesh <- inla.mesh.2d(
  loc      = train_coords,
  max.edge = c(1,2) * max.edge0,
  cutoff   = max.edge0/5,
  offset   = c(max.edge0, 6.0)
)

pdf(pdf_mesh_file, width=8, height=6)
plot(mesh, main="INLA Mesh + Training Locations")
points(train_coords, col="blue", pch=16, cex=0.5)
dev.off()

AMat_train <- inla.spde.make.A(mesh, loc = as.matrix(train_data[, c("lon","lat")]))
AMat_test  <- inla.spde.make.A(mesh, loc = as.matrix(test_data[,  c("lon","lat")]))

lap <- mesh_to_laplacian(mesh)
W_adj <- lap$W
Q_lap <- as.matrix(lap$Q) + diag(1e-6, lap$Nnew)

mBase <- build_moran_basis(W_adj, num_eigen = p_full)

tBase_train <- bs(train_data$year_mo, df=q_full, intercept=FALSE)
tBase_test  <- bs(test_data$year_mo,  df=q_full, intercept=FALSE)

# 10) BUILD FULL ST ONCE + SIMULATE Y USING Kbasis_truth = max
full_ST <- build_full_st_once(AMat_train, AMat_test, Q_lap, mBase, tBase_train, tBase_test,
                              p_full=p_full, q_full=q_full)

Kbasis_truth <- max(Kbasis_vals)

kb_truth <- build_Kbasis(full_ST$ST_train_full, full_ST$ST_test_full, full_ST$MQM_full,
                         Kbasis=Kbasis_truth, seed=1000 + Kbasis_truth)

Bc_truth <- center_basis_train_test(kb_truth$B_train, kb_truth$B_test)
B_train_truth <- Bc_truth$B_train
B_test_truth  <- Bc_truth$B_test

r_true <- ncol(B_train_truth)
tau_true <- 2
Qdelta <- tau_true * kb_truth$MQM_reduced + diag(1e-8, r_true)
Rchol <- chol(Qdelta)
delta_true <- backsolve(Rchol, rnorm(r_true))

W_true_train <- as.vector(B_train_truth %*% delta_true)
W_true_test  <- as.vector(B_test_truth  %*% delta_true)

eta_train <- matrix(0, nrow=nrow(train_data), ncol=Km1)
eta_test  <- matrix(0, nrow=nrow(test_data),  ncol=Km1)

Xtr <- as.matrix(train_data[, cov_names, drop=FALSE])
Xte <- as.matrix(test_data[,  cov_names, drop=FALSE])

for (k in 1:Km1) {
  eta_train[,k] <- alpha_true[k] + as.vector(Xtr %*% beta_true[k, ]) + W_true_train
  eta_test[,k]  <- alpha_true[k] + as.vector(Xte %*% beta_true[k, ]) + W_true_test
}

P_train <- make_pi_baseline(eta_train)
P_test  <- make_pi_baseline(eta_test)

Y_train <- matrix(0L, nrow=nrow(train_data), ncol=Kcat)
Y_test  <- matrix(0L, nrow=nrow(test_data),  ncol=Kcat)

for (j in 1:nrow(train_data)) Y_train[j,] <- as.vector(rmultinom(1, size=train_data$nsize[j], prob=P_train[j,]))
for (j in 1:nrow(test_data))  Y_test[j,]  <- as.vector(rmultinom(1, size=test_data$nsize[j],  prob=P_test[j,]))

colnames(Y_train) <- paste0("Y",1:Kcat)
colnames(Y_test)  <- paste0("Y",1:Kcat)
train_data[, paste0("Y",1:Kcat)] <- Y_train
test_data[,  paste0("Y",1:Kcat)] <- Y_test

# 11) LOOP OVER Kbasis: BUILD BASIS -> VGAM RANK SELECTION -> SAVE
rank_all <- data.frame()
best_by_Kbasis <- data.frame()

# 11b) GLOBAL ESS COLLECTORS
ess_all <- data.frame()
ess_summary_by_Kbasis <- data.frame()

for (Kbasis in Kbasis_vals) {

  cat("\n====================================================\n")
  cat("Kbasis =", Kbasis, "(constructed SVD basis size)\n")
  cat("====================================================\n")

  per_dir <- file.path(out_dir, "PerKbasis", paste0("Kbasis_", Kbasis))
  dir.create(per_dir, showWarnings = FALSE, recursive = TRUE)

  kb <- build_Kbasis(full_ST$ST_train_full, full_ST$ST_test_full, full_ST$MQM_full,
                     Kbasis=Kbasis, seed=1000 + Kbasis)

  Bc <- center_basis_train_test(kb$B_train, kb$B_test)
  B_train_full <- Bc$B_train
  B_test_full  <- Bc$B_test

  cov_train_df <- train_data[, cov_names, drop=FALSE]
  cov_test_df  <- test_data[,  cov_names, drop=FALSE]

# CONSISTENT rank cap
  r_cap <- min(
    ncol(B_train_full),
    safe_max_rank(n_train = nrow(train_data), Kcat = Kcat, pX = pX, frac = 0.90, buffer = 10)
  )

  rank_sequence <- make_rank_grid(r_cap)

  rs <- run_rank_selection_vgam(
    Kcat=Kcat, cov_names=cov_names,
    Y_train=as.data.frame(train_data[, paste0("Y",1:Kcat)]),
    Y_test =as.data.frame(test_data[,  paste0("Y",1:Kcat)]),
    cov_train_df=cov_train_df,
    cov_test_df=cov_test_df,
    nsize_test=test_data$nsize,
    B_train_full=B_train_full,
    B_test_full=B_test_full,
    rank_sequence=rank_sequence
  )

  tmp <- rs$results
  tmp$Kbasis <- Kbasis
  rank_all <- rbind(rank_all, tmp)

  best_rank <- rs$best_rank
  best_row <- tmp[tmp$ST_Rank == best_rank, , drop=FALSE]

  best_by_Kbasis <- rbind(best_by_Kbasis, data.frame(
    Kbasis=Kbasis,
    best_rank=best_rank,
    best_TestNLL=best_row$TestNLL,
    best_MSPE_count=best_row$MSPE_count,
    best_MSPE_prop=best_row$MSPE_prop,
    best_VGAM_TimeSec=best_row$TimeSec
  ))

  write.csv(tmp, file.path(per_dir, "RankSelection_AllRanks.csv"), row.names=FALSE)
  write.csv(best_row, file.path(per_dir, "RankSelection_BestRow.csv"), row.names=FALSE)

  pdf(file.path(per_dir, "RankSelection_Plots.pdf"), width=10, height=8)
  print(ggplot(tmp, aes(ST_Rank, TestNLL)) + geom_line() + geom_point() +
          theme_minimal() + ggtitle(paste("Test NLL vs Rank | Kbasis =",Kbasis)))
  print(ggplot(tmp, aes(ST_Rank, MSPE_count)) + geom_line() + geom_point() +
          theme_minimal() + ggtitle(paste("MSPE(count) vs Rank | Kbasis =",Kbasis)))
  print(ggplot(tmp, aes(ST_Rank, MSPE_prop)) + geom_line() + geom_point() +
          theme_minimal() + ggtitle(paste("MSPE(prop) vs Rank | Kbasis =",Kbasis)))
  print(ggplot(tmp, aes(ST_Rank, TimeSec)) + geom_line() + geom_point() +
          theme_minimal() + ggtitle(paste("Runtime vs Rank | Kbasis =",Kbasis)))
  dev.off()

# OPTIONAL: NIMBLE MCMC FOR THIS Kbasis AT ITS BEST RANK
  if (RUN_MCMC_EACH_KBASIS) {

    r_use <- best_rank
    B_train_rank <- B_train_full[, 1:r_use, drop=FALSE]
    B_test_rank  <- B_test_full[,  1:r_use, drop=FALSE]
    MQM_rank <- kb$MQM_reduced[1:r_use, 1:r_use, drop=FALSE]

    X_train_raw <- as.matrix(train_data[, cov_names, drop=FALSE])
    X_train_std <- scale(X_train_raw)
    muX <- attr(X_train_std, "scaled:center")
    sdX <- attr(X_train_std, "scaled:scale")

    X_test_raw <- as.matrix(test_data[, cov_names, drop=FALSE])
    X_test_std <- sweep(X_test_raw, 2, muX, "-")
    X_test_std <- sweep(X_test_std, 2, sdX, "/")

    alpha_truth_model <- alpha_true + as.vector(beta_true %*% muX)
    beta_truth_model  <- sweep(beta_true, 2, sdX, `*`)

    fit <- run_nimble_multinom(
      Kcat=Kcat, cov_names=cov_names,
      Y_train_mat=as.matrix(train_data[, paste0("Y",1:Kcat)]),
      X_train_std=X_train_std,
      B_train_rank=B_train_rank,
      MQM_rank=MQM_rank,
      nsize_train=train_data$nsize,
      niter=NITER_MCMC, burn_frac=BURN_FRAC,
      seed=500 + Kbasis + r_use
    )

    samples <- fit$samples
    post    <- fit$post

    runtime_sec   <- fit$runtime_sec
    time_per_iter <- fit$time_per_iter_sec

    save(samples, file=file.path(per_dir, "MCMC_Samples.RData"))

    cn <- colnames(post)
    summ <- function(x) c(mean=mean(x), lcl=unname(quantile(x,0.025)), ucl=unname(quantile(x,0.975)))

    rows <- list()
    for (k in 1:Km1) {
      a <- summ(post[, sprintf("alpha[%d]",k)])
      rows[[length(rows)+1]] <- data.frame(
        Param=sprintf("alpha[%d]",k), Truth=alpha_truth_model[k],
        Mean=a["mean"], LCL=a["lcl"], UCL=a["ucl"],
        Covered=(alpha_truth_model[k] >= a["lcl"] && alpha_truth_model[k] <= a["ucl"])
      )
      for (j in 1:pX) {
        nm <- grep(sprintf("^beta\\[%d,\\s*%d\\]$", k, j), cn, value=TRUE)[1]
        b <- summ(post[, nm])
        rows[[length(rows)+1]] <- data.frame(
          Param=sprintf("beta[%d,%s]", k, cov_names[j]), Truth=beta_truth_model[k,j],
          Mean=b["mean"], LCL=b["lcl"], UCL=b["ucl"],
          Covered=(beta_truth_model[k,j] >= b["lcl"] && beta_truth_model[k,j] <= b["ucl"])
        )
      }
    }
    check_tab <- do.call(rbind, rows)
    check_tab$Bias <- check_tab$Mean - check_tab$Truth
    write.csv(check_tab, file.path(per_dir, "Posterior_Coverage_Table.csv"), row.names=FALSE)

    ess_post <- coda::effectiveSize(coda::as.mcmc(post))
    rt_den   <- max(runtime_sec, 1e-8)

    ess_tab <- data.frame(
      Parameter   = names(ess_post),
      ESS         = as.numeric(ess_post),
      ESS_per_sec = as.numeric(ess_post) / rt_den
    )

    ess_tab$RuntimeSec     <- runtime_sec
    ess_tab$TimePerIterSec <- time_per_iter

    write.csv(ess_tab, file.path(per_dir, "ESS_and_ESSperSec.csv"), row.names=FALSE)

# ESS analysis for ALL parameters across ALL Kbasis

    ess_tab2 <- ess_tab
    ess_tab2$Kbasis <- Kbasis
    ess_tab2$Rank   <- r_use
    ess_tab2$RuntimeSec <- runtime_sec

    ess_tab2$Group <- dplyr::case_when(
      grepl("^alpha\\[", ess_tab2$Parameter) ~ "alpha",
      grepl("^beta\\[",  ess_tab2$Parameter) ~ "beta",
      grepl("^delta\\[", ess_tab2$Parameter) ~ "delta",
      ess_tab2$Parameter == "tau" ~ "tau",
      TRUE ~ "other"
    )

    ess_all <- rbind(ess_all, ess_tab2)

    ess_sum <- data.frame(
      Kbasis = Kbasis,
      Rank   = r_use,

      RuntimeSec     = runtime_sec,
      TimePerIterSec = time_per_iter,

      ESS_mean_all   = mean(ess_tab2$ESS, na.rm=TRUE),
      ESS_med_all    = median(ess_tab2$ESS, na.rm=TRUE),
      ESS_min_all    = min(ess_tab2$ESS, na.rm=TRUE),

      ESSps_mean_all = mean(ess_tab2$ESS_per_sec, na.rm=TRUE),
      ESSps_med_all  = median(ess_tab2$ESS_per_sec, na.rm=TRUE),
      ESSps_min_all  = min(ess_tab2$ESS_per_sec, na.rm=TRUE),

      ESSps_mean_fixed = mean(ess_tab2$ESS_per_sec[ess_tab2$Group %in% c("alpha","beta","tau")], na.rm=TRUE),
      ESSps_min_fixed  = min( ess_tab2$ESS_per_sec[ess_tab2$Group %in% c("alpha","beta","tau")], na.rm=TRUE),

      ESSps_mean_delta = mean(ess_tab2$ESS_per_sec[ess_tab2$Group == "delta"], na.rm=TRUE),
      ESSps_min_delta  = min( ess_tab2$ESS_per_sec[ess_tab2$Group == "delta"], na.rm=TRUE),

      n_params = nrow(ess_tab2)
    )
    ess_summary_by_Kbasis <- rbind(ess_summary_by_Kbasis, ess_sum)

    pdf(file.path(per_dir, "ESS_Plots.pdf"), width=11, height=8.5)

    print(
      ggplot(ess_tab2, aes(x = ESS)) +
        geom_histogram(bins = 40) +
        facet_wrap(~Group, scales = "free_y") +
        theme_minimal() +
        ggtitle(paste("ESS distribution by group | Kbasis", Kbasis, "| rank", r_use))
    )

    print(
      ggplot(ess_tab2, aes(x = ESS_per_sec)) +
        geom_histogram(bins = 40) +
        facet_wrap(~Group, scales = "free_y") +
        theme_minimal() +
        ggtitle(paste("ESS/sec distribution by group | Kbasis", Kbasis, "| rank", r_use))
    )

    print(
      ggplot(ess_tab2, aes(x = Group, y = ESS_per_sec)) +
        geom_boxplot() +
        theme_minimal() +
        ggtitle(paste("ESS/sec boxplots (ALL params) | Kbasis", Kbasis, "| rank", r_use))
    )

    dev.off()

    pdf(file.path(per_dir, "MCMC_Diagnostics.pdf"), width=10, height=8)
    par(mfrow=c(3,2))
    plot(samples[,"tau"], type="l", main=paste("Trace tau | Kbasis",Kbasis,"rank",r_use), xlab="iter", ylab="tau")
    for (k in 1:Km1) plot(samples[,sprintf("alpha[%d]",k)], type="l", main=paste("Trace alpha[",k,"]",sep=""), xlab="iter", ylab="")
    for (k in 1:min(Km1,2)) {
      for (j in 1:min(pX,2)) {
        nm <- grep(sprintf("^beta\\[%d,\\s*%d\\]$",k,j), colnames(samples), value=TRUE)[1]
        plot(samples[,nm], type="l", main=paste("Trace",nm), xlab="iter", ylab="")
      }
    }
    par(mfrow=c(3,2))
    acf(samples[,"tau"], main="ACF tau")
    for (k in 1:Km1) acf(samples[,sprintf("alpha[%d]",k)], main=paste("ACF alpha[",k,"]",sep=""))
    dev.off()

    thin_pp <- 10
    post_pp <- post[seq(1, nrow(post), by=thin_pp), , drop=FALSE]
    S <- nrow(post_pp)

    alpha_draws <- sapply(1:Km1, function(k) post_pp[, sprintf("alpha[%d]",k)])
    beta_draws <- array(NA_real_, dim=c(S, Km1, pX))
    for (k in 1:Km1) {
      for (j in 1:pX) {
        nm <- grep(sprintf("^beta\\[%d,\\s*%d\\]$",k,j), colnames(post_pp), value=TRUE)[1]
        beta_draws[,k,j] <- post_pp[, nm]
      }
    }
    delta_draws <- sapply(1:r_use, function(j) post_pp[, sprintf("delta[%d]",j)])
    W_test_draws <- B_test_rank %*% t(delta_draws)

    n_test <- nrow(X_test_std)
    P_bar <- matrix(0, nrow=n_test, ncol=Kcat)
    for (s in 1:S) {
      eta_s <- matrix(0, nrow=n_test, ncol=Km1)
      for (k in 1:Km1) {
        eta_s[,k] <- alpha_draws[s,k] + as.vector(X_test_std %*% beta_draws[s,k,]) + W_test_draws[,s]
      }
      P_bar <- P_bar + make_pi_baseline(eta_s)
    }
    P_bar <- P_bar / S

    Y_test_mat <- as.matrix(test_data[, paste0("Y",1:Kcat)])
    nsize_test <- test_data$nsize

    testNLL_postmean <- multinom_nll(Y_test_mat, P_bar)
    pred_counts_bar <- sweep(P_bar, 1, nsize_test, `*`)
    MSPE_count_postmean <- mean((Y_test_mat - pred_counts_bar)^2)

    Yprop_test <- sweep(Y_test_mat, 1, rowSums(Y_test_mat), `/`)
    MSPE_prop_postmean <- mean((Yprop_test - P_bar)^2)

    pp_metrics <- data.frame(
      Kbasis=Kbasis, rank=r_use,
      TestNLL_postmean=testNLL_postmean,
      MSPE_count_postmean=MSPE_count_postmean,
      MSPE_prop_postmean=MSPE_prop_postmean,
      draws_used=S, thin=thin_pp, runtime_sec=runtime_sec
    )
    write.csv(pp_metrics, file.path(per_dir, "PosteriorPredictive_TestMetrics.csv"), row.names=FALSE)

    pdf(file.path(per_dir, "PosteriorPredictive_Plots.pdf"), width=10, height=7)
    prop_hat <- P_bar
    prop_obs <- Yprop_test
    pp_long <- data.frame(
      obs=as.vector(prop_obs),
      pred=as.vector(prop_hat),
      cat=factor(rep(paste0("Y",1:Kcat), each=nrow(prop_obs)))
    )
    print(
      ggplot(pp_long, aes(x=pred, y=obs)) +
        geom_point(alpha=0.4) +
        facet_wrap(~cat, scales="free") +
        geom_abline(slope=1, intercept=0) +
        theme_minimal() +
        ggtitle(paste("Posterior mean predicted vs observed proportions | Kbasis",Kbasis,"rank",r_use))
    )
    dev.off()
  }
}

write.csv(rank_all, file.path(out_dir, "RankSelection_AllKbasis_AllRanks.csv"), row.names=FALSE)
write.csv(best_by_Kbasis, file.path(out_dir, "RankSelection_BestByKbasis.csv"), row.names=FALSE)

# 11c) GLOBAL ESS TABLES + PLOTS
if (RUN_MCMC_EACH_KBASIS && nrow(ess_all) > 0) {

  write.csv(ess_all, csv_ess_all_file, row.names = FALSE)

  write.csv(ess_summary_by_Kbasis, csv_ess_sum_file, row.names = FALSE)

  ess_group_summary <- ess_all %>%
    group_by(Kbasis, Group) %>%
    summarise(
      Rank = dplyr::first(Rank),
      RuntimeSec = dplyr::first(RuntimeSec),
      n_params = dplyr::n(),
      ESS_mean = mean(ESS, na.rm=TRUE),
      ESS_median = median(ESS, na.rm=TRUE),
      ESS_min = min(ESS, na.rm=TRUE),
      ESSps_mean = mean(ESS_per_sec, na.rm=TRUE),
      ESSps_median = median(ESS_per_sec, na.rm=TRUE),
      ESSps_min = min(ESS_per_sec, na.rm=TRUE),
      .groups="drop"
    )
  write.csv(ess_group_summary, csv_ess_group_file, row.names = FALSE)

  best_metrics <- merge(best_by_Kbasis, ess_summary_by_Kbasis, by="Kbasis", all.x=TRUE)
  write.csv(best_metrics, csv_best_metrics_file, row.names = FALSE)

  pdf(pdf_ess_all_file, width=11, height=8.5)

  print(
    ggplot(ess_summary_by_Kbasis, aes(x=Kbasis, y=ESSps_mean_all)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      ggtitle("Mean ESS/sec (ALL parameters) vs Kbasis (best rank per Kbasis)")
  )

  print(
    ggplot(ess_summary_by_Kbasis, aes(x=Kbasis, y=ESSps_min_all)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      ggtitle("Min ESS/sec (ALL parameters) vs Kbasis (best rank per Kbasis)")
  )

  print(
    ggplot(ess_summary_by_Kbasis, aes(x=Kbasis, y=ESSps_mean_fixed)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      ggtitle("Mean ESS/sec (alpha+beta+tau only) vs Kbasis (best rank per Kbasis)")
  )

  print(
    ggplot(ess_all, aes(x=factor(Kbasis), y=ESS_per_sec)) +
      geom_boxplot() +
      facet_wrap(~Group, scales="free_y") +
      theme_minimal() +
      ggtitle("ESS/sec distribution (ALL parameters) by Kbasis and parameter group")
  )

  dev.off()

  pdf(pdf_runtime_ess_file, width=10, height=7)

  print(
    ggplot(best_metrics, aes(x=RuntimeSec, y=ESSps_mean_all)) +
      geom_point(size=2) +
      geom_text(aes(label=Kbasis), vjust=-0.8, size=3) +
      theme_minimal() +
      ggtitle("Best-rank MCMC runtime vs mean ESS/sec (ALL params) | labels = Kbasis")
  )

  print(
    ggplot(best_metrics, aes(x=RuntimeSec, y=ESSps_min_all)) +
      geom_point(size=2) +
      geom_text(aes(label=Kbasis), vjust=-0.8, size=3) +
      theme_minimal() +
      ggtitle("Best-rank MCMC runtime vs min ESS/sec (ALL params) | labels = Kbasis")
  )

  print(
    ggplot(best_metrics, aes(x=RuntimeSec, y=ESSps_mean_fixed)) +
      geom_point(size=2) +
      geom_text(aes(label=Kbasis), vjust=-0.8, size=3) +
      theme_minimal() +
      ggtitle("Best-rank MCMC runtime vs mean ESS/sec (alpha+beta+tau) | labels = Kbasis")
  )

  dev.off()
}

# 12) GLOBAL PLOTS ACROSS ALL Kbasis
pdf(pdf_rank_all_file, width=11, height=8.5)

print(
  ggplot(rank_all, aes(x=ST_Rank, y=TestNLL)) +
    geom_line() + geom_point(size=1.2) +
    facet_wrap(~Kbasis, scales="free_x") +
    theme_minimal() +
    ggtitle("Test NLL vs Rank (VGAM) for each Kbasis")
)

print(
  ggplot(rank_all, aes(x=ST_Rank, y=MSPE_count)) +
    geom_line() + geom_point(size=1.2) +
    facet_wrap(~Kbasis, scales="free_x") +
    theme_minimal() +
    ggtitle("MSPE(count) vs Rank for each Kbasis")
)

print(
  ggplot(rank_all, aes(x=ST_Rank, y=TimeSec)) +
    geom_line() + geom_point(size=1.2) +
    facet_wrap(~Kbasis, scales="free_x") +
    theme_minimal() +
    ggtitle("Runtime vs Rank for each Kbasis")
)

dev.off()

# Coefficient vs rank across all Kbasis
coef_cols <- c(
  grep("^alpha[0-9]+$", names(rank_all), value=TRUE),
  grep("^b_", names(rank_all), value=TRUE)
)
coef_long <- rank_all %>%
  select(Kbasis, ST_Rank, all_of(coef_cols)) %>%
  pivot_longer(cols = -c(Kbasis, ST_Rank), names_to="param", values_to="est")

truth_df <- data.frame(param=character(), truth=numeric(), stringsAsFactors=FALSE)

for (k in 1:Km1) truth_df <- rbind(truth_df, data.frame(param=paste0("alpha",k), truth=alpha_true[k]))

for (k in 1:Km1) {
  for (j in 1:pX) {
    truth_df <- rbind(truth_df, data.frame(
      param=paste0("b_", cov_names[j], "_", k),
      truth=beta_true[k,j]
    ))
  }
}

coef_long <- merge(coef_long, truth_df, by="param", all.x=TRUE)

pdf(pdf_coef_all_file, width=12, height=9)
print(
  ggplot(coef_long, aes(x=ST_Rank, y=est)) +
    geom_line() + geom_point(size=1.1) +
    geom_hline(aes(yintercept=truth), linetype=2) +
    facet_grid(param ~ Kbasis, scales="free_y") +
    theme_minimal() +
    ggtitle("VGAM estimates vs Rank (dashed = true) for alphas and all betas")
)
dev.off()

# Best-by-Kbasis plots
pdf(pdf_bestK_file, width=10, height=7)
print(ggplot(best_by_Kbasis, aes(x=Kbasis, y=best_TestNLL)) + geom_line() + geom_point(size=2) +
        theme_minimal() + ggtitle("Best TestNLL vs Kbasis"))
print(ggplot(best_by_Kbasis, aes(x=Kbasis, y=best_MSPE_count)) + geom_line() + geom_point(size=2) +
        theme_minimal() + ggtitle("Best MSPE(count) vs Kbasis"))

if (RUN_MCMC_EACH_KBASIS && exists("best_metrics") && nrow(best_metrics) > 0) {
  print(ggplot(best_metrics, aes(x=Kbasis, y=RuntimeSec)) +
          geom_line() + geom_point(size=2) +
          theme_minimal() + ggtitle("Best-rank TOTAL MCMC runtime vs Kbasis (NIMBLE)"))
} else {
  print(ggplot(best_by_Kbasis, aes(x=Kbasis, y=best_TimeSec)) +
          geom_line() + geom_point(size=2) +
          theme_minimal() + ggtitle("Best-rank VGAM runtime vs Kbasis (MCMC not run)"))
}

print(ggplot(best_by_Kbasis, aes(x=Kbasis, y=best_rank)) + geom_line() + geom_point(size=2) +
        theme_minimal() + ggtitle("Best selected rank vs Kbasis"))

# ESS plots inside best metrics PDF
if (RUN_MCMC_EACH_KBASIS && exists("best_metrics") && nrow(best_metrics) > 0) {

  print(
    ggplot(best_metrics, aes(x=Kbasis, y=ESSps_mean_all)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() + ggtitle("Best-rank mean ESS/sec (ALL params) vs Kbasis (MCMC)")
  )

  print(
    ggplot(best_metrics, aes(x=Kbasis, y=ESSps_min_all)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() + ggtitle("Best-rank min ESS/sec (ALL params) vs Kbasis (MCMC)")
  )

  print(
    ggplot(best_metrics, aes(x=RuntimeSec, y=ESSps_mean_all)) +
      geom_point(size=2) +
      geom_text(aes(label=Kbasis), vjust=-0.8, size=3) +
      theme_minimal() +
      ggtitle("Best-rank runtime vs mean ESS/sec (ALL params) | labels=Kbasis")
  )
}
dev.off()

# 13) CV ON THE BEST Kbasis
if (RUN_CV_BEST_KBASIS) {
  bestK <- best_by_Kbasis$Kbasis[which.min(best_by_Kbasis$best_MSPE_count)]
  cat("\nRunning CV only for best Kbasis =", bestK, "\n")

  kb_best <- build_Kbasis(full_ST$ST_train_full, full_ST$ST_test_full, full_ST$MQM_full,
                          Kbasis=bestK, seed=1000 + bestK)
  Bc_best <- center_basis_train_test(kb_best$B_train, kb_best$B_test)
  B_train_best <- Bc_best$B_train

  Ytr_df <- as.data.frame(train_data[, paste0("Y",1:Kcat)])
  cov_tr_df <- train_data[, cov_names, drop=FALSE]
  rank_seq_cv <- make_rank_grid(ncol(B_train_best))

  cv_out <- run_kfold_cv_vgam(Kcat=Kcat, cov_names=cov_names,
                              Y=Ytr_df, cov_df=cov_tr_df, B_full=B_train_best,
                              Kfold=Kfold, rank_sequence=rank_seq_cv, seed=123)

  write.csv(cv_out$cv_results, file.path(out_dir, "CV_BestKbasis_Results.csv"), row.names=FALSE)

  pdf(pdf_cv_bestK_file, width=9, height=6)
  print(
    ggplot(cv_out$cv_results, aes(x=ST_Rank, y=CV_NLL)) +
      geom_line() + geom_point() +
      theme_minimal() +
      ggtitle(paste("K-fold CV NLL vs Rank | best Kbasis =", bestK))
  )
  dev.off()
}

# 14) SAVE EVERYTHING
save(
  Kcat, Km1, Kbasis_vals, p_full, q_full, cov_spec, cov_names, alpha_true, beta_true,
  mesh, AMat_train, AMat_test, W_adj, Q_lap, mBase, tBase_train, tBase_test,
  full_ST, tau_true, delta_true, Kbasis_truth,
  train_data, test_data, rank_all, best_by_Kbasis,
  ess_all, ess_summary_by_Kbasis,
  file = file.path(out_dir, "All_Objects_Multinomial_SVD_PICAR_Kbasis200_600.RData")
)

overall_end <- Sys.time()
cat(sprintf("\nDONE. Total runtime: %.2f sec\nOutputs in: %s\n",
            as.numeric(difftime(overall_end, overall_start, units="secs")),
            out_dir))

cat("\nGlobal PDFs:\n",
    pdf_mesh_file, "\n",
    pdf_rank_all_file, "\n",
    pdf_coef_all_file, "\n",
    pdf_bestK_file, "\n",
    if (RUN_CV_BEST_KBASIS) pdf_cv_bestK_file else "", "\n",
    if (RUN_MCMC_EACH_KBASIS) pdf_ess_all_file else "", "\n",
    if (RUN_MCMC_EACH_KBASIS) pdf_runtime_ess_file else "", "\n")

# 15) OPTION 1

pdf08_file        <- file.path(out_dir, "08_SVD_PICAR_Multinomial_Oversmoothing_and_MeshDecoupling.pdf")
csv_meshsens_file <- file.path(out_dir, "MeshDensity_Sensitivity_SVD_PICAR_Multinomial.csv")
csv_ranksweep     <- file.path(out_dir, "Roughness_ByRank_SVD_PICAR_Multinomial.csv")

# Packages used in this block
need_pkgs <- c("grid", "gridExtra", "ggplot2")
for (pp in need_pkgs) {
  if (!requireNamespace(pp, quietly = TRUE)) install.packages(pp)
}
library(grid)
library(gridExtra)
library(ggplot2)

# USER CONTROLS
t_star <- median(train_data$year_mo)

NITER_RANKSWEEP <- 30000
BURN_FRAC_SWEEP <- 0.60

Kbasis_sweep   <- sort(unique(Kbasis_vals))

RUN_MESH_PILOT      <- TRUE
maxedge_grid        <- c(0.8, 1.0, 1.2, 1.5, 1.8, 2.2)
NITER_PILOT_MESH    <- 5000
BURN_FRAC_PILOT     <- 0.50

Kbasis_pilot <- best_by_Kbasis$Kbasis[which.min(best_by_Kbasis$best_MSPE_count)]
rank_pilot   <- best_by_Kbasis$best_rank[best_by_Kbasis$Kbasis == Kbasis_pilot][1]

# Helpers
get_mesh_edges <- function(mesh) {
  tri <- mesh$graph$tv
  e <- rbind(tri[, c(1,2)], tri[, c(1,3)], tri[, c(2,3)])
  e <- t(apply(e, 1, function(x) sort(x)))
  e <- unique(e)
  colnames(e) <- c("i","j")
  e
}

roughness_R_L2 <- function(u, edges) {
  d <- u[edges[,1]] - u[edges[,2]]
  R  <- mean(d^2, na.rm=TRUE)
  L2 <- sqrt(R)
  list(R=R, L2=L2)
}

bs_at <- function(tval, tBase_train, q_full) {
  knots <- attr(tBase_train, "knots")
  bnd   <- attr(tBase_train, "Boundary.knots")
  as.numeric(splines::bs(tval, df=q_full, intercept=FALSE, knots=knots, Boundary.knots=bnd))
}

make_rank_grid <- function(Kb) {
  base <- unique(sort(c(2, 3, 5, 8, 10, 15, 20, 30, 40, 50, 75, 100)))
  base <- base[base <= Kb]
  if (length(base) == 0) base <- unique(sort(c(2, Kb)))
  unique(sort(c(base, Kb)))
}

build_kbasis_cached <- function(Kbasis, full_ST) {
  kb <- build_Kbasis(full_ST$ST_train_full, full_ST$ST_test_full, full_ST$MQM_full,
                     Kbasis=Kbasis, seed=1000 + Kbasis)

  Bc <- center_basis_train_test(kb$B_train, kb$B_test)

  list(
    Kbasis  = Kbasis,
    kb      = kb,
    B_train = Bc$B_train,
    B_test  = Bc$B_test,
    center  = Bc$center
  )
}

compute_u_vertices_at_t_cached <- function(kbc, r_use, delta_mean, t_star,
                                           mesh, full_ST, tBase_train, p_full, q_full) {

  t_sc <- bs_at(t_star, tBase_train=tBase_train, q_full=q_full)
  Nnew <- nrow(mesh$loc)
  t_mat <- matrix(rep(t_sc, each=Nnew), nrow=Nnew, byrow=FALSE)

  if (is.null(full_ST$mBase_sub)) stop("full_ST$mBase_sub is missing (vertex Moran basis).")

  ST_vert_full <- make_st_scores(full_ST$mBase_sub, t_mat)

  B_vert <- ST_vert_full %*% kbc$kb$V
  B_vert <- sweep(B_vert, 2, kbc$center, "-")

  u_vert <- as.vector(B_vert[, 1:r_use, drop=FALSE] %*% as.numeric(delta_mean[1:r_use]))
  u_vert
}

pdf_text_page <- function(title, bullets) {
  grid::grid.newpage()
  grid::grid.text(title, x=0.02, y=0.95, just=c("left","top"),
                  gp=grid::gpar(fontsize=16, fontface="bold"))
  grid::grid.text(paste0("\u2022 ", bullets, collapse="\n"),
                  x=0.02, y=0.88, just=c("left","top"),
                  gp=grid::gpar(fontsize=12))
}

get_best_rank_for_Kbasis <- function(Kb, Kb_eff) {
  if (exists("best_by_Kbasis") && nrow(best_by_Kbasis) > 0) {
    ii <- which(best_by_Kbasis$Kbasis == Kb)
    if (length(ii) >= 1) {
      r0 <- as.integer(best_by_Kbasis$best_rank[ii[1]])
      r0 <- max(1L, min(r0, as.integer(Kb_eff)))
      return(r0)
    }
  }
  return(max(2L, min(20L, as.integer(Kb_eff))))
}

# Prepare common inputs
edges0 <- get_mesh_edges(mesh)

if (!exists("Kcat")) stop("Kcat not found. (Number of multinomial categories.)")
Y_train_mat <- as.matrix(train_data[, paste0("Y",1:Kcat), drop=FALSE])
nsize_train <- train_data$nsize

X_train_raw <- as.matrix(train_data[, cov_names, drop=FALSE])
X_train_std <- scale(X_train_raw)

# (1) OPTION 1: Rank sweep for each Kbasis
cat("\n[OPTION 1 | MULTINOMIAL SVD–PICAR] Rank sweep starting...\n")

kb_cache <- list()
for (Kb in Kbasis_sweep) {
  cat("  caching Kbasis =", Kb, "\n")
  kb_cache[[as.character(Kb)]] <- build_kbasis_cached(Kb, full_ST)
}

rank_sweep_res <- data.frame()

delta_best_byKbasis <- list()

for (Kb in Kbasis_sweep) {

  kbc <- kb_cache[[as.character(Kb)]]
  if (is.null(kbc)) stop("kb_cache missing Kbasis=", Kb)

  Kb_eff <- ncol(kbc$B_train)

  r_bestK <- get_best_rank_for_Kbasis(Kb, Kb_eff)

  r_grid <- unique(sort(c(make_rank_grid(Kb_eff), r_bestK)))

  cat("\n  --- Kbasis =", Kb, " (K_eff=", Kb_eff, ") | best rank=", r_bestK,
      " | ranks: {", paste(r_grid, collapse=", "), "} ---\n")

  for (rnk in r_grid) {

    r_use <- min(as.integer(rnk), as.integer(Kb_eff))
    cat("    fitting rank r =", r_use, "...\n")

    B_train_rank <- kbc$B_train[, 1:r_use, drop=FALSE]
    MQM_rank <- kbc$kb$MQM_reduced[1:r_use, 1:r_use, drop=FALSE]
    MQM_rank <- 0.5*(MQM_rank + t(MQM_rank)) + diag(1e-6, nrow(MQM_rank))

    fit <- run_nimble_multinom(
      Kcat=Kcat, cov_names=cov_names,
      Y_train_mat=Y_train_mat,
      X_train_std=X_train_std,
      B_train_rank=B_train_rank,
      MQM_rank=MQM_rank,
      nsize_train=nsize_train,
      niter=NITER_RANKSWEEP,
      burn_frac=BURN_FRAC_SWEEP,
      seed= 70000 + 100*Kb + r_use
    )

    if (is.null(fit$samples)) stop("run_nimble_multinom must return $samples for Kbasis=",Kb," rank=",r_use)

    samples <- fit$samples
    burn <- floor(BURN_FRAC_SWEEP * nrow(samples))
    post <- samples[(burn+1):nrow(samples), , drop=FALSE]

    dcols <- paste0("delta[", 1:r_use, "]")
    miss <- setdiff(dcols, colnames(post))
    if (length(miss) > 0) {
      stop("Missing delta columns in posterior for Kbasis=",Kb," rank=",r_use,
           " | missing: ", paste(miss, collapse=", "))
    }
    delta_mean <- colMeans(post[, dcols, drop=FALSE])

    uV <- compute_u_vertices_at_t_cached(
      kbc=kbc, r_use=r_use, delta_mean=delta_mean, t_star=t_star,
      mesh=mesh, full_ST=full_ST, tBase_train=tBase_train, p_full=p_full, q_full=q_full
    )

    rr <- roughness_R_L2(uV, edges0)
    sf <- rr$L2 / max(sd(uV), 1e-12)

    rank_sweep_res <- rbind(rank_sweep_res, data.frame(
      Kbasis       = Kb,
      rank         = r_use,
      R_hat        = rr$R,
      RMS_hat      = rr$L2,
      sf_hat       = sf,
      runtime_sec  = if (!is.null(fit$runtime_sec)) fit$runtime_sec else NA_real_,
      time_per_iter= if (!is.null(fit$runtime_sec)) fit$runtime_sec / max(NITER_RANKSWEEP,1) else NA_real_,
      is_best_rank = (r_use == r_bestK),
      stringsAsFactors=FALSE
    ))

    if (r_use == r_bestK) {
      delta_best_byKbasis[[as.character(Kb)]] <- list(rank=r_bestK, delta_mean=delta_mean)
    }
  }
}

write.csv(rank_sweep_res, csv_ranksweep, row.names=FALSE)
cat("\nSaved rank-sweep CSV:", csv_ranksweep, "\n")

# (2) Mesh-density decoupling pilot
mesh_sens <- NULL

if (RUN_MESH_PILOT) {

  train_coords <- as.matrix(train_data[, c("lon","lat")])
  test_coords  <- as.matrix(test_data[,  c("lon","lat")])

  mesh_sens <- data.frame()

  for (me0 in maxedge_grid) {

    bt <- system.time({

      mesh_i <- INLA::inla.mesh.2d(
        loc      = train_coords,
        max.edge = c(1,2) * me0,
        cutoff   = me0/5,
        offset   = c(me0, 6.0)
      )

      AMat_train_i <- INLA::inla.spde.make.A(mesh_i, loc=train_coords)
      AMat_test_i  <- INLA::inla.spde.make.A(mesh_i, loc=test_coords)

      lap_i <- mesh_to_laplacian(mesh_i)
      W_i   <- lap_i$W
      Q_i   <- as.matrix(lap_i$Q) + diag(1e-6, lap_i$Nnew)

      mBase_i <- build_moran_basis(W_i, num_eigen = p_full)

      tBase_train_i <- splines::bs(train_data$year_mo, df=q_full, intercept=FALSE)
      tBase_test_i  <- splines::bs(test_data$year_mo,  df=q_full, intercept=FALSE)

      full_ST_i <- build_full_st_once(
        AMat_train=AMat_train_i, AMat_test=AMat_test_i,
        Q_lap=Q_i, mBase=mBase_i,
        tBase_train=tBase_train_i, tBase_test=tBase_test_i,
        p_full=p_full, q_full=q_full
      )

      kb_i <- build_Kbasis(full_ST_i$ST_train_full, full_ST_i$ST_test_full, full_ST_i$MQM_full,
                           Kbasis=Kbasis_pilot, seed=1000 + Kbasis_pilot)

      Bc_i <- center_basis_train_test(kb_i$B_train, kb_i$B_test)
      B_train_full_i <- Bc_i$B_train

      r_use_i <- min(rank_pilot, ncol(B_train_full_i))
      B_train_rank_i <- B_train_full_i[, 1:r_use_i, drop=FALSE]
      MQM_rank_i <- kb_i$MQM_reduced[1:r_use_i, 1:r_use_i, drop=FALSE]
      MQM_rank_i <- 0.5*(MQM_rank_i + t(MQM_rank_i)) + diag(1e-6, nrow(MQM_rank_i))
    })

    fit_i <- run_nimble_multinom(
      Kcat=Kcat, cov_names=cov_names,
      Y_train_mat=Y_train_mat,
      X_train_std=X_train_std,
      B_train_rank=B_train_rank_i,
      MQM_rank=MQM_rank_i,
      nsize_train=nsize_train,
      niter=NITER_PILOT_MESH,
      burn_frac=BURN_FRAC_PILOT,
      seed=900 + round(100*me0)
    )

    Nnew_i <- nrow(mesh_i$loc)
    ntri_i <- nrow(mesh_i$graph$tv)

    mesh_sens <- rbind(mesh_sens, data.frame(
      max_edge0 = me0,
      n_vertices = Nnew_i,
      n_triangles = ntri_i,
      Kbasis = Kbasis_pilot,
      Rank   = min(rank_pilot, ncol(B_train_rank_i)),
      build_time_sec = unname(bt["elapsed"]),
      mcmc_runtime_sec = fit_i$runtime_sec,
      mcmc_time_per_iter = fit_i$runtime_sec / max(NITER_PILOT_MESH, 1)
    ))

    cat(sprintf("Mesh pilot | max.edge0=%.2f | vertices=%d | build=%.2fs | MCMC/iter=%.6fs\n",
                me0, Nnew_i, unname(bt["elapsed"]), fit_i$runtime_sec / max(NITER_PILOT_MESH,1)))
  }

  write.csv(mesh_sens, csv_meshsens_file, row.names=FALSE)
  cat("\nSaved mesh-sensitivity CSV:", csv_meshsens_file, "\n")
}

# WRITE THE PDF
pdf(pdf08_file, width=11, height=8.5)

pdf_text_page(
  "08 — SVD–PICAR (Multinomial): Option 1 Rank Sweep + Mesh-density Decoupling",
  c(
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Time slice: t* = ", sprintf("%.4f", t_star)),
    paste0("Kbasis swept: {", paste(Kbasis_sweep, collapse=", "), "}"),
    paste0("Rank-sweep MCMC: niter=", NITER_RANKSWEEP, ", burn=", BURN_FRAC_SWEEP),
    paste0("Rank-sweep CSV: ", csv_ranksweep),
    paste0("Mesh pilot CSV: ", csv_meshsens_file)
  )
)

# Roughness vs rank
if (nrow(rank_sweep_res) > 0) {

  rank_sweep_res$logRMS <- log10(rank_sweep_res$RMS_hat)

  print(
    ggplot(rank_sweep_res, aes(x=rank, y=logRMS)) +
      geom_line() + geom_point(size=1.6) +
      facet_wrap(~Kbasis, scales="free_x") +
      theme_minimal() +
      labs(
        title="log10( RMS roughness ) vs rank r  (higher = finer-scale / less smooth)",
        x="Rank r",
        y="log10(RMS roughness)"
      )
  )

  print(
    ggplot(rank_sweep_res, aes(x=rank, y=sf_hat)) +
      geom_line() + geom_point(size=1.6) +
      facet_wrap(~Kbasis, scales="free_x") +
      theme_minimal() +
      labs(
        title="Scale-free roughness (RMS/sd) vs rank r",
        x="Rank r",
        y="RMS(u)/sd(u)"
      )
  )
}

# LATENT MAPS (per your instruction)
pdf_text_page(
  "UPDATED MAPS — One panel per Kbasis at its best-MSPE/TestNLL rank (no eigenvector labels)",
  c(
    "Replaces the old 'Kbasis_display + a few ranks' map section.",
    "Each panel shows u(v,t*) on mesh vertices, using posterior-mean delta at r*(Kbasis).",
    "All panels share a common color scale for direct comparison across Kbasis."
  )
)

latent_byK <- data.frame()

for (Kb in Kbasis_sweep) {

  if (is.null(delta_best_byKbasis[[as.character(Kb)]])) {
    warning(paste0("Missing stored delta_mean at best rank for Kbasis=", Kb,
                   ". (Check that best rank was included in the sweep grid.)"))
    next
  }

  kbc <- kb_cache[[as.character(Kb)]]
  r_bestK <- delta_best_byKbasis[[as.character(Kb)]]$rank
  delta_mean_best <- delta_best_byKbasis[[as.character(Kb)]]$delta_mean

  uV_best <- compute_u_vertices_at_t_cached(
    kbc=kbc, r_use=r_bestK, delta_mean=delta_mean_best, t_star=t_star,
    mesh=mesh, full_ST=full_ST, tBase_train=tBase_train, p_full=p_full, q_full=q_full
  )

  latent_byK <- rbind(
    latent_byK,
    data.frame(
      lon   = mesh$loc[,1],
      lat   = mesh$loc[,2],
      u_hat = uV_best,
      Kbasis= Kb,
      r_best= r_bestK,
      stringsAsFactors = FALSE
    )
  )
}

latent_byK$K_label <- paste0("K = ", latent_byK$Kbasis, "   (r* = ", latent_byK$r_best, ")")
latent_byK$K_label <- factor(latent_byK$K_label,
                             levels = unique(latent_byK$K_label[order(latent_byK$Kbasis)]))

lims_u <- range(latent_byK$u_hat, finite = TRUE)

p_panels <- ggplot(latent_byK, aes(x=lon, y=lat, color=u_hat)) +
  geom_point(size=0.55) +
  coord_equal() +
  facet_wrap(~K_label, ncol=3) +
  theme_minimal() +
  labs(
    title    = "Multinomial SVD–PICAR: posterior-mean latent field u(v,t*) at best rank per Kbasis",
    subtitle = paste0("t* = ", sprintf("%.4f", t_star), "   |   one panel per Kbasis (common color scale)"),
    x="Longitude", y="Latitude", color="u_hat"
  ) +
  scale_color_gradientn(
    colours = c("navy","skyblue","yellow","orange","red"),
    limits  = lims_u
  )

print(p_panels)

# Mesh-density decoupling pages
if (RUN_MESH_PILOT && !is.null(mesh_sens) && nrow(mesh_sens) > 0) {

  pdf_text_page(
    "Mesh-density decoupling pilot",
    c(
      paste0("Pilot fixed: Kbasis=", Kbasis_pilot, ", rank=", rank_pilot),
      paste0("Pilot MCMC: niter=", NITER_PILOT_MESH, ", burn=", BURN_FRAC_PILOT),
      paste0("max.edge0 grid: {", paste(maxedge_grid, collapse=", "), "}")
    )
  )

  print(
    ggplot(mesh_sens, aes(x=n_vertices, y=build_time_sec)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(
        title=paste0("One-time basis build time vs mesh density | Kbasis=", Kbasis_pilot, ", rank=", rank_pilot),
        x="# mesh vertices", y="Build time (sec)"
      )
  )

  print(
    ggplot(mesh_sens, aes(x=n_vertices, y=mcmc_time_per_iter)) +
      geom_line() + geom_point(size=2) +
      theme_minimal() +
      labs(
        title=paste0("MCMC time per iteration vs mesh density | Kbasis=", Kbasis_pilot, ", rank=", rank_pilot),
        x="# mesh vertices", y="Seconds per MCMC iteration"
      )
  )

  tab_show <- mesh_sens
  tab_show$build_time_sec <- round(tab_show$build_time_sec, 2)
  tab_show$mcmc_runtime_sec <- round(tab_show$mcmc_runtime_sec, 2)
  tab_show$mcmc_time_per_iter <- signif(tab_show$mcmc_time_per_iter, 4)
  grid::grid.newpage()
  gridExtra::grid.table(tab_show)
}

dev.off()

cat("\n[SECTION 15 OPTION 1 COMPLETE]\n",
    "PDF: ", pdf08_file, "\n",
    "CSV (rank sweep): ", csv_ranksweep, "\n",
    if (RUN_MESH_PILOT) paste0("CSV (mesh pilot): ", csv_meshsens_file, "\n") else "",
    sep="")

# END Section 15
