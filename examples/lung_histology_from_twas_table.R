#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript examples/lung_histology_from_twas_table.R [--input=FILE] --out_dir=DIR\n\n",
    "If --input is omitted, a deterministic demo table is generated in memory.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list(input = NA_character_, out_dir = "example_output/lung_histology_demo")
  for (arg in args) {
    if (arg %in% c("-h", "--help")) {
      usage()
      quit(status = 0)
    } else if (grepl("^--input=", arg)) {
      out$input <- sub("^--input=", "", arg)
    } else if (grepl("^--out_dir=", arg)) {
      out$out_dir <- sub("^--out_dir=", "", arg)
    } else {
      stop("Unknown argument: ", arg)
    }
  }
  out
}

load_harmoni <- function() {
  if (requireNamespace("HARMONI", quietly = TRUE)) {
    suppressPackageStartupMessages(library(HARMONI))
    return(invisible(TRUE))
  }

  script <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script, value = TRUE)
  start <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
  root <- normalizePath(file.path(start, ".."), mustWork = FALSE)
  desc <- file.path(root, "DESCRIPTION")
  r_dir <- file.path(root, "R")
  if (!file.exists(desc) || !dir.exists(r_dir)) {
    stop("HARMONI is not installed and the script could not locate the package source tree.")
  }

  for (f in sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
  invisible(TRUE)
}

make_demo_lung_table <- function(n_features = 120L, seed = 1001L) {
  set.seed(seed)
  id <- sprintf("GENE%03d", seq_len(n_features))

  z <- matrix(rnorm(n_features * 4L), ncol = 4L)
  colnames(z) <- c("lungoverall", "luad", "sqc", "sclc")

  shared <- 81:92
  luad_only <- 93:102
  opposite <- 103:112

  z[shared, ] <- z[shared, ] + matrix(rep(c(3.0, 2.8, 2.7, 2.5), length(shared)), ncol = 4L, byrow = TRUE)
  z[luad_only, "luad"] <- z[luad_only, "luad"] + 3.4
  z[luad_only, "lungoverall"] <- z[luad_only, "lungoverall"] + 1.1
  z[opposite, "luad"] <- z[opposite, "luad"] + 3.2
  z[opposite, "sqc"] <- z[opposite, "sqc"] - 3.0
  z[opposite, "lungoverall"] <- z[opposite, "lungoverall"] + 0.2

  p <- 2 * pnorm(-abs(z))
  fdr <- apply(p, 2L, p.adjust, method = "BH")

  data.frame(
    ID = id,
    TWAS.Z_lungoverall = z[, "lungoverall"],
    TWAS.Z_luad = z[, "luad"],
    TWAS.Z_sqc = z[, "sqc"],
    TWAS.Z_sclc = z[, "sclc"],
    BH_FDR_lungoverall = fdr[, "lungoverall"],
    BH_FDR_luad = fdr[, "luad"],
    BH_FDR_sqc = fdr[, "sqc"],
    BH_FDR_sclc = fdr[, "sclc"],
    stringsAsFactors = FALSE
  )
}

read_lung_table <- function(input) {
  if (is.na(input)) return(make_demo_lung_table())
  read.delim(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

select_null_indices <- function(z_matrix, threshold = 2.5, min_n = 30L) {
  max_abs <- apply(abs(z_matrix), 1L, max, na.rm = TRUE)
  idx <- which(max_abs < threshold)
  if (length(idx) >= min_n) return(idx)
  order(max_abs)[seq_len(min(min_n, nrow(z_matrix)))]
}

metric_qr <- function(B, Sigma) {
  k <- ncol(B)
  Q <- matrix(0, nrow(B), k)
  colnames(Q) <- colnames(B)
  rownames(Q) <- rownames(B)

  for (j in seq_len(k)) {
    v <- B[, j]
    if (j > 1L) {
      for (i in seq_len(j - 1L)) {
        qi <- Q[, i]
        denom <- as.numeric(t(qi) %*% Sigma %*% qi)
        if (denom > 0) {
          v <- v - as.numeric((t(qi) %*% Sigma %*% v) / denom) * qi
        }
      }
    }
    norm_v <- sqrt(as.numeric(t(v) %*% Sigma %*% v))
    if (!is.finite(norm_v) || norm_v <= 0) stop("Could not normalize contrast column: ", colnames(B)[j])
    Q[, j] <- v / norm_v
  }
  Q
}

normalize_columns <- function(B, Sigma) {
  Q <- B
  for (j in seq_len(ncol(Q))) {
    norm_j <- sqrt(as.numeric(t(Q[, j]) %*% Sigma %*% Q[, j]))
    if (!is.finite(norm_j) || norm_j <= 0) stop("Could not normalize contrast column: ", colnames(B)[j])
    Q[, j] <- Q[, j] / norm_j
  }
  Q
}

run_parameterization <- function(name, dat, z_orig, Sigma0, shared_vec, het_basis,
                                 primary_baseline = c("union_all", "overall", "union_sub"),
                                 mode = c("orth", "simple")) {
  primary_baseline <- match.arg(primary_baseline)
  mode <- match.arg(mode)

  B <- cbind(shared = shared_vec, het_basis)
  rownames(B) <- colnames(z_orig)

  if (mode == "orth") {
    Q <- metric_qr(B, Sigma0)
    colnames(Q) <- c("shared", paste0("het_", seq_len(ncol(Q) - 1L)))
  } else {
    Q <- normalize_columns(B, Sigma0)
    colnames(Q) <- c("shared", paste0("het_", seq_len(ncol(Q) - 1L)))
  }

  z_matrix <- z_orig %*% Q
  z_matrix <- as.matrix(z_matrix)
  rownames(z_matrix) <- dat$FEATURE_KEY
  colnames(z_matrix) <- colnames(Q)

  null_idx <- select_null_indices(z_matrix)
  R_null <- cor(z_matrix[null_idx, , drop = FALSE], use = "pairwise.complete.obs")
  R_null <- make_psd(R_null)
  colnames(R_null) <- rownames(R_null) <- colnames(z_matrix)

  cfgs <- build_configs(ncol(z_matrix), colnames(z_matrix))
  ctrl <- ms_bpwas_control(method = "vb", tol = 1e-5, max_iter = 500L, verbose = FALSE)
  fit <- fit_vb(z_matrix, R_null, cfgs, ctrl)
  der <- compute_derived(fit, z_matrix, R_null, cfgs, colnames(z_matrix), ctrl)

  res <- data.frame(
    ID = dat$ID,
    FEATURE_KEY = dat$FEATURE_KEY,
    TWAS.Z_lungoverall = dat$TWAS.Z_lungoverall,
    TWAS.Z_luad = dat$TWAS.Z_luad,
    TWAS.Z_sqc = dat$TWAS.Z_sqc,
    TWAS.Z_sclc = dat$TWAS.Z_sclc,
    BH_FDR_lungoverall = dat$BH_FDR_lungoverall,
    BH_FDR_luad = dat$BH_FDR_luad,
    BH_FDR_sqc = dat$BH_FDR_sqc,
    BH_FDR_sclc = dat$BH_FDR_sclc,
    p_null = der$p_null[dat$FEATURE_KEY],
    bfdr = der$bfdr[dat$FEATURE_KEY],
    map_config = der$map_config[dat$FEATURE_KEY],
    stringsAsFactors = FALSE
  )

  for (nm in colnames(z_matrix)) res[[nm]] <- z_matrix[, nm]
  for (nm in colnames(z_matrix)) res[[paste0("pip_", nm)]] <- der$pip[dat$FEATURE_KEY, nm]

  sig_harmoni <- res$FEATURE_KEY[res$bfdr < 0.05]
  sig_overall <- res$FEATURE_KEY[res$BH_FDR_lungoverall < 0.05]
  sig_union_sub <- res$FEATURE_KEY[
    (res$BH_FDR_luad < 0.05) |
      (res$BH_FDR_sqc < 0.05) |
      (res$BH_FDR_sclc < 0.05)
  ]
  sig_union_all <- res$FEATURE_KEY[
    (res$BH_FDR_lungoverall < 0.05) |
      (res$BH_FDR_luad < 0.05) |
      (res$BH_FDR_sqc < 0.05) |
      (res$BH_FDR_sclc < 0.05)
  ]
  primary <- switch(
    primary_baseline,
    overall = sig_overall,
    union_sub = sig_union_sub,
    union_all = sig_union_all
  )

  summary_tab <- data.frame(
    parameterization = name,
    n_features = nrow(res),
    n_null_cov = attr(Sigma0, "n_null_cov"),
    n_null_transformed = length(null_idx),
    vb_iter = fit$convergence$n_iter,
    vb_converged = fit$convergence$converged,
    pi0 = fit$hyperparams$pi0,
    xi = paste(sprintf("%s=%.3f", names(fit$hyperparams$xi), fit$hyperparams$xi), collapse = ";"),
    bfdr_sig = sum(res$bfdr < 0.05, na.rm = TRUE),
    baseline_primary = primary_baseline,
    baseline_primary_sig = length(primary),
    overlap_primary = length(intersect(sig_harmoni, primary)),
    harmoni_only_primary = length(setdiff(sig_harmoni, primary)),
    primary_only = length(setdiff(primary, sig_harmoni)),
    overall_sig = length(sig_overall),
    union_sub_sig = length(sig_union_sub),
    union_all_sig = length(sig_union_all),
    stringsAsFactors = FALSE
  )

  cfg_counts <- as.data.frame(table(res$map_config[res$bfdr < 0.05]), stringsAsFactors = FALSE)
  names(cfg_counts) <- c("map_config", "n_sig")
  if (!nrow(cfg_counts)) {
    cfg_counts <- data.frame(map_config = character(), n_sig = integer(), stringsAsFactors = FALSE)
  }

  loadings <- data.frame(coefficient = rownames(Q), Q, stringsAsFactors = FALSE)

  list(
    results = res,
    summary = summary_tab,
    cfg_counts = cfg_counts,
    loadings = loadings,
    fit = fit,
    derived = der,
    z_matrix = z_matrix,
    R_null = R_null,
    cfgs = cfgs
  )
}

write_parameterization <- function(obj, out_dir, name) {
  write.table(obj$results, file.path(out_dir, paste0(name, "_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(obj$summary, file.path(out_dir, paste0(name, "_summary.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(obj$cfg_counts, file.path(out_dir, paste0(name, "_config_counts.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(obj$loadings, file.path(out_dir, paste0(name, "_loadings.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(obj, file.path(out_dir, paste0(name, ".rds")))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  load_harmoni()
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- read_lung_table(args$input)
  required <- c(
    "ID",
    "TWAS.Z_lungoverall", "TWAS.Z_luad", "TWAS.Z_sqc", "TWAS.Z_sclc",
    "BH_FDR_lungoverall", "BH_FDR_luad", "BH_FDR_sqc", "BH_FDR_sclc"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Input is missing required columns: ", paste(missing, collapse = ", "))

  z_orig <- as.matrix(dat[, c("TWAS.Z_lungoverall", "TWAS.Z_luad", "TWAS.Z_sqc", "TWAS.Z_sclc")])
  storage.mode(z_orig) <- "double"
  colnames(z_orig) <- c("overall", "luad", "sqc", "sclc")
  dat$FEATURE_KEY <- make.unique(as.character(dat$ID))
  rownames(z_orig) <- dat$FEATURE_KEY

  keep <- complete.cases(z_orig)
  dat <- dat[keep, , drop = FALSE]
  z_orig <- z_orig[keep, , drop = FALSE]
  if (nrow(z_orig) < 40L) stop("At least 40 complete features are recommended for this example.")

  null_idx_orig <- select_null_indices(z_orig)
  Sigma0 <- cov(z_orig[null_idx_orig, , drop = FALSE], use = "pairwise.complete.obs")
  Sigma0 <- make_psd(Sigma0)
  attr(Sigma0, "n_null_cov") <- length(null_idx_orig)

  het_luad_vs_sqc <- c(0, 1, -1, 0)
  het_sclc_vs_nsclc <- c(0, -0.5, -0.5, 1)

  specs <- list(
    current_std = list(
      shared = c(1, 0, 0, 0),
      het = cbind(het_1 = het_luad_vs_sqc, het_2 = het_sclc_vs_nsclc),
      primary = "union_all",
      mode = "simple"
    ),
    orth_overall = list(
      shared = c(1, 0, 0, 0),
      het = cbind(het_1 = het_luad_vs_sqc, het_2 = het_sclc_vs_nsclc),
      primary = "union_all",
      mode = "orth"
    ),
    subtype_shared = list(
      shared = c(0, 1 / 3, 1 / 3, 1 / 3),
      het = cbind(het_1 = het_luad_vs_sqc, het_2 = het_sclc_vs_nsclc),
      primary = "union_sub",
      mode = "orth"
    )
  )

  summaries <- list()
  for (name in names(specs)) {
    spec <- specs[[name]]
    cat("Running ", name, "\n", sep = "")
    obj <- run_parameterization(
      name = name,
      dat = dat,
      z_orig = z_orig,
      Sigma0 = Sigma0,
      shared_vec = spec$shared,
      het_basis = spec$het,
      primary_baseline = spec$primary,
      mode = spec$mode
    )
    write_parameterization(obj, args$out_dir, name)
    summaries[[name]] <- obj$summary
  }

  summary_all <- do.call(rbind, summaries)
  write.table(summary_all, file.path(args$out_dir, "lung_histology_param_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(Sigma0, file.path(args$out_dir, "lung_histology_Sigma0.tsv"),
              sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

  print(summary_all)
}

main()
