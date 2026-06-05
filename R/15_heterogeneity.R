###############################################################################
# 15_heterogeneity.R
#
# Heterogeneity testing for the multi-subtype framework.
#
# Three levels of testing:
#   1. Qualitative:   do subtypes have different active/null status?
#   2. Quantitative:  given both active, are effect sizes equal?
#   3. Directional:   do effects go in opposite directions?
#
# Works with both parameterizations:
#   A. Direct subtype z-matrix (cc, pap, ...)
#   B. Shared + heterogeneity z-matrix (shared, het_1, het_2, ...)
#
# Also provides:
#   - Global heterogeneity Bayes factor per gene
#   - Omnibus Q-statistic (Cochran-like) from subtype z-scores
#   - Contrast-specific heterogeneity tests
#   - Posterior heterogeneity probabilities from configuration posteriors
#
###############################################################################


#' Comprehensive heterogeneity testing.
#'
#' @param z_original  P x S_orig matrix of ORIGINAL subtype z-statistics
#'                    (cc, pap, ...), even if inference was run on the
#'                    shared+het parameterization. Needed for quantitative
#'                    and directional tests.
#' @param R_original  S_orig x S_orig null correlation among original subtypes
#' @param fit         output from fit_vb / fit_em / fit_mcmc
#' @param der         output from compute_derived
#' @param cfgs        configuration matrix from build_configs
#' @param sn          subtype/column names of the z-matrix used for inference
#' @param parameterization  "direct" if inference was on per-subtype z-matrix,
#'                          "shared_het" if on shared+heterogeneity z-matrix
#' @param contrast_matrix   (S-1) x S contrast matrix (required if
#'                          parameterization = "shared_het")
#' @return list with heterogeneity results (see details)
#' @export
test_heterogeneity <- function(z_original,
                                R_original,
                                fit,
                                der,
                                cfgs,
                                sn,
                                parameterization = c("direct", "shared_het"),
                                contrast_matrix  = NULL) {

  parameterization <- match.arg(parameterization)
  P <- nrow(z_original)
  S_orig <- ncol(z_original)
  sn_orig <- colnames(z_original)

  stopifnot(nrow(R_original) == S_orig,
            ncol(R_original) == S_orig)

  q <- extract_fit_q(fit)   # P x K, with K matching the supplied cfgs

  if (nrow(q) != P || ncol(q) != nrow(cfgs)) {
    stop(sprintf(
      "Dimension mismatch: q is %d x %d, but expected %d x %d from z_original and cfgs.",
      nrow(q), ncol(q), P, nrow(cfgs)
    ))
  }

  results <- list()

  # ===========================================================================
  # 1. POSTERIOR HETEROGENEITY (from configuration posteriors)
  # ===========================================================================

  results$posterior <- posterior_heterogeneity(q, cfgs, sn, parameterization)

  # ===========================================================================
  # 2. QUANTITATIVE HETEROGENEITY (from original subtype z-statistics)
  # ===========================================================================

  results$quantitative <- quantitative_heterogeneity(
    z_original, R_original, sn_orig)

  # ===========================================================================
  # 3. DIRECTIONAL HETEROGENEITY (from original subtype z-statistics)
  # ===========================================================================

  results$directional <- directional_heterogeneity(
    z_original, R_original, sn_orig)

  # ===========================================================================
  # 4. CONTRAST-SPECIFIC TESTS (if using shared+het parameterization)
  # ===========================================================================

  if (parameterization == "shared_het" && !is.null(contrast_matrix)) {
    results$contrast <- contrast_heterogeneity(
      fit, der, cfgs, sn, contrast_matrix, sn_orig)
  }

  # ===========================================================================
  # 5. COMBINED SUMMARY
  # ===========================================================================

  results$summary <- combine_heterogeneity(results, P, sn_orig, der)

  results
}


###############################################################################
# Component functions
###############################################################################


#' Posterior heterogeneity from configuration probabilities.
#'
#' For direct parameterization:
#'   heterogeneous = configs where not all subtypes have the same status
#'   i.e., excluding null (0,0,...,0) and all-active (1,1,...,1)
#'
#' For shared+het parameterization:
#'   heterogeneous = any het component is active
#'
#' @keywords internal
posterior_heterogeneity <- function(q, cfgs, sn, parameterization) {

  P <- nrow(q)
  S <- ncol(cfgs)

  if (parameterization == "direct") {

    # Safe rowSums that handles length-0 and length-1 index vectors
    safe_rowSums_q <- function(indices) {
      if (length(indices) == 0) return(rep(0, P))
      if (length(indices) == 1) return(q[, indices])
      rowSums(q[, indices, drop = FALSE])
    }

    # Global: P(not all subtypes same)
    null_cfg <- which(rowSums(cfgs) == 0)
    all_cfg  <- which(rowSums(cfgs) == S)
    p_homogeneous <- safe_rowSums_q(null_cfg) + safe_rowSums_q(all_cfg)
    p_het_global <- 1 - p_homogeneous

    # Bayes factor for heterogeneity
    bf_het <- p_het_global / pmax(p_homogeneous, 1e-300)

    # Pairwise: P(subtype s1 and s2 have different status)
    pairs <- combn(S, 2, simplify = FALSE)
    pair_names <- sapply(pairs, function(p) paste(sn[p[1]], "vs", sn[p[2]]))
    p_het_pair <- matrix(NA_real_, P, length(pairs),
                          dimnames = list(NULL, pair_names))

    for (k in seq_along(pairs)) {
      s1 <- pairs[[k]][1]
      s2 <- pairs[[k]][2]
      diff_cfgs <- which(cfgs[, s1] != cfgs[, s2])
      p_het_pair[, k] <- safe_rowSums_q(diff_cfgs)
    }

    # Per-subtype specificity: P(only this subtype active)
    p_specific <- matrix(NA_real_, P, S, dimnames = list(NULL, sn))
    for (s in 1:S) {
      specific_cfgs <- which(cfgs[, s] == 1 & rowSums(cfgs) == 1)
      if (length(specific_cfgs) > 0) {
        p_specific[, s] <- safe_rowSums_q(specific_cfgs)
      } else {
        p_specific[, s] <- 0
      }
    }

  } else {
    # shared_het parameterization
    # Het components are columns 2:S
    # Heterogeneous = any het column active in the config

    het_cols <- 2:S

    # Guard: cfgs may need drop=FALSE for single het column
    cfg_het_part <- cfgs[, het_cols, drop = FALSE]
    has_het <- rowSums(cfg_het_part) > 0

    het_cfgs <- which(has_het)
    hom_cfgs <- which(!has_het)

    # Safe rowSums that handles length-0 and length-1 index vectors
    safe_rowSums_q <- function(indices) {
      if (length(indices) == 0) return(rep(0, P))
      if (length(indices) == 1) return(q[, indices])
      rowSums(q[, indices, drop = FALSE])
    }

    p_het_global  <- safe_rowSums_q(het_cfgs)
    p_homogeneous <- safe_rowSums_q(hom_cfgs)
    bf_het <- p_het_global / pmax(p_homogeneous, 1e-300)

    # Per het-component: P(this contrast is active)
    p_het_pair <- matrix(NA_real_, P, length(het_cols),
                          dimnames = list(NULL, sn[het_cols]))
    for (k in seq_along(het_cols)) {
      active_cfgs <- which(cfgs[, het_cols[k]] == 1)
      p_het_pair[, k] <- safe_rowSums_q(active_cfgs)
    }

    # Specificity: het without shared
    p_specific <- matrix(NA_real_, P, S, dimnames = list(NULL, sn))
    for (s in 1:S) {
      specific_cfgs <- which(cfgs[, s] == 1 & rowSums(cfgs) == 1)
      if (length(specific_cfgs) > 0) {
        p_specific[, s] <- safe_rowSums_q(specific_cfgs)
      } else {
        p_specific[, s] <- 0
      }
    }
  }

  list(p_het_global  = p_het_global,
       p_homogeneous = p_homogeneous,
       bf_het        = bf_het,
       p_het_pair    = p_het_pair,
       p_specific    = p_specific)
}


#' Quantitative heterogeneity from subtype z-statistics.
#'
#' Tests whether effect sizes are EQUAL across subtypes,
#' given that multiple subtypes are active. Uses Cochran's Q
#' and pairwise contrasts.
#'
#' @keywords internal
quantitative_heterogeneity <- function(z, R, sn) {

  P <- nrow(z)
  S <- ncol(z)

  # ---- Pairwise contrast tests ----
  # z_s1 - z_s2 ~ N(0, R[s1,s1] + R[s2,s2] - 2*R[s1,s2]) under H0: gamma_s1 = gamma_s2

  pairs <- combn(S, 2, simplify = FALSE)
  pair_names <- sapply(pairs, function(p) paste(sn[p[1]], "vs", sn[p[2]]))

  z_contrast <- matrix(NA_real_, P, length(pairs),
                        dimnames = list(NULL, pair_names))
  p_contrast <- matrix(NA_real_, P, length(pairs),
                        dimnames = list(NULL, pair_names))

  for (k in seq_along(pairs)) {
    s1 <- pairs[[k]][1]
    s2 <- pairs[[k]][2]
    v <- R[s1, s1] + R[s2, s2] - 2 * R[s1, s2]
    v <- max(v, 1e-6)   # guard against numerical issues
    z_contrast[, k] <- (z[, s1] - z[, s2]) / sqrt(v)
    p_contrast[, k] <- 2 * pnorm(-abs(z_contrast[, k]))
  }

  # ---- Omnibus Q-statistic (Cochran-like) ----
  # Q = z' (R^{-1} - R^{-1} 1 1' R^{-1} / (1' R^{-1} 1)) z
  #
  # This tests H0: all gamma_s equal vs H1: at least one differs.
  # Under H0, Q ~ chi-squared(S-1).
  #
  # Equivalent to: fit a fixed-effect meta-analysis, then test residual.

  R_inv <- tryCatch(solve(R), error = function(e) {
    # Regularize if singular
    solve(R + 0.01 * diag(S))
  })

  ones <- rep(1, S)
  denom <- as.numeric(t(ones) %*% R_inv %*% ones)
  # Projection matrix that removes the common mean
  H <- R_inv - (R_inv %*% ones %*% t(ones) %*% R_inv) / denom

  Q <- rep(NA_real_, P)
  p_Q <- rep(NA_real_, P)

  for (i in 1:P) {
    zi <- z[i, ]
    Q[i] <- as.numeric(t(zi) %*% H %*% zi)
    p_Q[i] <- pchisq(Q[i], df = S - 1, lower.tail = FALSE)
  }

  # ---- Fixed-effect meta-analytic mean ----
  # mu_hat = (1' R^{-1} z) / (1' R^{-1} 1)
  # z_meta = mu_hat * sqrt(1' R^{-1} 1)

  z_meta <- rep(NA_real_, P)
  for (i in 1:P) {
    mu_hat <- as.numeric(t(ones) %*% R_inv %*% z[i, ]) / denom
    z_meta[i] <- mu_hat * sqrt(denom)
  }

  list(z_contrast   = z_contrast,
       p_contrast   = p_contrast,
       Q            = Q,
       p_Q          = p_Q,
       df_Q         = S - 1,
       z_meta       = z_meta)
}


#' Directional heterogeneity: opposite effects across subtypes.
#'
#' For each pair, tests whether the effects have opposite signs,
#' accounting for null correlation.
#'
#' @keywords internal
directional_heterogeneity <- function(z, R, sn) {

  P <- nrow(z)
  S <- ncol(z)

  pairs <- combn(S, 2, simplify = FALSE)
  pair_names <- sapply(pairs, function(p) paste(sn[p[1]], "vs", sn[p[2]]))

  # P(opposite signs) using bivariate normal probability
  # Under alternative, if z_s1 and z_s2 have true effects gamma_s1 and gamma_s2
  # with opposite signs, the product z_s1 * z_s2 is expected to be negative.
  #
  # Simple approach: flag if z_s1 * z_s2 < 0 and both |z| > threshold
  #
  # More principled: posterior probability of opposite signs from the
  # model. Here we use the z-values directly since we don't have the
  # full posterior over effect sizes.

  opposite <- matrix(FALSE, P, length(pairs),
                      dimnames = list(NULL, pair_names))
  p_opposite <- matrix(NA_real_, P, length(pairs),
                        dimnames = list(NULL, pair_names))

  for (k in seq_along(pairs)) {
    s1 <- pairs[[k]][1]
    s2 <- pairs[[k]][2]
    r12 <- R[s1, s2]

    # Flag opposite signs with both |z| > 1.96
    opposite[, k] <- (z[, s1] * z[, s2] < 0) &
                      (abs(z[, s1]) > 1.96) &
                      (abs(z[, s2]) > 1.96)

    # Quantify: P(opposite signs) under bivariate normal
    # P(Z1 > 0, Z2 < 0 | rho) + P(Z1 < 0, Z2 > 0 | rho)
    # For observed z-values, compute the contrast p-value
    # conditioning on both being non-null.
    #
    # Use the sum z_s1 + z_s2: under opposite effects, the sum
    # is close to zero while the difference is large.
    # Test: |z_s1 + z_s2| << |z_s1 - z_s2|
    v_sum  <- 2 * (1 + r12)
    v_diff <- 2 * (1 - r12)
    z_sum  <- (z[, s1] + z[, s2]) / sqrt(max(v_sum, 1e-6))
    z_diff <- (z[, s1] - z[, s2]) / sqrt(max(v_diff, 1e-6))

    # Ratio: if truly opposite, |z_diff| >> |z_sum|
    # P-value for sum being null (conditioned on diff being large)
    p_sum_null <- 2 * pnorm(-abs(z_sum))
    p_diff_sig <- 2 * pnorm(-abs(z_diff))

    # Joint evidence for opposite effects:
    # difference is significant AND sum is not
    p_opposite[, k] <- p_diff_sig  # conservative: just the contrast p-value
  }

  list(opposite     = opposite,
       p_opposite   = p_opposite,
       n_opposite   = colSums(opposite))
}


#' Contrast-specific heterogeneity for shared+het parameterization.
#'
#' Reports, for each contrast, the posterior probability that it is
#' active, the Bayes factor, and which original subtypes it separates.
#'
#' @keywords internal
contrast_heterogeneity <- function(fit, der, cfgs, sn,
                                    contrast_matrix, sn_orig) {

  q <- extract_fit_q(fit)
  P <- nrow(q)
  S <- ncol(cfgs)

  # Het columns are 2:S in the shared+het parameterization
  het_cols <- 2:S
  n_contrasts <- length(het_cols)

  # Per-contrast PIP (already available from der$pip)
  pip_het <- der$pip[, het_cols, drop = FALSE]

  # Per-contrast Bayes factor
  bf_het <- matrix(NA_real_, P, n_contrasts,
                    dimnames = list(NULL, sn[het_cols]))
  for (k in 1:n_contrasts) {
    col <- het_cols[k]
    active_cfgs   <- which(cfgs[, col] == 1)
    inactive_cfgs <- which(cfgs[, col] == 0)
    p_active   <- rowSums(q[, active_cfgs, drop = FALSE])
    p_inactive <- rowSums(q[, inactive_cfgs, drop = FALSE])
    bf_het[, k] <- p_active / pmax(p_inactive, 1e-300)
  }

  # Map each contrast back to original subtype comparison
  contrast_interpretation <- character(n_contrasts)
  for (k in 1:n_contrasts) {
    c_k <- contrast_matrix[k, ]
    pos <- sn_orig[c_k > 0]
    neg <- sn_orig[c_k < 0]
    contrast_interpretation[k] <- paste0(
      paste(pos, collapse = "+"), " vs ",
      paste(neg, collapse = "+"))
  }

  list(pip_het     = pip_het,
       bf_het      = bf_het,
       interpretation = contrast_interpretation)
}


#' Combine all heterogeneity measures into a single summary table.
#'
#' @keywords internal
combine_heterogeneity <- function(results, P, sn_orig, der) {

  protein_ids <- derive_protein_ids(results, der, P)

  # Build summary data.frame
  summary_df <- data.frame(
    gene  = protein_ids,

    # From model: posterior probability of heterogeneity
    p_het_posterior = results$posterior$p_het_global,
    bf_het          = results$posterior$bf_het,

    # From Cochran Q: frequentist heterogeneity test
    Q_stat          = results$quantitative$Q,
    Q_pval          = results$quantitative$p_Q,

    # Fixed-effect meta z (for reference)
    z_meta          = results$quantitative$z_meta,

    stringsAsFactors = FALSE
  )

  # Add pairwise contrast p-values
  if (ncol(results$quantitative$p_contrast) > 0) {
    pairwise_df <- as.data.frame(results$quantitative$p_contrast)
    colnames(pairwise_df) <- paste0("p_", colnames(pairwise_df))
    summary_df <- cbind(summary_df, pairwise_df)
  }

  # Add opposite-effects flags
  if (!is.null(results$directional$opposite)) {
    opp_df <- as.data.frame(results$directional$opposite)
    colnames(opp_df) <- paste0("opposite_", colnames(opp_df))
    summary_df <- cbind(summary_df, opp_df)
  }

  # Add per-subtype specificity probabilities
  if (!is.null(results$posterior$p_specific)) {
    spec_df <- as.data.frame(results$posterior$p_specific)
    colnames(spec_df) <- paste0("p_specific_", colnames(spec_df))
    summary_df <- cbind(summary_df, spec_df)
  }

  # Classification
  summary_df$het_class <- classify_heterogeneity(
    p_het  = results$posterior$p_het_global,
    Q_pval = results$quantitative$p_Q,
    bfdr   = der$bfdr
  )

  summary_df
}


#' Extract posterior configuration probabilities from a fit object.
#' @keywords internal
extract_fit_q <- function(fit) {

  q <- fit$q %||% fit$responsibilities
  if (is.null(q))
    stop("fit must contain posterior configuration probabilities in 'q'.")

  as.matrix(q)
}


#' Derive stable protein identifiers for heterogeneity summaries.
#' @keywords internal
derive_protein_ids <- function(results, der, P) {

  ids <- NULL

  if (!is.null(der$proteins))
    ids <- der$proteins
  if (is.null(ids) && !is.null(names(der$bfdr)))
    ids <- names(der$bfdr)
  if (is.null(ids) && !is.null(rownames(results$posterior$p_het_pair)))
    ids <- rownames(results$posterior$p_het_pair)

  if (is.null(ids))
    ids <- paste0("gene_", seq_len(P))

  as.character(ids)
}


#' Classify genes into heterogeneity categories.
#'
#' @keywords internal
classify_heterogeneity <- function(p_het, Q_pval, bfdr,
                                    het_threshold = 0.8,
                                    q_threshold   = 0.05,
                                    bfdr_threshold = 0.05) {

  P <- length(p_het)
  classes <- character(P)

  for (i in 1:P) {
    if (bfdr[i] > bfdr_threshold) {
      classes[i] <- "null"
    } else if (p_het[i] > het_threshold && Q_pval[i] < q_threshold) {
      classes[i] <- "strong_heterogeneity"
    } else if (p_het[i] > het_threshold || Q_pval[i] < q_threshold) {
      classes[i] <- "moderate_heterogeneity"
    } else {
      classes[i] <- "homogeneous"
    }
  }

  classes
}


###############################################################################
# Convenience wrappers
###############################################################################


#' Quick heterogeneity test for the direct subtype parameterization.
#'
#' @param z_matrix  P x S subtype z-matrix (same as used for inference)
#' @param R_null    S x S null correlation
#' @param fit       model fit
#' @param der       derived results
#' @param cfgs      configurations
#' @param sn        subtype names
#' @return heterogeneity results
#' @export
test_het_direct <- function(z_matrix, R_null, fit, der, cfgs, sn) {
  test_heterogeneity(
    z_original       = z_matrix,
    R_original       = R_null,
    fit              = fit,
    der              = der,
    cfgs             = cfgs,
    sn               = sn,
    parameterization = "direct"
  )
}


#' Quick heterogeneity test for the shared+het parameterization.
#'
#' @param z_original   P x S_orig ORIGINAL subtype z-statistics
#' @param R_original   S_orig x S_orig null correlation among original subtypes
#' @param z_inference  P x S z-matrix used for inference (shared + hets)
#' @param R_inference  S x S null correlation of inference z-matrix
#' @param fit          model fit
#' @param der          derived results
#' @param cfgs         configurations
#' @param sn           column names of z_inference
#' @param C            contrast matrix used to build z_inference
#' @return heterogeneity results
#' @export
test_het_shared <- function(z_original, R_original,
                             z_inference, R_inference,
                             fit, der, cfgs, sn, C) {
  test_heterogeneity(
    z_original       = z_original,
    R_original       = R_original,
    fit              = fit,
    der              = der,
    cfgs             = cfgs,
    sn               = sn,
    parameterization = "shared_het",
    contrast_matrix  = C
  )
}


###############################################################################
# USAGE EXAMPLES
###############################################################################
#
# ---- Example 1: Direct parameterization (RCC, S=2) ----
#
#   z_matrix <- cbind(clear_cell = z_cc, papillary = z_pap)
#   R <- estimate_null_correlation(z_matrix, NULL, ctrl)
#   cfgs <- build_configs(2, colnames(z_matrix))
#   fit  <- fit_vb(z_matrix, R, cfgs, ctrl)
#   der  <- compute_derived(fit, z_matrix, R, cfgs, colnames(z_matrix), ctrl)
#
#   het <- test_het_direct(z_matrix, R, fit, der, cfgs, colnames(z_matrix))
#
#   # Key outputs:
#   het$summary$het_class           # "null", "homogeneous", "moderate_het", "strong_het"
#   het$summary$p_het_posterior      # P(qualitative heterogeneity)
#   het$summary$Q_pval              # Cochran Q p-value
#   het$summary$p_clear_cell.vs.papillary  # pairwise contrast p-value
#   het$posterior$p_specific         # P(only cc active), P(only pap active)
#
# ---- Example 2: Shared+het parameterization (Lung cancer, S=3) ----
#
#   z_sub <- cbind(LUAD = z_luad, SQC = z_sqc, SCLC = z_sclc)
#   R_sub <- cor(z_sub)   # null correlation among subtypes
#
#   # Build shared+het z-matrix
#   sh <- construct_shared_het_z(z_lc, z_sub)
#
#   # Run inference
#   cfgs <- build_configs(3, colnames(sh$z_matrix))
#   fit  <- fit_vb(sh$z_matrix, sh$R_null, cfgs, ctrl)
#   der  <- compute_derived(fit, sh$z_matrix, sh$R_null, cfgs,
#                           colnames(sh$z_matrix), ctrl)
#
#   # Heterogeneity test
#   het <- test_het_shared(
#     z_original  = z_sub,
#     R_original  = R_sub,
#     z_inference = sh$z_matrix,
#     R_inference = sh$R_null,
#     fit = fit, der = der, cfgs = cfgs,
#     sn = colnames(sh$z_matrix),
#     C  = sh$C
#   )
#
#   # Key outputs:
#   het$summary$het_class                    # classification
#   het$summary$p_het_posterior               # P(any contrast active)
#   het$contrast$pip_het                      # PIP per contrast
#   het$contrast$interpretation               # "LUAD vs SQC", "SCLC vs LUAD+SQC"
#   het$summary$Q_pval                        # Cochran Q across all subtypes
#   het$directional$opposite                  # opposite-effect flags
#
# ---- Example 3: Visualize heterogeneity ----
#
#   # Volcano-like plot: meta z (shared signal) vs Q (heterogeneity)
#   plot(het$summary$z_meta, -log10(het$summary$Q_pval),
#        xlab = "Meta-analytic z (shared effect)",
#        ylab = "-log10(Q p-value) [heterogeneity]",
#        col  = ifelse(het$summary$het_class == "strong_heterogeneity",
#                      "red", "grey60"),
#        pch  = 16, cex = 0.6)
#   abline(h = -log10(0.05), lty = 2)
#   # Upper-left: heterogeneous without shared signal (opposite effects)
#   # Upper-right: heterogeneous with shared signal (different magnitudes)
#   # Lower-right: homogeneous shared signal
#   # Lower-left: null
#
###############################################################################
