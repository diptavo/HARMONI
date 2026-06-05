###############################################################################
# MS-BPWAS  -  02_control.R
# Control-parameter constructor.  One function, no side-effects.
###############################################################################

#' Build a control-parameter list for \code{ms_bpwas()}
#'
#' @param method       Inference engine: \code{"vb"} (default),
#'                     \code{"em"}, or \code{"mcmc"}.
#' @param prior_type   Effect-size prior on active mu:
#'                     \code{"scale_mixture"} (default, adaptive shrinkage),
#'                     \code{"horseshoe"}, or \code{"soft_spike_slab"}.
#' @param sigma2_grid  Numeric vector of variance components for the scale
#'                     mixture.  \code{NULL} -> auto-generate.
#' @param n_grid       Grid size when auto-generating (default 8).
#' @param sigma2_small Small-component variance for soft spike-slab.
#' @param sigma2_large Large-component variance for soft spike-slab.
#' @param pi0_init     Starting null probability (default 0.80).
#' @param xi_init      Starting per-subtype inclusion probability (default 0.20).
#' @param subtype_prevalences Named numeric vector giving the fraction of
#'   mixture-GWAS cases belonging to each explicit subtype (e.g.,
#'   \code{c(subtype1 = 0.10, subtype2 = 0.50)}).  Must sum to less than 1;
#'   the residual subtype prevalence is computed as \code{1 - sum(pi)}.
#'   Determines how the residual-subtype GWAS is derived from the mixture.
#'   If \code{NULL}, prevalences are approximated from GWAS sample sizes.
#' @param overlap_correction How to handle sample overlap between GWAS:
#'   \code{"ldsc"}, \code{"known"}, \code{"empirical"}, \code{"none"}.
#' @param N_overlap_matrix  S x S matrix of shared sample counts
#'   (required when \code{overlap_correction = "known"}).
#' @param max_iter     Max VB / EM iterations.
#' @param tol          Relative convergence tolerance.
#' @param n_chains     MCMC chains.
#' @param n_iter       MCMC iterations per chain.
#' @param n_warmup     MCMC warm-up iterations.
#' @param thin         MCMC thinning interval.
#' @param min_cv_R2    Minimum protein-model R^2 to keep.
#' @param min_n_snps   Minimum SNPs in prediction model.
#' @param ld_window    cis window (bp) for LD computation.
#' @param ld_shrinkage Apply Ledoit-Wolf to LD matrix?
#' @param bfdr_threshold Discovery threshold.
#' @param verbose      Print progress?
#' @param n_cores      Parallel cores (future use).
#' @return A list of class \code{"msbpwas_control"}.
#' @export
ms_bpwas_control <- function(
   method             = c("vb", "em", "mcmc"),
   prior_type         = c("scale_mixture", "horseshoe", "soft_spike_slab"),
   sigma2_grid        = NULL,
   n_grid             = 10L,
   sigma2_small       = 0.1,
   sigma2_large       = 9.0,
   sigma2_null        = 0.005^2,
   pi0_init           = 0.80,
   xi_init            = 0.20,
   pi0_prior_a        = 10.0,
   pi0_prior_b        = 2.0,
   subtype_prevalences = NULL,
   overlap_correction = c("ldsc", "known", "empirical", "none"),
   N_overlap_matrix   = NULL,
   max_iter           = 200L,
   tol                = 1e-6,
   n_chains           = 4L,
   n_iter             = 5000L,
   n_warmup           = 2000L,
   thin               = 5L,
   min_cv_R2          = 0.01,
   min_n_snps         = 2L,
   ld_window          = 1e6,
   ld_shrinkage       = TRUE,
   bfdr_threshold     = 0.05,
   verbose            = TRUE,
   n_cores            = 1L
) {
 method             <- match.arg(method)
 prior_type         <- match.arg(prior_type)
 overlap_correction <- match.arg(overlap_correction)

 ## Variance grid for scale mixture: covers z-score effect range
 ## PWAS Z-statistics have effects in z-score units (typically 0.5-5)
 ## Grid spans sigma from 0.1 to 4.0 (sigma^2 from 0.01 to 16)
 if (is.null(sigma2_grid) && prior_type == "scale_mixture")
   sigma2_grid <- exp(seq(log(0.01), log(16.0), length.out = n_grid))

 structure(list(
   method = method, prior_type = prior_type,
   sigma2_grid = sigma2_grid, n_grid = n_grid,
   sigma2_small = sigma2_small, sigma2_large = sigma2_large,
   sigma2_null = sigma2_null,
   pi0_init = pi0_init, xi_init = xi_init,
   pi0_prior_a = pi0_prior_a, pi0_prior_b = pi0_prior_b,
   subtype_prevalences = subtype_prevalences,
   overlap_correction = overlap_correction,
   N_overlap_matrix = N_overlap_matrix,
   max_iter = max_iter, tol = tol,
   n_chains = n_chains, n_iter = n_iter,
   n_warmup = n_warmup, thin = thin,
   min_cv_R2 = min_cv_R2, min_n_snps = min_n_snps,
   ld_window = ld_window, ld_shrinkage = ld_shrinkage,
   bfdr_threshold = bfdr_threshold,
   verbose = verbose, n_cores = n_cores
 ), class = "msbpwas_control")
}
