###############################################################################
# MS-BPWAS  -  12_ms_bpwas.R
# Main entry point.  Sources nothing (all functions assumed loaded).
###############################################################################

#' @title Multi-Subtype Bayesian Proteome-Wide Association Study
#'
#' @description Jointly models protein-disease associations across S
#'   subtypes using pretrained protein prediction weights and GWAS
#'   summary statistics.  Returns posterior probabilities over 2^S
#'   association configurations, marginal PIPs, conditional p-values,
#'   posterior effect sizes, and Bayesian FDR.
#'
#' @param gwas_list Named list: \code{"mixture"}, \code{"subtype1"},
#'   \code{"subtype2"}, ....
#'   Each element is a file path or data.frame with columns
#'   \code{SNP, A1, A2, BETA, SE, P, N}.
#'   The residual subtype is computed automatically.
#' @param protein_models Path to \code{.RData} file or a list.
#'   Each element must have \code{protein_id, gene, chr,
#'   weights} (data.frame: \code{SNP, A1, weight}),
#'   \code{cv_R2}.
#' @param plink_prefix PLINK binary-file prefix (\code{.bed/.bim/.fam}).
#' @param control Output of \code{\link{ms_bpwas_control}()}.
#'
#' @return Object of class \code{"msbpwas"} with components:
#' \describe{
#'   \item{posterior_configs}{P x 2^S matrix of configuration posteriors}
#'   \item{pip}{P x S posterior inclusion probabilities per subtype}
#'   \item{p_null}{P-vector of posterior null probabilities}
#'   \item{map_config}{P-vector of MAP configuration labels}
#'   \item{p_conditional}{P-vector of conditional chi-squared p-values}
#'   \item{effect_sizes}{P x S posterior mean effects (model-averaged)}
#'   \item{effect_ses}{P x S posterior SDs}
#'   \item{bfdr}{P-vector of Bayesian FDR values}
#'   \item{hyperparams}{Estimated global hyperparameters}
#'   \item{z_matrix}{P x S PWAS Z-statistics}
#'   \item{R_overlap}{S x S null correlation matrix}
#'   \item{convergence}{Convergence diagnostics}
#' }
#'
#' @export
ms_bpwas <- function(gwas_list,
                     protein_models,
                     plink_prefix,
                     control = ms_bpwas_control()) {

 cl <- match.call(); t0 <- proc.time()
 V <- control$verbose

 if(V) {
   msg("============================================================")
   msg("  MS-BPWAS: Multi-Subtype Bayesian PWAS")
   msg("============================================================")
 }

 ## 0. Validate inputs
 if(V) msg("\n--- Stage 0: input validation ---")
 inp <- validate_inputs(gwas_list, protein_models, plink_prefix, control)
 S <- inp$S; sn <- inp$subtype_names
 if(V) msg("  S=%d: %s | %d protein models", S, toString(sn), length(inp$pmodels))

 ## 1. PWAS Z-statistics
 if(V) msg("\n--- Stage 1: PWAS Z-statistics ---")
 pw  <- compute_all_pwas_z(inp$gwas, inp$pmodels, plink_prefix, sn, control)
 zmat <- pw$z_matrix; P <- nrow(zmat)
 if(V) msg("  %d proteins x %d subtypes", P, S)

 ## 2. Null correlation
 if(V) msg("\n--- Stage 2: null correlation ---")
 R <- estimate_null_correlation(zmat, inp$gwas, control)
 if(V) { msg("  R:"); print(round(R,4)) }

 ## 3. Configuration space
 cfgs <- build_configs(S, sn)
 if(V) msg("\n--- %d configurations ---", nrow(cfgs))

 ## 4. Inference
 if(V) msg("\n--- Stage 3: %s inference (prior=%s) ---",
            toupper(control$method), control$prior_type)
 fit <- switch(control$method,
               vb   = fit_vb(zmat, R, cfgs, control),
               em   = fit_em(zmat, R, cfgs, control),
               mcmc = fit_mcmc(zmat, R, cfgs, control),
               stop("Unknown method"))

 ## 5. Derived quantities
 if(V) msg("\n--- Stage 4: derived quantities ---")
 der <- compute_derived(fit, zmat, R, cfgs, sn, control)

 ## 6. Assemble
 elapsed <- (proc.time()-t0)[3]
 out <- structure(list(
   posterior_configs = fit$q,
   pip              = der$pip,
   p_null           = der$p_null,
   map_config       = der$map_config,
   p_conditional    = der$p_cond,
   effect_sizes     = fit$post_mean,
   effect_ses       = sqrt(pmax(fit$post_var_diag, 0)),
   bfdr             = der$bfdr,
   z_matrix         = zmat,
   R_overlap        = R,
   configs          = cfgs,
   hyperparams      = fit$hyperparams,
   convergence      = fit$convergence,
   protein_info     = pw$protein_info,
   protein_ids      = rownames(zmat),
   subtype_names    = sn,
   S = S, P = P,
   expected_n_active = der$expected_n_active,
   elapsed_sec      = elapsed,
   control          = control,
   call             = cl
 ), class = "msbpwas")

 if(V) {
   msg("\n  Done in %.1f s. BFDR<%.2f: %d discoveries.",
       elapsed, control$bfdr_threshold, sum(der$bfdr<control$bfdr_threshold))
 }
 out
}
