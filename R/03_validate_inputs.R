###############################################################################
# MS-BPWAS  -  03_validate_inputs.R
# Read / validate all inputs; compute the residual subtype.
###############################################################################

#' Master input validator called once at the start of ms_bpwas().
#' @return list(gwas, pmodels, subtype_names, S, prevalences)
#' @keywords internal
validate_inputs <- function(gwas_list, protein_models, plink_prefix, ctrl) {

 ## ---- 1. GWAS summary statistics --------------------------------
 stopifnot(is.list(gwas_list), !is.null(names(gwas_list)))
 if (!"mixture" %in% names(gwas_list))
   stop("gwas_list must contain an element named 'mixture'.")
 explicit <- setdiff(names(gwas_list), "mixture")
 if (length(explicit) < 1L)
   stop("Provide at least one named subtype besides 'mixture'.")

 req <- c("SNP","A1","A2","BETA","SE","P","N")
 gwas <- list()
 for (nm in names(gwas_list)) {
   gwas[[nm]] <- read_gwas(gwas_list[[nm]], nm, ctrl$verbose)
   miss <- setdiff(req, colnames(gwas[[nm]]))
   if (length(miss)) stop(sprintf("'%s' missing cols: %s", nm, toString(miss)))
   if (!"Z" %in% colnames(gwas[[nm]]))
     gwas[[nm]]$Z <- gwas[[nm]]$BETA / gwas[[nm]]$SE
 }

 ## ---- 2. Resolve subtype prevalences ----------------------------
 ##
 ## subtype_prevalences: named numeric vector giving the fraction
 ## of mixture-GWAS cases belonging to each explicit subtype.
 ## Must NOT include "residual" - that is 1 - sum(pi_s).
 ##
 ## If NULL, we fall back to GWAS sample-size ratios as a proxy.
 prev <- ctrl$subtype_prevalences
 if (!is.null(prev)) {
   stopifnot(is.numeric(prev), length(prev) == length(explicit))
   if (is.null(names(prev))) names(prev) <- explicit
   miss <- setdiff(explicit, names(prev))
   if (length(miss)) stop("subtype_prevalences missing: ", toString(miss))
   pi_sum <- sum(prev)
   if (pi_sum >= 1 || pi_sum <= 0)
     stop("subtype_prevalences must be positive and sum to < 1 ",
          "(residual = 1 - sum).  Got sum = ", pi_sum)
   pi_res <- 1 - pi_sum
   if (ctrl$verbose)
     msg("  Prevalences (user): %s,  residual=%.3f",
         paste(sprintf("%s=%.3f", names(prev), prev), collapse=", "),
         pi_res)
 } else {
   ## Fallback: N_s / N_mix as proxy
   N_mix <- median(gwas[["mixture"]]$N, na.rm = TRUE)
   prev <- setNames(numeric(length(explicit)), explicit)
   for (nm in explicit)
     prev[nm] <- median(gwas[[nm]]$N, na.rm = TRUE) / N_mix
   pi_sum <- sum(prev)
   if (pi_sum >= 1) { prev <- prev / (pi_sum + 0.05); pi_sum <- sum(prev) }
   pi_res <- 1 - pi_sum
   if (ctrl$verbose)
     msg("  Prevalences (from N): %s,  residual=%.3f",
         paste(sprintf("%s=%.3f", names(prev), prev), collapse=", "),
         pi_res)
 }

 ## ---- 3. Residual subtype ---------------------------------------
 if (ctrl$verbose) msg("  Computing residual subtype ...")
 gwas[["residual"]] <- compute_residual(
   gwas[["mixture"]], gwas[explicit], prev, ctrl$verbose)
 subtype_names <- c(explicit, "residual")

 ## ---- 4. Protein models -----------------------------------------
 pmodels <- load_protein_models(protein_models, ctrl)

 ## ---- 5. PLINK files --------------------------------------------
 for (ext in c(".bed",".bim",".fam")) {
   f <- paste0(plink_prefix, ext)
   if (!file.exists(f)) stop("Missing: ", f)
 }

 list(gwas = gwas, pmodels = pmodels,
      subtype_names = subtype_names, S = length(subtype_names),
      prevalences = c(prev, residual = pi_res))
}


# ---------- helpers ---------------------------------------------------

#' Read a single GWAS (path -> data.frame, or pass through)
#' @keywords internal
read_gwas <- function(x, label, verbose = TRUE) {
 if (is.data.frame(x)) return(x)
 stopifnot(is.character(x))
 if (verbose) msg("  Reading %s: %s", label, x)
 if (grepl("\\.rds$",x,TRUE))  return(readRDS(x))
 if (grepl("\\.gz$",x))        return(read.table(gzfile(x), header=TRUE,
                                                  stringsAsFactors=FALSE))
 read.table(x, header = TRUE, stringsAsFactors = FALSE)
}


#' Load and QC protein prediction models from .RData or list
#' @keywords internal
load_protein_models <- function(protein_models, ctrl) {
 if (is.character(protein_models) && length(protein_models) == 1L) {
   if (ctrl$verbose) msg("  Loading protein models: %s", protein_models)
   env <- new.env(parent = emptyenv())
   load(protein_models, envir = env)
   pm <- NULL
   for (on in ls(env)) {
     cand <- get(on, envir = env)
     if (is.list(cand) && length(cand) > 1) { pm <- cand; break }
   }
   if (is.null(pm)) stop("No suitable list found in ", protein_models)
 } else if (is.list(protein_models)) {
   pm <- protein_models
 } else stop("protein_models: expected file path or list.")

 n0 <- length(pm)
 pm <- Filter(function(m) {
   ok_r2  <- !is.null(m$cv_R2)  && m$cv_R2 >= ctrl$min_cv_R2
   ok_snp <- !is.null(m$weights) && nrow(m$weights) >= ctrl$min_n_snps
   ok_r2 && ok_snp
 }, pm)
 if (ctrl$verbose) msg("  Protein models: %d / %d pass QC", length(pm), n0)
 pm
}


#' Compute the residual subtype using prevalence-based decomposition.
#'
#' The mixture GWAS estimates a prevalence-weighted average:
#'   beta_mix approx Sigma_s pi_s beta_s  +  pi_res beta_res
#'
#' We solve for the residual:
#'   beta_res = (beta_mix - Sigma_s pi_s beta_s) / pi_res
#'
#' SE propagated under independence of subtype GWAS.
#'
#' @param gwas_mix   data.frame for the mixture (overall) GWAS
#' @param gwas_subs  named list of subtype GWAS data.frames
#' @param prev       named numeric, prevalence of each explicit subtype
#' @param verbose    logical
#' @return data.frame with residual GWAS columns
#' @keywords internal
compute_residual <- function(gwas_mix, gwas_subs, prev, verbose = TRUE) {

 out <- gwas_mix[, c("SNP","A1","A2")]
 pi_res <- 1 - sum(prev)
 if (pi_res < 0.01) {
   warning("Residual prevalence < 1%.  Residual estimates will be noisy.")
   pi_res <- max(pi_res, 0.01)
 }

 wt_beta     <- rep(0, nrow(out))
 var_contrib <- rep(0, nrow(out))

 for (nm in names(gwas_subs)) {
   gs   <- gwas_subs[[nm]]
   pi_s <- prev[nm]
   idx  <- match(out$SNP, gs$SNP)
   ok   <- !is.na(idx)
   if (!any(ok)) { warning("No SNP overlap with ", nm); next }

   beta_s <- se_s <- rep(NA_real_, nrow(out))
   beta_s[ok] <- gs$BETA[idx[ok]]
   se_s[ok]   <- gs$SE[idx[ok]]

   # flip if alleles swapped
   flip <- ok & toupper(out$A1) == toupper(gs$A2[idx]) &
                toupper(out$A2) == toupper(gs$A1[idx])
   beta_s[flip] <- -beta_s[flip]

   wt_beta     <- wt_beta     + ifelse(ok, pi_s * beta_s, 0)
   var_contrib <- var_contrib  + ifelse(ok, pi_s^2 * se_s^2, 0)

   out         <- out[ok, ]
   wt_beta     <- wt_beta[ok]
   var_contrib <- var_contrib[ok]
 }

 beta_mix <- gwas_mix$BETA[match(out$SNP, gwas_mix$SNP)]
 se_mix   <- gwas_mix$SE[match(out$SNP, gwas_mix$SNP)]

 beta_res <- (beta_mix - wt_beta) / pi_res
 se_res   <- sqrt(se_mix^2 / pi_res^2 + var_contrib / pi_res^2)

 res <- data.frame(
   SNP = out$SNP, A1 = out$A1, A2 = out$A2,
   BETA = beta_res, SE = se_res,
   P = 2 * pnorm(-abs(beta_res / se_res)),
   N = median(gwas_mix$N, na.rm = TRUE) * pi_res,
   Z = beta_res / se_res,
   stringsAsFactors = FALSE)

 if (verbose) msg("    residual: pi_res=%.3f, N_eff=%.0f, %d SNPs",
                   pi_res, res$N[1], nrow(res))
 res
}
