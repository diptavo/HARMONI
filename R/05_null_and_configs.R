###############################################################################
# MS-BPWAS  -  05_null_and_configs.R
# (a) Estimate the SxS null correlation R
# (b) Enumerate the 2^S configuration space
# (c) Compute the configuration prior
###############################################################################

# ===========================================================================
# A.  NULL CORRELATION MATRIX
# ===========================================================================

#' Estimate R so that z_p ~ N(0, R) under the global null.
#' @param z_matrix  P x S
#' @param gwas      named list of GWAS data.frames (unused for "empirical"/"ldsc")
#' @param ctrl      control list
#' @return S x S positive-definite matrix
#' @keywords internal
estimate_null_correlation <- function(z_matrix, gwas, ctrl) {

 S <- ncol(z_matrix); sn <- colnames(z_matrix)
 R <- diag(S); dimnames(R) <- list(sn, sn)
 method <- ctrl$overlap_correction

 # ---- none ----
 if (method == "none") {
   if (ctrl$verbose) msg("  Overlap: NONE (identity)")
   return(R)
 }

 # ---- known ----
 if (method == "known") {
   if (ctrl$verbose) msg("  Overlap: KNOWN matrix")
   Nm <- ctrl$N_overlap_matrix
   stopifnot(!is.null(Nm), nrow(Nm)==S, ncol(Nm)==S)
   dg <- diag(Nm)
   for (i in 1:(S-1)) for (j in (i+1):S) {
     R[i,j] <- Nm[i,j] / sqrt(dg[i]*dg[j]); R[j,i] <- R[i,j]
   }
   return(make_psd(R))
 }

 # ---- empirical ----
 if (method == "empirical") {
   if (ctrl$verbose) msg("  Overlap: EMPIRICAL null")
   thr <- 1.5
   for (t in c(1.5, 2.0, 2.5)) {
     nmask <- apply(abs(z_matrix), 1, max, na.rm=TRUE) < t
     if (sum(nmask) >= max(3*S, 30)) { thr <- t; break }
   }
   Z0 <- z_matrix[nmask, , drop=FALSE]
   if (ctrl$verbose) msg("    %d null proteins (|z|<%g)", nrow(Z0), thr)
   if (nrow(Z0) < S+2) { warning("Too few nulls -> identity R"); return(R) }
   Re <- cor(Z0, use="pairwise.complete.obs")
   lam <- min(0.5, S / nrow(Z0))
   R <- (1-lam)*Re + lam*diag(S); dimnames(R) <- list(sn,sn)
   return(make_psd(R))
 }

 # ---- ldsc intercept approximation ----
 if (method == "ldsc") {
   if (ctrl$verbose) msg("  Overlap: LDSC intercept approx")
   nmask <- apply(abs(z_matrix), 1, max, na.rm=TRUE) < 2.5
   Z0 <- z_matrix[nmask, , drop=FALSE]
   if (ctrl$verbose) msg("    %d proteins for intercept", nrow(Z0))
   if (nrow(Z0) < 30) { warning("Too few -> identity R"); return(R) }
   for (i in 1:(S-1)) for (j in (i+1):S) {
     R[i,j] <- mean(Z0[,i]*Z0[,j], na.rm=TRUE); R[j,i] <- R[i,j]
   }
   return(make_psd(R))
 }

 stop("Unknown overlap method: ", method)
}


# ===========================================================================
# B.  CONFIGURATION ENUMERATION
# ===========================================================================

#' Binary matrix of all 2^S configurations (null = row 1).
#' @keywords internal
build_configs <- function(S, subtype_names) {
 nc <- 2^S
 if (nc > 1e6) stop("2^S = ", nc, " too large; reduce S or prune.")
 cfgs <- matrix(0L, nc, S, dimnames = list(NULL, subtype_names))
 for (i in 0:(nc-1)) for (s in 1:S)
   cfgs[i+1L, s] <- bitwAnd(bitwShiftR(i, s-1L), 1L)
 labs <- apply(cfgs, 1, function(r) {
   a <- subtype_names[r==1]; if(!length(a)) "null" else paste(a, collapse="+")
 })
 rownames(cfgs) <- labs
 attr(cfgs, "size") <- as.integer(rowSums(cfgs))
 cfgs
}


# ===========================================================================
# C.  CONFIGURATION PRIOR
# ===========================================================================

#' P(c) = pi0*I(c=null) + (1-pi0)*prod_s xi_s^{c_s}(1-xi_s)^{1-c_s} (normalised)
#' @return numeric vector length 2^S summing to 1
#' @keywords internal
config_prior <- function(cfgs, pi0, xi) {
 S <- ncol(cfgs); nc <- nrow(cfgs)
 if (length(xi) == 1) xi <- rep(xi, S)
 if (length(pi0) == 0 || is.null(pi0)) pi0 <- 0.80  # safety default
 pi0 <- max(min(pi0, 0.999), 0.001)
 cs <- rowSums(cfgs)

 lp <- rep(-Inf, nc)
 for (i in which(cs > 0))
   lp[i] <- sum(ifelse(cfgs[i,]==1, log(xi), log(1-xi)))

 nn <- cs > 0
 lp[nn] <- lp[nn] - logsumexp(lp[nn])   # normalise non-null part

 out <- numeric(nc)
 out[cs == 0] <- pi0
 out[nn]      <- (1 - pi0) * exp(lp[nn])
 out / sum(out)
}
