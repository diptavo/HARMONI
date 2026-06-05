###############################################################################
# MS-BPWAS  -  10_derived.R
# Derived quantities computed once after inference.
###############################################################################

#' Compute PIPs, MAP configs, conditional p-values, BFDR
#' @keywords internal
compute_derived <- function(fit, z_matrix, R, cfgs, sn, ctrl) {

 P<-nrow(z_matrix); S<-ncol(z_matrix); nc<-nrow(cfgs); q<-fit$q
 Ri <- solve(R)

 ## 1. P(null)
 p_null <- q[,1]; names(p_null) <- rownames(z_matrix)

 ## 2. Marginal PIP per subtype
 pip <- matrix(0, P, S, dimnames=dimnames(z_matrix))
 for (s in seq_len(S))
   pip[,s] <- rowSums(q[, which(cfgs[,s]==1), drop=FALSE])

 ## 3. MAP configuration
 mi <- apply(q, 1, which.max)
 map_config <- rownames(cfgs)[mi]; names(map_config) <- rownames(z_matrix)

 ## 4. Conditional p-value  (chi-squared projection test)
 ##    T = z' R^{-1} A (A'R^{-1}A)^{-1} A' R^{-1} z  ~  chi2(|c|)
 p_cond <- rep(NA_real_, P); names(p_cond) <- rownames(z_matrix)
 for (p in seq_len(P)) {
   act <- which(cfgs[mi[p],]==1); k <- length(act)
   if (k==0) { p_cond[p] <- 1; next }
   zp <- z_matrix[p,]; if(anyNA(zp)) next
   Ac <- matrix(0,S,k); for(j in seq_len(k)) Ac[act[j],j] <- 1
   AtRiz <- as.numeric(crossprod(Ac, Ri %*% zp))
   AtRiA <- crossprod(Ac, Ri %*% Ac)
   Tstat <- as.numeric(crossprod(AtRiz, solve_safe(AtRiA) %*% AtRiz))
   p_cond[p] <- pchisq(Tstat, df=k, lower.tail=FALSE)
 }

 ## 5. Bayesian FDR
 bfdr <- bfdr_from_pnull(p_null)

 ## 6. expected number of active subtypes
 en <- rowSums(pip); names(en) <- rownames(z_matrix)

 ## print top hits
 if (ctrl$verbose) {
   msg("  BFDR<0.05: %d  |  P(null)<0.10: %d",
       sum(bfdr<0.05), sum(p_null<0.10))
   top <- order(p_null)[1:min(10,P)]
   msg("\n  Top proteins:")
   msg("  %-15s P(null)  MAP%-28s  p_cond     PIPs", "Protein", "")
   for(i in top)
     msg("  %-15s %.4f   %-28s  %.2e  %s",
         rownames(z_matrix)[i], p_null[i], map_config[i], p_cond[i],
         paste(sprintf("%.2f",pip[i,]), collapse=" "))
 }

 list(pip=pip, p_null=p_null, map_config=map_config,
      p_cond=p_cond, bfdr=bfdr, expected_n_active=en)
}


#' Bayesian FDR from posterior null probabilities
#' @keywords internal
bfdr_from_pnull <- function(pn) {
 P <- length(pn); ord <- order(pn)
 ca <- cumsum(pn[ord]) / seq_len(P)
 # enforce monotonicity from bottom
 for(k in (P-1):1) ca[k] <- min(ca[k], ca[k+1])
 out <- numeric(P); out[ord] <- ca; names(out) <- names(pn); out
}
