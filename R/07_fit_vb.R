###############################################################################
# MS-BPWAS  -  07_fit_vb.R
# Variational Bayes (default).  E-step is exact (no approximation needed
# because proteins are conditionally independent given hyperparameters).
###############################################################################

#' @keywords internal
fit_vb <- function(z_matrix, R, cfgs, ctrl) {

 P <- nrow(z_matrix); S <- ncol(z_matrix); nc <- nrow(cfgs)
 L <- length(ctrl$sigma2_grid); sn <- colnames(z_matrix)
 csize <- rowSums(cfgs)
 ## Use effective R (absorbs sigma^2_null background into the correlation)
 Rpre <- precompute_Reff(R, ctrl$sigma2_null %||% 0.005^2)
 cfg_cache <- prepare_config_cache(cfgs, Rpre)

 ## initialise hyper-parameters
 hp <- list(pi0 = ctrl$pi0_init,
            xi  = setNames(rep(ctrl$xi_init, S), sn),
            pi_mix = rep(1/L, L),
            pi_small = 0.5, global_scale = 1)

 q   <- matrix(0, P, nc, dimnames=list(rownames(z_matrix), rownames(cfgs)))
 rsp <- if (ctrl$prior_type %in% c("scale_mixture","soft_spike_slab"))
          array(0, c(P, nc, L)) else NULL
 pm  <- matrix(0, P, S, dimnames=dimnames(z_matrix))
 pv  <- matrix(0, P, S, dimnames=dimnames(z_matrix))
 elbo_h <- numeric(ctrl$max_iter)

 if (ctrl$verbose) msg("  VB: P=%d S=%d cfgs=%d L=%d  sigma2_null=%.2e",
                        P, S, nc, L, ctrl$sigma2_null %||% 0.005^2)

 for (it in seq_len(ctrl$max_iter)) {

   lpr  <- log(config_prior(cfgs, hp$pi0, hp$xi))
   estep <- run_exact_estep(z_matrix, cfgs, Rpre, hp, ctrl, lpr,
                            cfg_cache = cfg_cache, objective = "elbo")
   q   <- estep$q
   rsp <- estep$rsp
   pm  <- estep$post_mean
   pv  <- estep$post_var_diag
   elbo <- estep$objective
   elbo_h[it] <- elbo

   ## ---- M-step -------------------------------------------------
   ## pi0 with Beta(a0, b0) conjugate prior to prevent collapse
   n_null_post <- sum(q[,1])
   a0 <- ctrl$pi0_prior_a %||% 10; b0 <- ctrl$pi0_prior_b %||% 2
   hp$pi0 <- clamp((a0 - 1 + n_null_post) / (a0 + b0 - 2 + P), 0.01, 0.999)

   tn <- sum(1 - q[,1])
   if (tn > 1) for (s in seq_len(S)) {
     cs <- which(cfgs[,s]==1)
     hp$xi[s] <- clamp(sum(q[,cs]) / tn, 0.005, 0.995)
   }

   if (ctrl$prior_type == "scale_mixture" && !is.null(rsp)) {
     np <- numeric(L); tw <- 0
     for (p in seq_len(P)) for (ci in 2:nc) {
       if (q[p,ci] < 1e-12) next
       np <- np + q[p,ci]*rsp[p,ci,]; tw <- tw + q[p,ci]
     }
     if (tw > 0) { hp$pi_mix <- pmax(np/tw, 1e-8)
                    hp$pi_mix <- hp$pi_mix / sum(hp$pi_mix) }
   }
   if (ctrl$prior_type == "soft_spike_slab" && !is.null(rsp)) {
     ns<-0; tw<-0
     for(p in seq_len(P)) for(ci in 2:nc) {
       if(q[p,ci]<1e-12) next
       ns<-ns+q[p,ci]*rsp[p,ci,1]; tw<-tw+q[p,ci]
     }
     if(tw>0) hp$pi_small <- clamp(ns/tw, 0.01, 0.99)
   }

   ## ---- convergence check --------------------------------------
   if (it >= 3L) {
     rel <- abs(elbo_h[it]-elbo_h[it-1]) / max(abs(elbo_h[it-1]),1)
     if (ctrl$verbose && (it%%10==0 || it<=3))
       msg("    VB %3d | ELBO %.4f | Delta %.1e | pi0 %.3f | xi [%s]",
           it, elbo, rel, hp$pi0,
           paste(sprintf("%.2f",hp$xi), collapse=","))
     if (rel < ctrl$tol) { msg("    VB converged iter %d", it); break }
   }
 }

 list(q=q, hyperparams=hp, post_mean=pm, post_var_diag=pv,
      convergence=list(converged = it<ctrl$max_iter, n_iter=it,
                       elbo_final=elbo_h[it], elbo_hist=elbo_h[1:it]))
}
