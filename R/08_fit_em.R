###############################################################################
# MS-BPWAS  -  08_fit_em.R
# EM algorithm.  Structurally identical to VB (the E-step is already
# exact), but monitors the observed-data log-likelihood Q rather than
# the ELBO and guarantees Q is non-decreasing.
###############################################################################

#' @keywords internal
fit_em <- function(z_matrix, R, cfgs, ctrl) {

 P <- nrow(z_matrix); S <- ncol(z_matrix); nc <- nrow(cfgs)
 L <- length(ctrl$sigma2_grid); sn <- colnames(z_matrix)
 Rpre <- precompute_Reff(R, ctrl$sigma2_null %||% 0.005^2)
 cfg_cache <- prepare_config_cache(cfgs, Rpre)

 hp <- list(pi0=ctrl$pi0_init, xi=setNames(rep(ctrl$xi_init,S),sn),
            pi_mix=rep(1/L,L), pi_small=0.5, global_scale=1)

 q   <- matrix(0, P, nc)
 rsp <- if (ctrl$prior_type %in% c("scale_mixture","soft_spike_slab"))
          array(0, c(P,nc,L)) else NULL
 pm  <- matrix(0, P, S, dimnames=dimnames(z_matrix))
 pv  <- matrix(0, P, S, dimnames=dimnames(z_matrix))
 Qh  <- numeric(ctrl$max_iter)

 if (ctrl$verbose) msg("  EM: P=%d S=%d cfgs=%d", P, S, nc)

 for (it in seq_len(ctrl$max_iter)) {
   lpr <- log(config_prior(cfgs, hp$pi0, hp$xi))
   estep <- run_exact_estep(z_matrix, cfgs, Rpre, hp, ctrl, lpr,
                            cfg_cache = cfg_cache, objective = "loglik")
   q   <- estep$q
   rsp <- estep$rsp
   pm  <- estep$post_mean
   pv  <- estep$post_var_diag
   Qobs <- estep$objective
   Qh[it] <- Qobs

   ## M-step (same as VB, with Beta prior on pi0)
   n_null_post <- sum(q[,1])
   a0 <- ctrl$pi0_prior_a %||% 10; b0 <- ctrl$pi0_prior_b %||% 2
   hp$pi0 <- clamp((a0 - 1 + n_null_post) / (a0 + b0 - 2 + P), 0.01, 0.999)
   tn <- sum(1-q[,1])
   if (tn>1) for(s in seq_len(S)) {
     cs<-which(cfgs[,s]==1); hp$xi[s]<-clamp(sum(q[,cs])/tn, 0.005, 0.995)
   }
   if (ctrl$prior_type=="scale_mixture"&&!is.null(rsp)) {
     np<-numeric(L); tw<-0
     for(p in seq_len(P)) for(ci in 2:nc) {
       if(q[p,ci]<1e-12) next; np<-np+q[p,ci]*rsp[p,ci,]; tw<-tw+q[p,ci]
     }
     if(tw>0){hp$pi_mix<-pmax(np/tw,1e-8); hp$pi_mix<-hp$pi_mix/sum(hp$pi_mix)}
   }
   if (ctrl$prior_type=="soft_spike_slab"&&!is.null(rsp)) {
     ns<-0; tw<-0
     for(p in seq_len(P)) for(ci in 2:nc) {
       if(q[p,ci]<1e-12) next; ns<-ns+q[p,ci]*rsp[p,ci,1]; tw<-tw+q[p,ci]
     }
     if(tw>0) hp$pi_small<-clamp(ns/tw,0.01,0.99)
   }

   ## convergence
   if (it>=3L) {
     rel <- abs(Qh[it]-Qh[it-1])/max(abs(Qh[it-1]),1)
     if (Qh[it] < Qh[it-1]-1e-4)
       warning("EM Q decreased at iter ", it)
     if (ctrl$verbose && (it%%10==0||it<=3))
       msg("    EM %3d | Q %.4f | Delta %.1e | pi0 %.3f", it,Qobs,rel,hp$pi0)
     if (rel < ctrl$tol) { msg("    EM converged iter %d",it); break }
   }
 }

 list(q=q, hyperparams=hp, post_mean=pm, post_var_diag=pv,
      convergence=list(converged=it<ctrl$max_iter, n_iter=it,
                       Q_final=Qh[it], Q_hist=Qh[1:it]))
}
