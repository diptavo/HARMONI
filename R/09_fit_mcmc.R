###############################################################################
# MS-BPWAS  -  09_fit_mcmc.R
# Full MCMC via Gibbs sampling.
#
# Blocks per iteration:
#   1  c_p      : collapsed (mu integrated out) multinomial draw
#   2  mu_p|c_p  : conjugate normal draw
#   3  l_p|mu_p  : component assignment (scale mixture only)
#   4  globals  : pi0, xi_s, pi_mix from conjugate full-conditionals
###############################################################################

#' @keywords internal
fit_mcmc <- function(z_matrix, R, cfgs, ctrl) {

 P<-nrow(z_matrix); S<-ncol(z_matrix); nc<-nrow(cfgs)
 L<-length(ctrl$sigma2_grid)
 n_ch<-ctrl$n_chains; ni<-ctrl$n_iter; nw<-ctrl$n_warmup; th<-ctrl$thin
 nk <- (ni-nw) %/% th
 Rpre <- precompute_Reff(R, ctrl$sigma2_null %||% 0.005^2)

 if(ctrl$verbose)
   msg("  MCMC: %d chains x %d iter (%d warmup, thin %d -> %d kept)",
       n_ch, ni, nw, th, nk)

 chains <- vector("list", n_ch)
 for (ch in seq_len(n_ch)) {
   if(ctrl$verbose) msg("    chain %d ...", ch)
   chains[[ch]] <- run_chain(z_matrix,Rpre,cfgs,ctrl,P,S,nc,L,ni,nw,th,nk,
                              seed_off=ch*1000)
 }

 ## combine
 q_c  <- Reduce("+", lapply(chains, `[[`, "q")) / n_ch
 pm_c <- Reduce("+", lapply(chains, `[[`, "pm")) / n_ch
 pv_c <- Reduce("+", lapply(chains, `[[`, "pv")) / n_ch
 hp <- avg_hyper(chains)

 rh <- rhat_check(chains)
 if(ctrl$verbose) msg("    R-hat(pi0)=%.3f  R-hat(xi)=[%s]",
   rh$rh_pi0, paste(sprintf("%.3f",rh$rh_xi),collapse=","))

 list(q=q_c, hyperparams=hp, post_mean=pm_c, post_var_diag=pv_c,
      convergence=list(converged=all(rh$rh_pi0<1.05), n_iter=ni,
                       rhat=rh, n_chains=n_ch, n_kept=nk),
      chains=chains)
}

# --- single chain ---
#' @keywords internal
run_chain <- function(z_matrix,Rpre,cfgs,ctrl,P,S,nc,L,ni,nw,th,nk,seed_off) {
 set.seed(42+seed_off)
 sn <- colnames(z_matrix)

 c_cur <- rep(1L, P)                    # all null
 mu_cur<- matrix(0, P, S)
 l_cur <- sample(L, P, replace=TRUE)
 pi0   <- ctrl$pi0_init
 xi    <- rep(ctrl$xi_init, S)
 pmix  <- rep(1/L, L)
 lam   <- 1

 q_acc  <- matrix(0, P, nc)
 pm_acc <- matrix(0, P, S)
 pv_acc <- matrix(0, P, S)
 pi0_s  <- numeric(nk); xi_s <- matrix(0, nk, S)
 pmix_s <- matrix(0, nk, L)
 kept <- 0L

 for (it in seq_len(ni)) {
   hp  <- list(pi0=pi0, xi=xi, pi_mix=pmix,
               pi_small=0.5, global_scale=lam)
   lpr <- log(config_prior(cfgs, pi0, xi))

   ## Block 1: sample c_p (collapsed) ---------
   for (p in seq_len(P)) {
     zp <- z_matrix[p,]; if(anyNA(zp)) next
     lp <- numeric(nc)
     for(ci in seq_len(nc))
       lp[ci] <- compute_marglik(zp, Rpre, cfgs[ci,], hp, ctrl)$log_ml + lpr[ci]
     mx<-max(lp[is.finite(lp)]); pr<-exp(lp-mx); pr<-pr/sum(pr)
     c_cur[p] <- sample.int(nc, 1, prob=pr)
   }

   ## Block 2: sample mu_p | c_p ---------------
   s2n <- ctrl$sigma2_null %||% 0.005^2
   for (p in seq_len(P)) {
     zp<-z_matrix[p,]; ci<-c_cur[p]
     act<-which(cfgs[ci,]==1)

     if(!length(act)) {
       ## Soft null: mu | z, c=null ~ N(sigma^2_null R_eff^{-1} z, V_null)
       m_null <- s2n * as.numeric(Rpre$Rinv %*% zp)
       V_null <- s2n * (diag(S) - s2n * Rpre$Rinv)
       mu_cur[p,] <- as.numeric(m_null + chol_safe(V_null) %*% rnorm(S))
       next
     }

     s2 <- if(ctrl$prior_type=="scale_mixture") ctrl$sigma2_grid[l_cur[p]]
           else if(ctrl$prior_type=="soft_spike_slab")
             ifelse(runif(1)<hp$pi_small, ctrl$sigma2_small, ctrl$sigma2_large)
           else lam^2*0.1

     Ac<-matrix(0,S,length(act)); for(j in seq_along(act)) Ac[act[j],j]<-1
     M  <- diag(1/s2, length(act)) + crossprod(Ac, Rpre$Rinv %*% Ac)
     V  <- solve_safe(M)
     m  <- V %*% crossprod(Ac, Rpre$Rinv %*% zp)
     sig_draw <- as.numeric(m + chol_safe(V) %*% rnorm(length(act)))

     ## Full S-dim: background + signal
     Mi_AtRiz <- V %*% crossprod(Ac, Rpre$Rinv %*% zp)
     Riz <- as.numeric(Rpre$Rinv %*% zp)
     RiAc <- Rpre$Rinv %*% Ac
     bg_mean <- s2n * (Riz - as.numeric(RiAc %*% Mi_AtRiz))
     bg_cov  <- s2n * (diag(S) - s2n * Rpre$Rinv)  # approximate
     bg_draw <- as.numeric(bg_mean + chol_safe(bg_cov) %*% rnorm(S))

     mu_cur[p,] <- bg_draw + as.numeric(Ac %*% sig_draw)
   }

   ## Block 3: sample l_p | mu_p (scale mix) ---
   if (ctrl$prior_type == "scale_mixture") {
     for (p in seq_len(P)) {
       act <- which(cfgs[c_cur[p],]==1)
       if(!length(act)) { l_cur[p]<-sample(L,1); next }
       ma <- mu_cur[p, act]
       lpl <- log(pmix) - 0.5*length(act)*log(ctrl$sigma2_grid) -
              0.5*sum(ma^2)/ctrl$sigma2_grid
       mx<-max(lpl); pr<-exp(lpl-mx); l_cur[p]<-sample.int(L,1,prob=pr/sum(pr))
     }
   }

   ## Block 4: globals -------------------------
   nn <- sum(c_cur==1)
   a0 <- ctrl$pi0_prior_a %||% 10; b0 <- ctrl$pi0_prior_b %||% 2
   pi0 <- clamp(rbeta(1, a0+nn, b0+P-nn), 0.01, 0.999)

   nn_mask <- c_cur != 1; n_nn <- sum(nn_mask)
   for (s in seq_len(S)) {
     ms <- sum(cfgs[c_cur[nn_mask], s])
     xi[s] <- clamp(rbeta(1, 1+ms, 1+max(n_nn-ms,0)), 0.005, 0.995)
   }
   if (ctrl$prior_type=="scale_mixture") {
     ct <- tabulate(l_cur[nn_mask], nbins=L)
     pmix <- rgamma(L, 1+ct, 1); pmix <- pmax(pmix,1e-8)
     pmix <- pmix/sum(pmix)
   }

   ## accumulate post-warmup --------------------
   if (it > nw && (it-nw)%%th == 0) {
     kept <- kept + 1L
     for(p in seq_len(P)) q_acc[p, c_cur[p]] <- q_acc[p, c_cur[p]] + 1
     pm_acc <- pm_acc + mu_cur
     pv_acc <- pv_acc + mu_cur^2
     pi0_s[kept] <- pi0; xi_s[kept,] <- xi
     if(ctrl$prior_type=="scale_mixture") pmix_s[kept,] <- pmix
   }
   if (ctrl$verbose && it%%500==0)
     msg("      iter %d: %d null, pi0=%.3f", it, nn, pi0)
 }

 list(q  = q_acc / nk,
      pm = pm_acc / nk,
      pv = pv_acc / nk - (pm_acc/nk)^2,
      pi0_s = pi0_s, xi_s = xi_s, pmix_s = pmix_s)
}

# --- combine / diagnose helpers ---
#' @keywords internal
avg_hyper <- function(chs) {
 list(pi0    = mean(sapply(chs, function(ch) mean(ch$pi0_s))),
      xi     = colMeans(do.call(rbind, lapply(chs,
                 function(ch) colMeans(ch$xi_s)))),
      pi_mix = colMeans(do.call(rbind, lapply(chs,
                 function(ch) colMeans(ch$pmix_s)))))
}

#' @keywords internal
rhat_check <- function(chs) {
 rh1 <- function(x) {
   if(!is.matrix(x)) x<-as.matrix(x); n<-nrow(x); m<-ncol(x)
   if(n<2||m<2) return(NA_real_)
   B<-n*var(colMeans(x)); W<-mean(apply(x,2,var))
   sqrt(((1-1/n)*W + B/n) / W)
 }
 pi0_mat <- sapply(chs, `[[`, "pi0_s")
 S <- ncol(chs[[1]]$xi_s)
 rh_xi <- numeric(S)
 for(s in seq_len(S)) rh_xi[s] <- rh1(sapply(chs, function(ch) ch$xi_s[,s]))
 list(rh_pi0 = rh1(pi0_mat), rh_xi = rh_xi)
}
