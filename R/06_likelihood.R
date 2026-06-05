###############################################################################
# MS-BPWAS  -  06_likelihood.R
# Marginal likelihood  L(c | z_p)  for one protein / one configuration.
#
# CONTINUOUS NULL (no point mass at zero):
#
#   Under the null configuration c = (0,...,0):
#     mu | c=null  ~  N(0, sigma^2_null * I_S)
#     z | c=null  ~  N(0, sigma^2_null I  +  R)  =  N(0, R_eff)
#
#   Posterior:
#     E[mu | z, c=null]  =  sigma^2_null * R_eff^{-1} * z
#   which is small but never exactly zero.
#
#   Under non-null config c with active set S_c:
#     mu  =  mu_bg  +  A_c * mu_signal
#     mu_bg     ~  N(0, sigma^2_null I_S)       [background, all subtypes]
#     mu_signal ~  Sigma_l pi_l N(0, sigma^2_l I)    [signal, active only]
#
#   Marginal:
#     z | c, l  ~  N(0,  sigma^2_l A_c A_c' + sigma^2_null I + R)
#                = N(0,  sigma^2_l A_c A_c' + R_eff)
#
#   All Woodbury algebra works identically with R_eff replacing R.
#
#   Post_mean is ALWAYS length S.  Post_cov is ALWAYS S x S.
#   No protein ever has exactly zero estimated effect.
###############################################################################


# ===================================================================
# EFFECTIVE-R PRECOMPUTATION  (absorbs sigma^2_null into R)
# ===================================================================

#' Precompute R_eff = sigma^2_null I + R, its inverse, and log|R_eff|.
#' @keywords internal
precompute_Reff <- function(R, sigma2_null) {
 S <- nrow(R)
 Reff <- R + sigma2_null * diag(S)
 list(R       = Reff,
      Rinv    = solve(Reff),
      logdetR = as.numeric(determinant(Reff, logarithm = TRUE)$modulus),
      sigma2_null = sigma2_null,
      R_bare  = R)
}


#' Precompute configuration-specific selector algebra for a fixed R.
#' @keywords internal
prepare_config_cache <- function(cfgs, Rpre) {

 S <- ncol(cfgs)
 Ri <- Rpre$Rinv
 out <- vector("list", nrow(cfgs))

 for (ci in seq_len(nrow(cfgs))) {
   active <- which(cfgs[ci, ] == 1L)
   n_active <- length(active)

   if (n_active == 0L) {
     out[[ci]] <- list(active = integer(0), n_active = 0L)
     next
   }

   out[[ci]] <- list(
     active   = active,
     n_active = n_active,
     AtRiA    = Ri[active, active, drop = FALSE],
     RiAc     = Ri[, active, drop = FALSE]
   )
 }

 out
}


#' One exact E-step shared by VB and EM.
#' @keywords internal
run_exact_estep <- function(z_matrix, cfgs, Rpre, hp, ctrl, lpr,
                            cfg_cache = NULL,
                            objective = c("elbo", "loglik")) {

 objective <- match.arg(objective)
 P <- nrow(z_matrix)
 S <- ncol(z_matrix)
 nc <- nrow(cfgs)
 L <- length(ctrl$sigma2_grid)

 q <- matrix(0, P, nc, dimnames = list(rownames(z_matrix), rownames(cfgs)))
 rsp <- if (ctrl$prior_type %in% c("scale_mixture", "soft_spike_slab"))
          array(0, c(P, nc, L)) else NULL
 pm <- matrix(0, P, S, dimnames = dimnames(z_matrix))
 pv <- matrix(0, P, S, dimnames = dimnames(z_matrix))
 objective_value <- 0

 for (p in seq_len(P)) {
   zp <- z_matrix[p, ]
   if (anyNA(zp)) next

   lqp <- numeric(nc)
   ml_cache <- vector("list", nc)

   for (ci in seq_len(nc)) {
     ml <- compute_marglik(zp, Rpre, cfgs[ci, ], hp, ctrl,
                           cfg_info = if (is.null(cfg_cache)) NULL else cfg_cache[[ci]])
     ml_cache[[ci]] <- ml
     lqp[ci] <- ml$log_ml + lpr[ci]
     if (!is.null(rsp) && !is.null(ml$resp))
       rsp[p, ci, ] <- ml$resp
   }

   mx <- max(lqp[is.finite(lqp)])
   qp <- exp(lqp - mx)
   qp <- qp / sum(qp)
   q[p, ] <- qp

   if (objective == "elbo") {
     pos <- qp > 1e-300
     objective_value <- objective_value + sum(qp[pos] * (lqp[pos] - log(qp[pos])))
   } else {
     objective_value <- objective_value + logsumexp(lqp)
   }

   for (ci in seq_len(nc)) {
     if (qp[ci] < 1e-12) next
     ml <- ml_cache[[ci]]
     pm[p, ] <- pm[p, ] + qp[ci] * ml$post_mean
     pv[p, ] <- pv[p, ] + qp[ci] * (diag(ml$post_cov) + ml$post_mean^2)
   }
   pv[p, ] <- pv[p, ] - pm[p, ]^2
 }

 list(q = q, rsp = rsp, post_mean = pm, post_var_diag = pv,
      objective = objective_value)
}


# ===================================================================
# 1.  SCALE MIXTURE OF NORMALS  (primary / default)
# ===================================================================

#' @return list(log_ml, log_ml_l, resp, post_mean [length S], post_cov [SxS])
#' @keywords internal
marglik_scale_mixture <- function(z, Rpre, config, s2grid, pi_mix,
                                 cfg_info = NULL) {

 S  <- length(z); L <- length(s2grid)
 Ri <- Rpre$Rinv; ldR <- Rpre$logdetR
 s2n <- Rpre$sigma2_null

 Riz <- as.numeric(Ri %*% z)   # used by both null and non-null

 if (is.null(cfg_info)) {
   ac <- which(config == 1)
   nc <- length(ac)
   if (nc > 0L) {
     AtRiA <- Ri[ac, ac, drop = FALSE]
     RiAc  <- Ri[, ac, drop = FALSE]
   }
 } else {
   ac <- cfg_info$active
   nc <- cfg_info$n_active
   if (nc > 0L) {
     AtRiA <- cfg_info$AtRiA
     RiAc  <- cfg_info$RiAc
   }
 }

 ## ==============================================================
 ## NULL configuration: z ~ N(0, R_eff)
 ## ==============================================================
 if (nc == 0L) {
   lml    <- dmvn_log(z, rep(0, S), Rpre$R, Ri, ldR)
   m_null <- s2n * Riz                          # small but non-zero
   V_null <- s2n * (diag(S) - s2n * Ri)

   return(list(log_ml = lml, log_ml_l = rep(lml, L), resp = rep(1/L, L),
               post_mean = m_null, post_cov = V_null))
 }

 ## ==============================================================
 ## NON-NULL configuration:  Sigma_{c,l} = sigma^2_l A A' + R_eff
 ## ==============================================================
 AtRiz <- Riz[ac]
 ztRiz <- as.numeric(crossprod(z, Ri %*% z))

 log_ml_l <- numeric(L)
 pm_sig_l <- matrix(0, L, nc)
 pV_sig_l <- array(0, c(L, nc, nc))

 for (l in seq_len(L)) {
   s2 <- s2grid[l]
   M  <- diag(1/s2, nc) + AtRiA
   Mi <- solve_safe(M)

   ldf  <- as.numeric(determinant(diag(nc) + s2 * AtRiA, logarithm = TRUE)$modulus)
   quad <- ztRiz - as.numeric(crossprod(AtRiz, Mi %*% AtRiz))

   log_ml_l[l]    <- -0.5*S*log(2*pi) - 0.5*(ldR + ldf) - 0.5*quad
   pV_sig_l[l,,]  <- Mi
   pm_sig_l[l, ]  <- Mi %*% AtRiz
 }

 ## mixture integration
 lt   <- log(pi_mix) + log_ml_l
 lml  <- logsumexp(lt)
 resp <- exp(lt - lml)

 ## ---------- Full S-dim posterior (background + signal) ----------
 ##
 ## E[mu | z, c, l]  =  background(all S)  +  Ac * signal(|c|)
 ##
 ## background  =  sigma^2_null * Sigma_{c,l}^{-1} z
 ##             =  sigma^2_null * (Ri z - Ri Ac Mi AtRiz)     [Woodbury]
 ##
 ## signal projected to S  =  Ac * pm_sig_l[l, ]

 pm_full_l <- matrix(0, L, S)
 for (l in seq_len(L)) {
   Mi_l <- matrix(pV_sig_l[l,,], nc, nc)
   bg   <- s2n * (Riz - as.numeric(RiAc %*% Mi_l %*% AtRiz))
   sig  <- numeric(S)
   sig[ac] <- pm_sig_l[l, ]
   pm_full_l[l, ] <- bg + sig
 }

 ## Mixture-averaged posterior mean
 pm <- as.numeric(crossprod(resp, pm_full_l))

 ## Mixture-averaged posterior covariance (full S x S)
 pV <- matrix(0, S, S)
 for (l in seq_len(L))
   pV <- pV + resp[l] * tcrossprod(pm_full_l[l, ])
 pV <- pV - tcrossprod(pm)

 ## Add within-component variance contribution
 for (l in seq_len(L)) {
   Mi_l <- matrix(pV_sig_l[l,,], nc, nc)
   ## Signal covariance projected to S dimensions
   sig_cov <- matrix(0, S, S)
   sig_cov[ac, ac] <- Mi_l
   ## Background covariance: sigma^2_null(I - sigma^2_null R_eff^{-1})  (diagonal approx)
   bg_cov <- s2n * (diag(S) - s2n * Ri)
   pV <- pV + resp[l] * (sig_cov + bg_cov)
 }

 list(log_ml = lml, log_ml_l = log_ml_l, resp = resp,
      post_mean = pm, post_cov = pV)
}


# ===================================================================
# 2.  SOFT SPIKE-AND-SLAB  (two-component special case)
# ===================================================================

#' @keywords internal
marglik_soft_ss <- function(z, Rpre, config,
                             s2_small, s2_large, pi_small = 0.5,
                             cfg_info = NULL) {
 marglik_scale_mixture(z, Rpre, config,
                        s2grid = c(s2_small, s2_large),
                        pi_mix = c(pi_small, 1 - pi_small),
                        cfg_info = cfg_info)
}


# ===================================================================
# 3.  HORSESHOE  (numerical quadrature)
# ===================================================================

#' @keywords internal
marglik_horseshoe <- function(z, Rpre, config, lambda = 1, nq = 20L,
                              cfg_info = NULL) {

 S  <- length(z)
 Ri <- Rpre$Rinv; ldR <- Rpre$logdetR
 s2n <- Rpre$sigma2_null
 Riz <- as.numeric(Ri %*% z)

 if (is.null(cfg_info)) {
   ac <- which(config == 1)
   nc <- length(ac)
   if (nc > 0L) {
     AtRiA <- Ri[ac, ac, drop = FALSE]
     RiAc  <- Ri[, ac, drop = FALSE]
   }
 } else {
   ac <- cfg_info$active
   nc <- cfg_info$n_active
   if (nc > 0L) {
     AtRiA <- cfg_info$AtRiA
     RiAc  <- cfg_info$RiAc
   }
 }

 ## Null: same soft-null as scale mixture
 if (nc == 0L) {
   lml    <- dmvn_log(z, rep(0, S), Rpre$R, Ri, ldR)
   m_null <- s2n * Riz
   V_null <- s2n * (diag(S) - s2n * Ri)
   return(list(log_ml = lml, post_mean = m_null, post_cov = V_null))
 }

 AtRiz <- Riz[ac]
 ztRiz <- as.numeric(crossprod(z, Ri %*% z))

 log_t2 <- seq(log(1e-4), log(100), length.out = nq)
 t2g    <- exp(log_t2)
 dlog   <- diff(log_t2)[1]
 log_hc <- function(t2) -log(pi) - 0.5*log(t2) - log(1 + t2)

 if (nc <= 3L) {
   gl <- replicate(nc, seq_len(nq), simplify = FALSE)
   ig <- as.matrix(expand.grid(gl)); np <- nrow(ig)
   lv <- numeric(np)
   for (q in seq_len(np)) {
     t2v <- t2g[ig[q,]] * lambda^2
     Dt  <- diag(t2v, nc)
     M   <- solve_safe(Dt) + AtRiA; Mi <- solve_safe(M)
     ldf <- as.numeric(determinant(diag(nc) + Dt %*% AtRiA, logarithm = TRUE)$modulus)
     qf  <- ztRiz - as.numeric(crossprod(AtRiz, Mi %*% AtRiz))
     lv[q] <- -0.5*S*log(2*pi) - 0.5*(ldR+ldf) - 0.5*qf +
              sum(vapply(t2g[ig[q,]], log_hc, 0)) + nc*log(dlog)
   }
   lml  <- logsumexp(lv) - nc*log(nq)
   best <- which.max(lv)
   t2_map <- t2g[ig[best,]] * lambda^2
 } else {
   t2_cur <- rep(0.1*lambda^2, nc)
   for (mf in 1:10) for (d in seq_len(nc)) {
     lv <- numeric(nq)
     for (q in seq_len(nq)) {
       t2_try <- t2_cur; t2_try[d] <- t2g[q]*lambda^2
       Sig_c  <- Rpre$R
       Sig_c[ac, ac] <- Sig_c[ac, ac, drop = FALSE] + diag(t2_try, nc)
       lv[q]  <- dmvn_log(z, rep(0,S), Sig_c) + log_hc(t2g[q])
     }
     t2_cur[d] <- t2g[which.max(lv)] * lambda^2
   }
   t2_map <- t2_cur
   Sig_m  <- Rpre$R
   Sig_m[ac, ac] <- Sig_m[ac, ac, drop = FALSE] + diag(t2_map, nc)
   lml    <- dmvn_log(z, rep(0,S), Sig_m) +
             sum(vapply(t2_map/lambda^2, log_hc, 0))
 }

 Mf       <- solve_safe(diag(t2_map, nc)) + AtRiA
 Mi       <- solve_safe(Mf)
 sig_mean <- as.numeric(Mi %*% AtRiz)

 ## Full S-dim posterior
 bg   <- s2n * (Riz - as.numeric(RiAc %*% Mi %*% AtRiz))
 sig  <- numeric(S)
 sig[ac] <- sig_mean
 pm   <- bg + sig

 sig_cov <- matrix(0, S, S)
 sig_cov[ac, ac] <- Mi
 bg_cov  <- s2n * (diag(S) - s2n * Ri)
 pV      <- sig_cov + bg_cov

 list(log_ml = lml, post_mean = pm, post_cov = pV)
}


# ===================================================================
# 4.  UNIFIED DISPATCHER
# ===================================================================

#' @keywords internal
compute_marglik <- function(z, Rpre, config, hyper, ctrl, cfg_info = NULL) {
 switch(ctrl$prior_type,
   scale_mixture   = marglik_scale_mixture(z, Rpre, config,
                       ctrl$sigma2_grid, hyper$pi_mix, cfg_info = cfg_info),
   soft_spike_slab = marglik_soft_ss(z, Rpre, config,
                       ctrl$sigma2_small, ctrl$sigma2_large,
                       hyper$pi_small %||% 0.5, cfg_info = cfg_info),
   horseshoe       = marglik_horseshoe(z, Rpre, config,
                       hyper$global_scale %||% 1, cfg_info = cfg_info),
   stop("Unknown prior: ", ctrl$prior_type)
 )
}
