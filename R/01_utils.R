###############################################################################
# MS-BPWAS  -  01_utils.R
# Small self-contained helpers referenced by every other file.
###############################################################################

#' Numerically stable log-sum-exp
#' @keywords internal
logsumexp <- function(x) {
 x <- x[is.finite(x)]
 if (length(x) == 0L) return(-Inf)
 m <- max(x)
 m + log(sum(exp(x - m)))
}

#' Log multivariate-normal density  N(x; mu, Sigma)
#' Optionally accepts pre-computed inverse and log-determinant.
#' @keywords internal
dmvn_log <- function(x, mu, Sigma, Sinv = NULL, logdet = NULL) {
 d <- length(x)
 if (is.null(Sinv))   Sinv   <- solve(Sigma)
 if (is.null(logdet)) logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
 delta <- x - mu
 -0.5 * d * log(2 * pi) - 0.5 * logdet - 0.5 * sum(delta * (Sinv %*% delta))
}

#' Ledoit-Wolf shrinkage of a square matrix toward its diagonal
#' @keywords internal
lw_shrink <- function(S) {
 n <- nrow(S)
 target <- diag(diag(S))
 off <- S - target
 lam <- min(0.5, max(0.01, sum(off^2) / (sum(S^2) * max(n - 1, 1))))
 (1 - lam) * S + lam * target
}

#' Floor eigenvalues so a symmetric matrix becomes PSD
#' @keywords internal
make_psd <- function(M, eps = 1e-6) {
 e <- eigen(M, symmetric = TRUE)
 v <- pmax(e$values, eps)
 out <- e$vectors %*% diag(v, nrow = length(v)) %*% t(e$vectors)
 0.5 * (out + t(out))
}

#' Safe matrix inversion with ridge fall-back
#' @keywords internal
solve_safe <- function(M, ridge = 1e-10) {
 tryCatch(solve(M), error = function(e) solve(M + diag(ridge, nrow(M))))
}

#' Safe Cholesky (adds ridge on failure)
#' @keywords internal
chol_safe <- function(V) {
 tryCatch(chol(V), error = function(e) chol(V + diag(1e-8, nrow(V))))
}

#' Clamp scalar to [lo, hi]
#' @keywords internal
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' Console message (only when verbose)
#' @keywords internal
msg <- function(fmt, ...) cat(sprintf(paste0(fmt, "\n"), ...))

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Precompute R^{-1} and log|R| (done once, shared by all proteins)
#' @keywords internal
precompute_R <- function(R) {
 list(R = R,
      Rinv = solve(R),
      logdetR = as.numeric(determinant(R, logarithm = TRUE)$modulus))
}
