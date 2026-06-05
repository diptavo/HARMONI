#' Construct shared + heterogeneity z-matrix for S > 2 subtypes.
#'
#' For S subtypes + 1 overall GWAS, construct:
#'   - 1 shared component (the overall z-statistics)
#'   - S-1 orthogonal heterogeneity contrasts
#'
#' @param z_overall   P-vector of z-stats from overall GWAS
#' @param z_subtypes  P x S matrix of z-stats from subtype GWAS
#'                    (columns = subtypes)
#' @param R_subtypes  S x S null correlation matrix among subtype z-stats
#'                    (estimated empirically from the z_subtypes matrix)
#' @param contrasts   (S-1) x S contrast matrix.  If NULL, automatically
#'                    generates Helmert contrasts.
#' @return list with:
#'   z_matrix:   P x S matrix (1 shared + S-1 heterogeneity)
#'   R_null:     S x S null correlation for the transformed z-matrix
#'   C:          the contrast matrix used
#'   back_map:   function to recover subtype interpretation
construct_shared_het_z <- function(z_overall,
                                   z_subtypes,
                                   R_subtypes = NULL,
                                   contrasts  = NULL) {
  
  P <- nrow(z_subtypes)
  S <- ncol(z_subtypes)
  sn <- colnames(z_subtypes)
  
  stopifnot(length(z_overall) == P)
  
  # ---- 1. Estimate null correlation among subtype z-stats ----
  if (is.null(R_subtypes)) {
    # Use empirical correlation (most z's are null)
    R_subtypes <- cor(z_subtypes)
    cat("Estimated null correlation among subtypes:\n")
    print(round(R_subtypes, 3))
  }
  
  # ---- 2. Build contrast matrix ----
  if (is.null(contrasts)) {
    # Helmert contrasts: each compares one subtype to the
    # average of all preceding subtypes
    #
    # For S=3 (LUAD, SQC, SCLC):
    #   C1 = [1, -1,  0]        LUAD vs SQC
    #   C2 = [1,  1, -2] / ...  SCLC vs avg(LUAD, SQC)
    #
    # General Helmert:
    #   C_k compares subtype (k+1) to the mean of subtypes 1..k
    C <- matrix(0, S - 1, S)
    for (k in 1:(S - 1)) {
      C[k, 1:k]   <-  1 / k
      C[k, k + 1] <- -1
    }
    rownames(C) <- paste0("het_", 1:(S-1))
    colnames(C) <- sn
  } else {
    C <- contrasts
    stopifnot(nrow(C) == S - 1, ncol(C) == S)
  }
  
  cat("Contrast matrix:\n")
  print(round(C, 3))
  
  # ---- 3. Compute heterogeneity z-statistics ----
  # For contrast c_k (1 x S):
  #   raw_k = z_subtypes %*% c_k'   (P x 1)
  #   Var(raw_k) = c_k R_subtypes c_k'  (scalar, under the null)
  #   z_het_k = raw_k / sqrt(Var(raw_k))
  
  z_het <- matrix(NA_real_, P, S - 1)
  het_var <- numeric(S - 1)
  
  for (k in 1:(S - 1)) {
    c_k <- C[k, ]
    raw <- z_subtypes %*% c_k           # P x 1
    v   <- as.numeric(c_k %*% R_subtypes %*% c_k)  # null variance
    het_var[k] <- v
    z_het[, k] <- as.numeric(raw) / sqrt(v)
  }
  
  colnames(z_het) <- rownames(C)
  
  # ---- 4. Assemble full z-matrix ----
  z_matrix <- cbind(shared = z_overall, z_het)
  
  # ---- 5. Compute null correlation for the full z-matrix ----
  # R_null[shared, shared] = 1
  # R_null[shared, het_k]  = Cov(z_overall, z_het_k) / sqrt(Var * 1)
  #
  # z_overall approx Sigma_s pi_s z_s + pi_res z_res, so:
  #   Cov(z_overall, c_k' z_sub) depends on the prevalences and
  #   sample overlap structure
  #
  # Estimate empirically: safest approach
  
  R_null <- cor(z_matrix)
  
  # Regularize: force diagonal to 1
  d <- sqrt(diag(R_null))
  R_null <- R_null / outer(d, d)
  
  cat("Null correlation of transformed z-matrix:\n")
  print(round(R_null, 3))
  
  # ---- 6. Build interpretation helper ----
  back_map <- function(config, z_row) {
    # Given a configuration label and the z-values for one protein,
    # return a biological interpretation
    has_shared <- grepl("shared", config)
    het_active <- sapply(rownames(C), function(h) grepl(h, config))
    
    if (config == "null") return("null")
    if (has_shared && !any(het_active)) return("pan-cancer")
    
    # Determine which subtypes drive the heterogeneity
    drivers <- character(0)
    for (k in which(het_active)) {
      c_k <- C[k, ]
      # Positive z_het means the positive side of contrast is larger
      pos_side <- sn[c_k > 0]
      neg_side <- sn[c_k < 0]
      if (z_row[k + 1] > 0) {
        drivers <- c(drivers, paste(pos_side, "high"))
      } else {
        drivers <- c(drivers, paste(neg_side, "high"))
      }
    }
    
    if (has_shared) {
      paste0("shared + ", paste(drivers, collapse = ", "))
    } else {
      paste0("heterogeneous: ", paste(drivers, collapse = ", "))
    }
  }
  
  list(z_matrix = z_matrix,
       R_null   = R_null,
       C        = C,
       back_map = back_map,
       het_var  = het_var)
}