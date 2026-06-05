###############################################################################
# MS-BPWAS  -  04_pwas_z.R
# Compute the P x S matrix of PWAS Z-statistics.
#
#   z_ps  =  w_p' z_s  /  sqrt( w_p' V_p w_p )
#
# w_p = pretrained cis-prediction weights
# z_s = GWAS Z-scores at those SNPs for subtype s
# V_p = LD correlation among the SNPs (from PLINK reference)
###############################################################################

#' Compute PWAS Z-statistics for every protein x subtype
#' @return list(z_matrix, protein_info)
#' @keywords internal
compute_all_pwas_z <- function(gwas, pmodels, plink_prefix,
                                subtype_names, ctrl) {

 P <- length(pmodels)
 S <- length(subtype_names)
 pids <- vapply(pmodels, function(m) m$protein_id, "")

 zmat <- matrix(NA_real_, P, S, dimnames = list(pids, subtype_names))

 info <- data.frame(protein_id = pids,
   gene = vapply(pmodels, function(m) m$gene %||% "", ""),
   chr  = vapply(pmodels, function(m) as.integer(m$chr %||% 0L), 0L),
   n_snps = integer(P), cv_R2 = numeric(P), wVw = numeric(P),
   stringsAsFactors = FALSE)

 for (p in seq_len(P)) {
   ld <- get_ld(pmodels[[p]], plink_prefix, ctrl)
   if (is.null(ld) || length(ld$w) < ctrl$min_n_snps) next

   wVw <- as.numeric(crossprod(ld$w, ld$V %*% ld$w))
   if (wVw <= 0) next
   denom <- sqrt(wVw)

   for (s in seq_len(S)) {
     zs <- align_z(ld$snps, ld$a1, gwas[[ subtype_names[s] ]])
     if (is.null(zs) || length(zs) < ctrl$min_n_snps) next
     zmat[p, s] <- sum(ld$w * zs) / denom
   }
   info$n_snps[p] <- length(ld$w)
   info$cv_R2[p]  <- pmodels[[p]]$cv_R2 %||% NA_real_
   info$wVw[p]    <- wVw

   if (ctrl$verbose && p %% 200 == 0) msg("    %d / %d proteins", p, P)
 }

 ok <- rowSums(is.na(zmat)) == 0L
 list(z_matrix = zmat[ok, , drop = FALSE], protein_info = info[ok, ])
}


# ---- LD for one protein --------------------------------------------------

#' Call PLINK to get LD correlation among a protein's cis-SNPs
#' @return list(V, w, snps, a1) or NULL
#' @keywords internal
get_ld <- function(model, plink_prefix, ctrl) {

 snps <- model$weights$SNP
 chr  <- model$chr
 if (is.null(chr) || length(snps) < 2L) return(NULL)

 tmp_snp <- tempfile(fileext = ".snplist")
 writeLines(snps, tmp_snp)
 tmp_out <- tempfile()

 # try plink 1.9
 cmd <- sprintf(
   "plink --bfile %s --chr %s --extract %s --r square --out %s --silent 2>/dev/null",
   plink_prefix, chr, tmp_snp, tmp_out)
 system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

 ld_file <- paste0(tmp_out, ".ld")
 if (!file.exists(ld_file)) {
   # try plink 2
   cmd2 <- sprintf(
     "plink2 --bfile %s --chr %s --extract %s --r-unphased square --out %s 2>/dev/null",
     plink_prefix, chr, tmp_snp, tmp_out)
   system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)
   ld_file <- paste0(tmp_out, ".unphased.vcor")
 }
 V_raw <- tryCatch(as.matrix(read.table(ld_file)), error = function(e) NULL)
 unlink(c(tmp_snp, list.files(dirname(tmp_out),
          pattern = basename(tmp_out), full.names = TRUE)))
 if (is.null(V_raw) || nrow(V_raw) < 2L) return(NULL)

 bim <- read.table(paste0(plink_prefix, ".bim"),
                     header = FALSE, stringsAsFactors = FALSE)
 colnames(bim) <- c("CHR","SNP","CM","POS","A1","A2")
 bim <- bim[bim$SNP %in% snps & bim$CHR == chr, ]

 m <- min(nrow(V_raw), nrow(bim))
 bim <- bim[seq_len(m), ]; V_raw <- V_raw[seq_len(m), seq_len(m)]

 common <- intersect(bim$SNP, snps)
 if (length(common) < 2L) return(NULL)
 ib <- match(common, bim$SNP)
 im <- match(common, model$weights$SNP)

 V <- V_raw[ib, ib, drop = FALSE]
 w <- model$weights$weight[im]

 # allele alignment: flip weight when model A1 != ref A1
 flip <- toupper(model$weights$A1[im]) != toupper(bim$A1[ib])
 w[flip] <- -w[flip]

 if (ctrl$ld_shrinkage && nrow(V) > 3L) V <- lw_shrink(V)
 V <- make_psd(V)

 list(V = V, w = w, snps = common, a1 = toupper(bim$A1[ib]))
}


# ---- Allele-aligned Z extraction ----------------------------------------

#' Extract Z-scores from a GWAS, aligned to reference alleles
#' @return numeric vector (same length as snps) or NULL
#' @keywords internal
align_z <- function(snps, a1_ref, gwas) {

 idx <- match(snps, gwas$SNP)
 ok  <- !is.na(idx)
 if (sum(ok) < 2L) return(NULL)

 z <- rep(NA_real_, length(snps))
 z[ok] <- gwas$Z[idx[ok]]

 a1g <- toupper(gwas$A1[idx[ok]])
 a1r <- a1_ref[ok]
 flip <- a1g != a1r
 z[ok][flip] <- -z[ok][flip]

 # drop palindromic
 a2g <- toupper(gwas$A2[idx[ok]])
 pal <- (a1g %in% c("A","T") & a2g %in% c("A","T")) |
        (a1g %in% c("C","G") & a2g %in% c("C","G"))
 z[ok][pal] <- NA

 z <- z[!is.na(z)]
 if (length(z) < 2L) NULL else z
}
