#' HARMONI model fit
#'
#' Thin alias for \code{\link{ms_bpwas}}.
#'
#' @inheritParams ms_bpwas
#' @return The same object returned by \code{\link{ms_bpwas}}.
#' @export
harmoni <- function(gwas_list,
                    protein_models,
                    plink_prefix,
                    control = ms_bpwas_control()) {
  ms_bpwas(
    gwas_list = gwas_list,
    protein_models = protein_models,
    plink_prefix = plink_prefix,
    control = control
  )
}

#' HARMONI control constructor
#'
#' Thin alias for \code{\link{ms_bpwas_control}}.
#'
#' @inheritParams ms_bpwas_control
#' @return The same control object returned by \code{\link{ms_bpwas_control}}.
#' @export
harmoni_control <- function(...) {
  ms_bpwas_control(...)
}
