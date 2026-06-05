\name{HARMONI-exports}
\alias{HARMONI-exports}
\alias{bfdr_from_pnull}
\alias{build_configs}
\alias{clamp}
\alias{classify_heterogeneity}
\alias{combine_heterogeneity}
\alias{compute_derived}
\alias{compute_marglik}
\alias{config_prior}
\alias{construct_shared_het_z}
\alias{contrast_heterogeneity}
\alias{directional_heterogeneity}
\alias{estimate_null_correlation}
\alias{fit_em}
\alias{fit_mcmc}
\alias{fit_vb}
\alias{harmoni}
\alias{harmoni_control}
\alias{make_psd}
\alias{marglik_horseshoe}
\alias{marglik_scale_mixture}
\alias{marglik_soft_ss}
\alias{ms_bpwas}
\alias{ms_bpwas_control}
\alias{msg}
\alias{posterior_heterogeneity}
\alias{precompute_R}
\alias{precompute_Reff}
\alias{prepare_config_cache}
\alias{quantitative_heterogeneity}
\alias{run_exact_estep}
\alias{solve_safe}
\alias{test_het_direct}
\alias{test_het_shared}
\alias{test_heterogeneity}
\title{HARMONI exported functions}
\description{
Consolidated reference page for the exported HARMONI API. HARMONI fits
multi-axis Bayesian PWAS/TWAS models for pooled and subtype-specific association
statistics. The high-level entry point is \code{harmoni()}, with the original
implementation name \code{ms_bpwas()} retained as an alias.
}
\details{
The exported functions fall into four groups.

\strong{High-level workflow:}
\code{harmoni()}, \code{ms_bpwas()}, \code{harmoni_control()}, and
\code{ms_bpwas_control()} run the full analysis from GWAS summary statistics,
prediction models, and a PLINK LD reference panel.

\strong{Low-level fitting:}
\code{build_configs()}, \code{estimate_null_correlation()}, \code{fit_vb()},
\code{fit_em()}, \code{fit_mcmc()}, and \code{compute_derived()} support
workflows that already have feature-by-axis Z statistics.

\strong{Likelihood and numerical helpers:}
\code{make_psd()}, \code{solve_safe()}, \code{clamp()}, \code{precompute_R()},
\code{precompute_Reff()}, \code{prepare_config_cache()},
\code{run_exact_estep()}, \code{config_prior()}, \code{marglik_scale_mixture()},
\code{marglik_soft_ss()}, \code{marglik_horseshoe()},
\code{compute_marglik()}, \code{bfdr_from_pnull()}, and \code{msg()} are
used internally and are exported for reproducible lower-level workflows.

\strong{Heterogeneity helpers:}
\code{construct_shared_het_z()}, \code{test_heterogeneity()},
\code{posterior_heterogeneity()}, \code{quantitative_heterogeneity()},
\code{directional_heterogeneity()}, \code{contrast_heterogeneity()},
\code{combine_heterogeneity()}, \code{classify_heterogeneity()},
\code{test_het_direct()}, and \code{test_het_shared()} support shared and
heterogeneous signal summaries.
}
\value{
Return values depend on the function. High-level fits return an object of class
\code{"msbpwas"}. Lower-level fitting functions return fit lists containing
posterior probabilities, hyperparameters, and convergence diagnostics. See the
README and files in \code{docs/} for workflow-level examples and expected output
files.
}
\seealso{
\code{\link{harmoni}}, \code{\link{ms_bpwas}}, \code{\link{harmoni_control}},
\code{\link{ms_bpwas_control}}
}
\keyword{package}
