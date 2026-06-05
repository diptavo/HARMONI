# HARMONI

HARMONI is an R package for joint Bayesian analysis of pooled and
subtype-specific PWAS or TWAS association statistics. It was developed for
multi-subtype disease settings where a signal may be shared across disease
subtypes, restricted to one subtype, or heterogeneous across subtypes.

The package name expands to **Hierarchical Association Resolution via
Multivariate Orthogonalized Non-Null Inference**.

## What HARMONI Does

Standard PWAS or TWAS analyses test one molecular feature against one trait at a
time. In cancer and other heterogeneous diseases, that approach can miss signals
that are diluted in an all-cases analysis or that differ across subtypes.

HARMONI jointly models the vector of association Z statistics for each feature
across a set of pooled and subtype-specific traits. For each feature, it returns:

- posterior probabilities over association configurations
- marginal posterior inclusion probabilities for each modeled axis
- posterior null probabilities
- maximum a posteriori configuration labels
- Bayesian false-discovery-rate summaries
- posterior effect estimates and uncertainty summaries
- fitted global hyperparameters and convergence diagnostics

The most direct workflow starts from an already-computed feature-by-trait table
of PWAS or TWAS Z statistics. This is the common handoff point after per-trait
PWAS/TWAS scans have already been run. HARMONI then models those Z statistics
jointly across user-defined shared and heterogeneity axes.

The full `harmoni()` / `ms_bpwas()` entry point is also available for users who
want HARMONI to compute PWAS Z statistics directly from GWAS summary statistics,
prediction weights, and a PLINK LD reference panel.

## Repository Layout

```text
.
├── DESCRIPTION                 # R package metadata
├── NAMESPACE                   # exported functions and S3 methods
├── R/                          # package implementation
├── examples/                   # runnable examples
│   └── lung_histology_from_twas_table.R
├── docs/                       # detailed user and developer documentation
│   ├── prerequisites.md
│   ├── lung_histology_example.md
│   └── github_release_checklist.md
├── .Rbuildignore
├── .gitignore
└── README.md
```

This directory is intended to be a standalone GitHub repository. From this
directory, `git init`, commit, and push directly to a remote repository.

## Prerequisites

### Required for Installing the Package

- R version 4.0 or newer.
- Standard R build tools for source-package installation.
  - Linux: system compilers and development headers appropriate for building R packages.
  - macOS: Xcode command line tools are usually enough.
  - Windows: Rtools matching your R version.
- No mandatory CRAN package dependencies for the core package namespace.

### Required for Full PWAS/TWAS Runs from GWAS Summary Statistics

- PLINK 1.9 or PLINK 2.0 available on `PATH`.
- A matched LD reference panel in PLINK binary format:
  - `.bed`
  - `.bim`
  - `.fam`
- GWAS summary statistics with at least:
  - `SNP`
  - `A1`
  - `A2`
  - `BETA`
  - `SE`
  - `P`
  - `N`
- Protein or expression prediction models containing:
  - `protein_id`
  - `gene`
  - `chr`
  - `weights`, a data frame with `SNP`, `A1`, and `weight`
  - `cv_R2`
- Consistent genome build, variant identifiers, and allele coding between the
  GWAS files, prediction weights, and LD reference panel.

### Required for Simulation and Some Analysis Scripts

The package itself is intentionally light, but some simulation and workflow
scripts use additional R packages:

- `MASS`, for multivariate normal simulation
- `dplyr`, for scenario-grid manipulation in some simulation drivers

Install them when running those scripts:

```r
install.packages(c("MASS", "dplyr"))
```

### Running Larger Analyses

The examples can be run on a laptop or workstation. Larger PWAS/TWAS analyses
should be run wherever adequate memory, CPU time, and temporary disk space are
available. On a shared cluster, use the scheduler and resource requests
appropriate for that system. Keep logs and generated outputs outside the source
tree or in ignored directories such as `logs/`, `outputs/`, or
`example_output/`.

## Installation

### Install from a Local Checkout

From the parent directory of this repository:

```r
install.packages("HARMONI", repos = NULL, type = "source")
```

Or from anywhere:

```r
install.packages("/path/to/HARMONI", repos = NULL, type = "source")
```

### Build a Source Tarball

```bash
R CMD build HARMONI
R CMD INSTALL HARMONI_0.1.0.tar.gz
```

### Install from GitHub

After pushing this directory to GitHub:

```r
install.packages("remotes")
remotes::install_github("OWNER/HARMONI")
```

Replace `OWNER/HARMONI` with the actual GitHub organization or user and
repository name.

## Quick Start: Already-Computed TWAS/PWAS Z Statistics

Use this workflow when you already have one row per feature and one Z-statistic
column per trait or subtype. The lung histology example is the primary template:
it starts from a pooled lung cancer TWAS/PWAS column plus LUAD, SQC, and SCLC
histology-specific columns.

Required table shape:

```text
ID
TWAS.Z_lungoverall
TWAS.Z_luad
TWAS.Z_sqc
TWAS.Z_sclc
BH_FDR_lungoverall
BH_FDR_luad
BH_FDR_sqc
BH_FDR_sclc
```

Run HARMONI on an already-generated lung table:

```bash
cd /path/to/HARMONI

Rscript examples/lung_histology_from_twas_table.R \
  --input=/path/to/lung_histology_twas_table.tsv \
  --out_dir=example_output/lung_histology_real
```

For large tables, run the same command on a workstation or cluster job with
enough memory for the number of features and modeled axes.

The example writes one result set per contrast parameterization:

```text
current_std_results.tsv
current_std_summary.tsv
current_std_config_counts.tsv
current_std_loadings.tsv
current_std.rds

orth_overall_results.tsv
orth_overall_summary.tsv
orth_overall_config_counts.tsv
orth_overall_loadings.tsv
orth_overall.rds

subtype_shared_results.tsv
subtype_shared_summary.tsv
subtype_shared_config_counts.tsv
subtype_shared_loadings.tsv
subtype_shared.rds

lung_histology_param_summary.tsv
lung_histology_Sigma0.tsv
```

The top-level comparison file is:

```text
example_output/lung_histology_real/lung_histology_param_summary.tsv
```

Use that file to compare convergence, discovery counts, baseline overlap, and
HARMONI-only discoveries across contrast definitions.

### Contrast Design for Already-Computed Z Tables

For an already-computed Z table, the analysis has two explicit steps:

1. Arrange the original Z matrix with one column per observed trait.
2. Choose a contrast/loading matrix that maps those observed traits into
   modeled HARMONI axes.

In the lung example, the observed coefficients are:

```text
overall
luad
sqc
sclc
```

The example evaluates three biologically motivated parameterizations:

- `current_std`: uses the pooled all-lung Z statistic as the shared axis, then
  adds LUAD-vs-SQC and SCLC-vs-non-SCLC heterogeneity axes.
- `orth_overall`: uses the same conceptual axes, but orthogonalizes them under
  the empirical null covariance before fitting HARMONI.
- `subtype_shared`: uses the average of LUAD, SQC, and SCLC as the shared axis,
  then adds orthogonalized histology heterogeneity axes.

The core loading definitions are:

```r
coef_names <- c("overall", "luad", "sqc", "sclc")

shared_overall <- c(1, 0, 0, 0)
shared_subtype_average <- c(0, 1 / 3, 1 / 3, 1 / 3)

luad_vs_sqc <- c(0, 1, -1, 0)
sclc_vs_nsclc <- c(0, -0.5, -0.5, 1)

B <- cbind(
  shared = shared_overall,
  het_luad_vs_sqc = luad_vs_sqc,
  het_sclc_vs_nsclc = sclc_vs_nsclc
)
rownames(B) <- coef_names
```

The transformed HARMONI Z matrix is:

```r
z_harmoni <- z_original %*% Q
```

where `Q` is either the column-normalized version of `B` or an empirical-null
orthogonalized version of `B`. Orthogonalization is useful when axes are not
independent under the null, because HARMONI then fits configuration probabilities
on axes that are less correlated by construction.

To run a different contrast analysis, add another entry to the `specs` list in
`examples/lung_histology_from_twas_table.R`. For example, to compare the pooled
overall axis with a single LUAD-vs-all-other heterogeneity axis:

```r
specs[["luad_vs_other"]] <- list(
  shared = c(1, 0, 0, 0),
  het = cbind(het_1 = c(0, 1, -0.5, -0.5)),
  primary = "union_all",
  mode = "orth"
)
```

General rules for contrast construction:

- the vector length must equal the number of original Z columns
- a shared axis should encode the pooled or average disease signal you want to
  preserve
- heterogeneity axes should sum to approximately zero across the subtype columns
  when they are intended to compare subtypes
- use `mode = "orth"` when comparing multiple non-independent contrasts
- inspect `*_loadings.tsv` after each run to confirm the transformed axes match
  the intended scientific question

Run the built-in demo only as a smoke test when a real lung table is not
available:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --out_dir=example_output/lung_histology_demo
```

See `docs/lung_histology_example.md` for the full column specification, contrast
details, output files, and interpretation guidance.

## Advanced Workflow: Full PWAS/TWAS from GWAS Summary Statistics

Use `harmoni()` when you want the package to compute PWAS Z statistics from
subtype GWAS summary statistics, prediction weights, and an LD reference panel.

```r
library(HARMONI)

fit <- harmoni(
  gwas_list = list(
    mixture = "data/all_cases_vs_controls.tsv.gz",
    subtype1 = "data/subtype1_vs_controls.tsv.gz",
    subtype2 = "data/subtype2_vs_controls.tsv.gz"
  ),
  protein_models = "data/protein_prediction_models.RData",
  plink_prefix = "data/1000G_EUR",
  control = harmoni_control(
    method = "vb",
    prior_type = "scale_mixture",
    overlap_correction = "ldsc",
    subtype_prevalences = c(subtype1 = 0.65, subtype2 = 0.25)
  )
)

summary(fit)
head(fit$bfdr)
head(fit$pip)
```

The residual non-modeled subtype is computed from the mixture and named
subtypes when subtype prevalences are supplied.

## Main User-Facing Functions

- `harmoni()` / `ms_bpwas()`: full analysis from GWAS files, prediction models,
  and PLINK LD reference data.
- `harmoni_control()` / `ms_bpwas_control()`: create a control object specifying
  inference method, prior, convergence settings, and overlap correction.
- `build_configs()`: enumerate association configurations.
- `fit_vb()`, `fit_em()`, `fit_mcmc()`: fit lower-level models to an already
  constructed Z matrix and null correlation matrix.
- `compute_derived()`: compute PIPs, posterior null probabilities, MAP
  configurations, conditional p-values, and BFDR from a fitted model.
- `construct_shared_het_z()` and heterogeneity helpers: construct and test
  shared and heterogeneous association axes.

## Input Requirements for `harmoni()`

### `gwas_list`

A named list of GWAS files or data frames. It must include `mixture` and at
least two subtype entries for a multi-subtype analysis.

Each file or data frame should contain:

```text
SNP  A1  A2  BETA  SE  P  N
```

Recommended additional checks before running:

- remove duplicated variants or make duplicate handling explicit
- ensure effect alleles are harmonized to the prediction weights
- ensure all traits use the same genome build
- check that sample sizes and case definitions match the intended subtype model

### `protein_models`

Either a path to an `.RData` object or a list of model objects. Each model should
include:

```text
protein_id
gene
chr
weights: data.frame(SNP, A1, weight)
cv_R2
```

### `plink_prefix`

Prefix to a PLINK binary reference panel:

```text
/path/to/reference
/path/to/reference.bed
/path/to/reference.bim
/path/to/reference.fam
```

Pass only the prefix, not the file extension.

## Output Overview

A fitted `harmoni()` object has class `msbpwas` and includes:

- `posterior_configs`: posterior probability for each association configuration
- `pip`: marginal posterior inclusion probabilities
- `p_null`: posterior null probability
- `map_config`: maximum a posteriori configuration
- `p_conditional`: conditional frequentist p-value given the MAP configuration
- `effect_sizes`: model-averaged posterior mean effects
- `effect_ses`: posterior standard deviations
- `bfdr`: Bayesian FDR values
- `z_matrix`: feature-by-subtype Z-statistic matrix
- `R_overlap`: estimated null correlation matrix
- `configs`: configuration matrix
- `hyperparams`: fitted global hyperparameters
- `convergence`: iteration count and convergence flag

## Recommended Analysis Workflow

1. Harmonize GWAS files, prediction weights, and LD reference variants.
2. Run a small subset or smoke test first.
3. Inspect `R_overlap` to confirm the inferred null correlation is plausible.
4. Check convergence diagnostics before interpreting discoveries.
5. Use BFDR thresholds for discovery reporting.
6. Compare HARMONI discoveries with single-trait baseline analyses.
7. Report both shared and heterogeneous configuration summaries.

## Development Checks Before Publishing

From the parent directory:

```bash
R CMD build HARMONI
R CMD check HARMONI_0.1.0.tar.gz
```

If you run checks on a cluster, submit the same `R CMD build` and
`R CMD check` commands through your site's scheduler.

## Citation

If this package is used in a manuscript, cite the associated HARMONI/MS-BPWAS
method manuscript or preprint once available. Add the final citation to this
README and to `inst/CITATION` before public release.
