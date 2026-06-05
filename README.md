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

The main package entry point for full PWAS/TWAS workflows is `harmoni()` or its
original alias `ms_bpwas()`. Lower-level exported functions are also available
for workflows that already have a matrix of Z statistics, such as the lung
histology example in `examples/`.

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

- R version 4.0 or newer. The Biowulf readiness test was run with R 4.5.2.
- Standard R build tools for source-package installation.
  - Linux/Biowulf: compiler toolchain provided by the R module.
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

### Biowulf Notes

Biowulf blocks `module load R` on the login node for some R modules. Run package
checks and analysis jobs on a compute node through `sinteractive` or `sbatch`.

Interactive example:

```bash
sinteractive --mem=8g --time=02:00:00
module load R
R
```

Batch example:

```bash
sbatch --mem=8g --time=02:00:00 --wrap='module load R; R CMD check HARMONI'
```

Biowulf may print a reminder to allocate `lscratch` for R jobs. For small
package installation and toy examples this is usually not necessary, but large
PWAS/TWAS jobs should use appropriate scratch space and job resources.

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

## Quick Start: Full PWAS/TWAS Workflow

Use `harmoni()` when you have subtype GWAS summary statistics, prediction
weights, and an LD reference panel.

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

## Quick Start: Already-Computed TWAS Z Statistics

If you already have a feature-by-trait table of TWAS or PWAS Z statistics, use
the lower-level functions. The lung histology example is based on the existing
lung workflow and demonstrates this mode.

Run a built-in demo:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --out_dir=example_output/lung_histology_demo
```

Run the same workflow on an already-generated lung TWAS table:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --input=/path/to/lung_histology_twas_table.tsv \
  --out_dir=example_output/lung_histology_real
```

See `docs/lung_histology_example.md` for the required columns, output files, and
interpretation details.

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

On Biowulf, run the check on a compute node:

```bash
sbatch --mem=8g --time=01:00:00 \
  --output=harmoni_check_%j.log \
  --wrap='module load R; cd /data/Dutta_lab; R CMD build HARMONI; R CMD check HARMONI_0.1.0.tar.gz'
```

## Citation

If this package is used in a manuscript, cite the associated HARMONI/MS-BPWAS
method manuscript or preprint once available. Add the final citation to this
README and to `inst/CITATION` before public release.
