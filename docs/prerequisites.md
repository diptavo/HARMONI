# Prerequisites for Running HARMONI

This document lists the software, data, and analysis prerequisites needed for
HARMONI. The core R package has few software dependencies, but full PWAS/TWAS
analysis requires harmonized genetic inputs and an LD reference panel.

## 1. R and Build Tools

Required:

- R 4.0 or newer
- a working source-package build environment
- write access to an R library path

Tested environment:

- Biowulf compute node
- R 4.5.2
- HARMONI 0.1.0

Platform notes:

- Linux: use system compilers or the compiler toolchain loaded by the R module.
- macOS: install Xcode command line tools.
- Windows: install Rtools matching the installed R version.

Check your R version:

```r
getRversion()
```

Check library paths:

```r
.libPaths()
```

## 2. R Package Dependencies

The core HARMONI package namespace is designed to use base R functionality.

Optional packages used by simulation and workflow scripts:

- `MASS`
- `dplyr`

Install optional packages:

```r
install.packages(c("MASS", "dplyr"))
```

## 3. External Command-Line Tools

Full PWAS/TWAS runs from GWAS summary statistics require PLINK for LD
calculation.

Supported:

- PLINK 1.9
- PLINK 2.0

Check availability:

```bash
plink --version
```

or:

```bash
plink2 --version
```

The PLINK executable must be on `PATH` unless the workflow script exposes a
custom executable path.

## 4. LD Reference Panel

HARMONI expects a PLINK binary reference panel with files:

```text
reference.bed
reference.bim
reference.fam
```

Use the prefix without extension in R:

```r
plink_prefix <- "/path/to/reference"
```

The LD reference should match the GWAS ancestry as closely as possible. For
cross-ancestry or multi-ancestry analyses, prepare ancestry-appropriate
references and document how they were selected.

## 5. GWAS Summary Statistics

Each GWAS file or data frame should contain:

```text
SNP
A1
A2
BETA
SE
P
N
```

Strongly recommended:

- consistent genome build across all traits
- consistent variant IDs across traits and prediction weights
- no ambiguous or duplicated variants unless handled explicitly
- harmonized effect allele coding
- well-documented case/control definitions
- sample size columns that reflect the analyzed trait

For multi-subtype disease models, provide:

- one pooled or mixture GWAS
- one GWAS per modeled subtype
- subtype prevalence estimates if the residual subtype is to be inferred

## 6. Prediction Models

`protein_models` may be a list or an `.RData` file containing model objects.
Each model should include:

```text
protein_id
gene
chr
weights
cv_R2
```

The `weights` data frame should include:

```text
SNP
A1
weight
```

The SNP IDs and effect alleles must be compatible with the GWAS and LD
reference panel.

## 7. Already-Computed Z-Statistic Workflows

Some workflows already have feature-level TWAS or PWAS Z statistics. Those can
use lower-level functions such as `build_configs()`, `fit_vb()`, and
`compute_derived()` directly.

Required inputs for this mode:

- a numeric feature-by-axis Z matrix
- a positive semidefinite null correlation matrix among axes
- a configuration matrix from `build_configs()`

The lung histology example estimates the null correlation from near-null
features in an already-computed TWAS table.

## 8. Biowulf Requirements

Use compute nodes for R jobs. Biowulf may block `module load R` on login nodes.

Interactive:

```bash
sinteractive --mem=8g --time=02:00:00
module load R
```

Batch:

```bash
sbatch --mem=8g --time=02:00:00 --wrap='module load R; Rscript script.R'
```

For large runs:

- choose memory based on the number of features and modeled axes
- allocate local scratch if using large temporary files
- write logs to a dedicated `logs/` directory
- run a smoke test before submitting a full grid

## 9. Minimum Smoke Tests

After installation:

```r
library(HARMONI)
ctrl <- harmoni_control(method = "vb")
print(ctrl$method)
```

Run the included lower-level example:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --out_dir=example_output/lung_histology_demo
```

Expected files:

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
