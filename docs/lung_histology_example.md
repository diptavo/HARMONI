# Lung Histology Example

This example is adapted from the existing lung histology sensitivity workflow.
It is intended for users who already have a table of feature-level TWAS or PWAS
Z statistics across a pooled lung cancer analysis and histology-specific
analyses.

The example does not require PLINK or prediction weights because it starts from
already-computed feature-level Z statistics.

## Input Table

The script expects a tab-delimited table with these columns:

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

Column meaning:

- `ID`: feature identifier, such as a gene, transcript, protein, or model ID
- `TWAS.Z_lungoverall`: pooled lung cancer Z statistic
- `TWAS.Z_luad`: lung adenocarcinoma Z statistic
- `TWAS.Z_sqc`: squamous-cell carcinoma Z statistic
- `TWAS.Z_sclc`: small-cell lung cancer Z statistic
- `BH_FDR_*`: baseline single-trait Benjamini-Hochberg FDR values

The script removes rows with missing Z statistics.

## Run the Demo

From the repository root:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --out_dir=example_output/lung_histology_demo
```

With no `--input`, the script generates a deterministic demonstration table.
This confirms that the package and example workflow run end to end.

## Run on an Already-Generated Lung Table

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --input=/path/to/lung_histology_twas_table.tsv \
  --out_dir=example_output/lung_histology_real
```

On Biowulf:

```bash
module load R
Rscript examples/lung_histology_from_twas_table.R \
  --input=/path/to/lung_histology_twas_table.tsv \
  --out_dir=/data/Dutta_lab/HARMONI/example_output/lung_histology_real
```

If the R module cannot be loaded on the login node, submit a batch job:

```bash
sbatch --mem=8g --time=01:00:00 \
  --output=lung_histology_example_%j.log \
  --wrap='module load R; cd /data/Dutta_lab/HARMONI; Rscript examples/lung_histology_from_twas_table.R --input=/path/to/lung_histology_twas_table.tsv --out_dir=example_output/lung_histology_real'
```

## Modeled Axes

The script evaluates three parameterizations:

1. `current_std`: pooled overall axis plus two hand-defined histology contrasts.
2. `orth_overall`: same conceptual axes, but orthogonalized under the empirical
   null covariance.
3. `subtype_shared`: shared axis defined by the average of the subtype-specific
   Z statistics, with subtype contrasts orthogonalized under the empirical null
   covariance.

The histology contrasts are:

```text
LUAD vs SQC
SCLC vs average non-SCLC subtype signal
```

## Null Correlation Estimation

The script first estimates an empirical covariance matrix from features with no
large absolute original Z statistic. It then transforms the original Z matrix
into the selected shared/heterogeneity basis and estimates a null correlation
matrix in the transformed space.

For small demo tables, if too few features pass the strict near-null threshold,
the script falls back to the lowest-ranked features by maximum absolute Z. This
fallback is intended for examples and smoke tests; production analyses should
have enough near-null features for stable null-correlation estimation.

## Output Files

For each parameterization, the script writes:

```text
<parameterization>_results.tsv
<parameterization>_summary.tsv
<parameterization>_config_counts.tsv
<parameterization>_loadings.tsv
<parameterization>.rds
```

It also writes:

```text
lung_histology_param_summary.tsv
lung_histology_Sigma0.tsv
```

## Interpreting the Summary

Important columns in `<parameterization>_summary.tsv`:

- `n_features`: number of analyzed features
- `n_null_cov`: number of features used to estimate covariance in original axes
- `n_null_transformed`: number of features used for transformed null correlation
- `vb_iter`: variational-Bayes iteration count
- `vb_converged`: convergence flag
- `pi0`: fitted global null prior probability
- `xi`: fitted axis-level non-null probabilities
- `bfdr_sig`: number of HARMONI discoveries at BFDR < 0.05
- `baseline_primary_sig`: number of discoveries from the selected baseline
- `msbpwas_only_primary`: HARMONI discoveries not detected by the selected
  baseline

Important columns in `<parameterization>_results.tsv`:

- `p_null`: posterior null probability
- `bfdr`: Bayesian FDR value
- `map_config`: most likely association configuration
- `pip_*`: posterior inclusion probability for each modeled axis

## Production Recommendations

Before interpreting a real lung analysis:

- inspect `lung_histology_Sigma0.tsv`
- confirm `vb_converged` is `TRUE`
- inspect transformed-axis loadings in `*_loadings.tsv`
- compare `bfdr_sig` against baseline single-trait discoveries
- review top HARMONI-only findings manually
- check whether discoveries are driven by shared or heterogeneous axes
