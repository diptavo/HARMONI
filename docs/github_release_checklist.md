# GitHub Release Checklist

Use this checklist before pushing HARMONI to a public or private GitHub
repository.

## Repository Hygiene

- Confirm no private data files are included.
- Confirm no large result files, logs, `.rds` outputs, or compressed GWAS files
  are included.
- Confirm `.gitignore` excludes generated outputs.
- Confirm `.Rbuildignore` excludes GitHub-only and local-development files from
  the R source package.
- Confirm macOS resource-fork files named `._*` have been removed.

## Package Checks

Run from the parent directory:

```bash
R CMD build HARMONI
R CMD check HARMONI_0.1.0.tar.gz
```

On Biowulf, run these commands on a compute node:

```bash
sbatch --mem=8g --time=01:00:00 \
  --output=harmoni_check_%j.log \
  --wrap='module load R; cd /data/Dutta_lab; R CMD build HARMONI; R CMD check HARMONI_0.1.0.tar.gz'
```

## Example Check

From the repository root:

```bash
Rscript examples/lung_histology_from_twas_table.R \
  --out_dir=example_output/lung_histology_demo
```

Confirm that `lung_histology_param_summary.tsv` exists and that all
`vb_converged` values are `TRUE`.

## Suggested First Commit

```bash
git init
git add .
git status --short
git commit -m "Initial HARMONI package release"
git branch -M main
git remote add origin git@github.com:OWNER/HARMONI.git
git push -u origin main
```

Replace `OWNER/HARMONI` with the target repository.
