# Hyperaldosteronism GWAS workflow demonstration

This repository is a public-safe, synthetic-data demonstration of the summary-statistic harmonization and meta-analysis workflow used in a multi-biobank hyperaldosteronism genetics project.

It contains no participant-level data, controlled-access summary statistics, credentials, private paths, or unpublished numerical results. The example values are synthetic and must not be interpreted as scientific findings.

## What this demonstrates

- A manifest-driven separation of European and multi-ancestry inputs.
- Summary-statistic validation and conservative allele harmonization.
- Effect-direction and effect-allele-frequency flipping for swapped alleles.
- Conservative removal of palindromic or irreconcilable variants.
- Fixed-effect inverse-variance meta-analysis.
- Cochran's Q and I-squared heterogeneity diagnostics.
- Reproducible QC logs and deterministic outputs.

The production research used REGENIE for MVP/UK Biobank association testing, integrated FinnGen SAIGE and published French GWAS summary statistics, and used GWAMA for the final meta-analyses. The R code here reimplements the statistical core on synthetic data so that the workflow can run without restricted software or data.

This is a portfolio demonstration, not the primary research repository and not a release of the study's results.

## Run the demo

Requirements: R 4.1 or later. No external R packages are required.

```bash
bash tests/run_demo.sh
```

The command creates:

- `results/european_meta.tsv`
- `results/european_qc.tsv`
- `results/multi_ancestry_meta.tsv`
- `results/multi_ancestry_qc.tsv`

The same command runs automatically on every GitHub push and pull request.

## Input contract

Each cohort file must contain:

| Column | Meaning |
|---|---|
| `CHR` | Chromosome |
| `POS` | One-based position on a common genome build |
| `EA` | Effect allele |
| `OA` | Other allele |
| `BETA` | Effect on the log-odds or quantitative-trait scale |
| `SE` | Standard error of `BETA` |
| `P` | Cohort association P value |
| `EAF` | Effect-allele frequency |
| `N` | Variant-level or analysis sample size |

The manifest records cohort name, ancestry label, whether the input belongs in the European analysis, and the relative file path.

## Important limitations

- The demo aligns to the first valid allele representation at each chromosome-position key. Production work should use a locked reference allele convention and explicit build/variant crosswalks.
- All A/T and C/G variants are dropped conservatively. A production pipeline may retain some with trustworthy frequency-based resolution.
- The demo uses complete-case summary rows and does not replace cohort-specific genotype, imputation, sample, phenotype, or association-model QC.
- Fixed-effect pooling is shown for transparency; scientific use requires prespecified models and cohort-aware sensitivity analyses.

## Repository map

- `config/cohorts.example.tsv`: analysis manifest.
- `data/example/`: synthetic cohort summary statistics.
- `scripts/fixed_effect_meta.R`: validation, harmonization, pooling, and QC logging.
- `tests/run_demo.sh`: end-to-end smoke test.
- `docs/interview_refresher.md`: concise explanation of the workflow.
- `docs/publication_boundary.md`: what must never be committed from a controlled-access project.
