# Interview refresher

## Short explanation

The workflow makes cohort summary statistics comparable before pooling. It validates required fields, harmonizes effect alleles and direction, removes unresolved variants, runs separate European and multi-ancestry fixed-effect inverse-variance meta-analyses, and reports Cochran's Q and I-squared for heterogeneity.

## Why harmonization matters

A beta is meaningful only relative to its effect allele, genome build, phenotype, scale, and population. A silent allele swap reverses the sign. A build mismatch can join different variants. A duplicated marker can be counted twice. The QC layer therefore protects the scientific identity of each association before any pooling occurs.

## Core fixed-effect equations

For cohort estimate `beta_i` and standard error `se_i`:

```text
w_i       = 1 / se_i^2
beta_meta = sum(w_i * beta_i) / sum(w_i)
se_meta   = sqrt(1 / sum(w_i))
Q         = sum(w_i * (beta_i - beta_meta)^2)
I2        = max(0, (Q - (k - 1)) / Q)
```

For a binary outcome, `exp(beta_meta)` is the pooled odds ratio.

## What to emphasize

- European and multi-ancestry analyses answer related but not identical questions.
- Heterogeneity diagnostics should trigger inspection, not automatic storytelling.
- Lead SNP, independent significant SNP, conditional signal, and genomic risk locus are not synonyms.
- A public portfolio should demonstrate code and reasoning without exposing controlled data or unpublished collaborator results.

