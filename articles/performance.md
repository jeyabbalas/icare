# Performance

This distribution of iCARE has no compiled code. The official version
offers a C implementation of its slowest steps (`use.c.code = 1`, the
default); here, the same computations are written as vectorized R, and
every workload of the benchmark below runs faster than the official
version with its C code.

## Benchmark

The script
[`tools/benchmarks/run-benchmarks.R`](https://github.com/jeyabbalas/icare/tree/main/tools/benchmarks)
times seven workloads on the example data, with the same inputs in every
version:

| Workload | Description |
|:---|:---|
| `risk_complete` | [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md): 10,000 complete covariate profiles, start ages 30-70, intervals of 1-15 years |
| `risk_missing` | [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md): 2,000 profiles with 15% of covariate values missing |
| `risk_snps` | [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md): covariates and 72 SNPs with family history, 1,000 profiles with 30% of genotypes missing |
| `split_interval` | [`computeAbsoluteRiskSplitInterval()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRiskSplitInterval.md): 2,000 profiles from age 40 to 70, model changing at 50 |
| `validation_cohort` | [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md): the full cohort (49,767 subjects), risks supplied |
| `validation_ncc` | [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md): the nested case-control study (5,240 subjects), risks supplied |
| `validation_internal` | [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md): 10,000 cohort subjects, risks computed from the model |

The validation workloads use the subjects whose follow-up ends before
age 85 (49,767 of the 50,000 cohort members and 5,240 of the 5,285
case-control subjects). The example mortality rates stop at age 85, and
the official version gives `NA` or bin-averaged risks for intervals that
end there (see
[`vignette("differences")`](https://jeyabbalas.github.io/icare/articles/differences.md)),
so its timings would not be comparable.

Versions compared: this distribution; the official Bioconductor iCARE
1.40.0 with its C code (`use.c.code = 1`) and without it
(`use.c.code = 0`); and 1.33.1.9002, the patched development version
(with its C code). Times are the median elapsed seconds of three runs
(one run when the first takes longer than 20 seconds), on an Apple M5
Pro with R 4.6.2 and R’s reference BLAS.

| Workload | This distribution | Official 1.40.0, C code | Official 1.40.0, R code | 1.33.1.9002, C code | Speed-up vs official C |
|:---|---:|---:|---:|---:|---:|
| Risks, 10,000 complete profiles | 0.17 | 13.2 | 51.4 | 13.8 | 79.9× |
| Risks, 2,000 profiles with missing covariates | 0.94 | 2.97 | 7.63 | 3.06 | 3.2× |
| Risks with 72 SNPs, 1,000 profiles with missing genotypes | 7.43 | 12.2 | 28.1 | 11.6 | 1.6× |
| Split-interval risks, 2,000 profiles | 0.22 | 0.51 | 0.48 | 0.51 | 2.3× |
| Validation, full cohort (49,767), risks supplied | 0.048 | 0.96 | 0.82 | 1.18 | 20.0× |
| Validation, nested case-control (5,240), risks supplied | 0.010 | 0.21 | 0.21 | 0.22 | 20.9× |
| Validation, 10,000 subjects, risks computed | 0.17 | 0.39 | 0.38 | 0.40 | 2.3× |

The largest gains are for complete profiles, whose risks are computed
for all profiles and ages at once instead of profile by profile, and for
model validation, whose AUC is computed from sorted scores instead of a
matrix of all case and non-case pairs. The validation workloads with
supplied risks take a few hundredths of a second.

## Where the time goes

- **Complete profiles.** The cost is one matrix product for the linear
  predictors and a vectorized pass over the years of age of all
  intervals; it grows linearly with the number of profiles.
- **Missing risk factors.** Profiles are grouped by their pattern of
  missing risk factors. For each pattern, the reference population is
  sorted by the linear predictor of the observed risk factors (to form
  the percentile bins) and each profile’s risk is averaged over its
  donors. The cost therefore grows with the number of *distinct*
  missingness patterns rather than the number of profiles; profiles
  sharing a pattern cost little more than one.
- **SNPs.** With SNPs, the reference population holds `n.imp` copies of
  the reference dataset with simulated genotypes (5 × 14,137 rows in the
  example). Profiles with missing genotypes usually all have different
  patterns, so each costs one pass over that population: about 7 ms per
  profile here. Smaller `n.imp` makes this faster, at the price of more
  Monte Carlo error in the risks.

For large jobs, pass all profiles to one call: the model (including the
calibration of the baseline hazard and the simulated genotypes) is built
once per call.

## Reproducing the benchmark

From the root of a clone of the repository:

``` sh
# This distribution, installed into a library of its own
R CMD INSTALL --library=/path/to/lib-dev .
Rscript tools/benchmarks/run-benchmarks.R /path/to/lib-dev refactored

# The official version, installed with
# BiocManager::install("iCARE", lib = "/path/to/lib-official")
Rscript tools/benchmarks/run-benchmarks.R /path/to/lib-official official-c 1
Rscript tools/benchmarks/run-benchmarks.R /path/to/lib-official official-r 0

# A table of all results
Rscript tools/benchmarks/summarize-benchmarks.R
```

Each version must run in its own R process, because all versions are
named iCARE.
