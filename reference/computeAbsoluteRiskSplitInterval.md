# Absolute risk with a model that changes at a cut age

Estimates absolute risks over an age interval when the relative-risk
model, the covariate profiles or the reference population change at a
cut age (for example, different log relative risks before and after
menopause). The interval \\\[a, a + \tau)\\ is split at the cut age
\\c\\; the first model gives the risk \\AR_1\\ over \\\[a, c)\\ and the
probability \\S_1\\ of reaching age \\c\\ alive and free of the disease,
the second model gives the risk \\AR_2\\ over \\\[c, a + \tau)\\, and
the risk over the whole interval is \$\$AR = AR_1 + S_1 AR_2.\$\$ With
the same model on both sides and a whole-number cut age, this equals the
risk from
[`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
over the whole interval for profiles without missing risk factors. (A
fractional cut age splits a year of age, within which the yearly risk
model is not exactly additive, so the two can then differ very
slightly.) (Earlier versions combined the parts as \\AR_1 + (1 - AR_1)
AR_2\\, which ignores deaths from competing causes before the cut age
and overstates the risk.)

## Usage

``` r
computeAbsoluteRiskSplitInterval(
  apply.age.start,
  apply.age.interval.length,
  apply.cov.profile,
  model.formula,
  model.disease.incidence.rates,
  model.log.RR,
  model.ref.dataset,
  model.ref.dataset.weights = NULL,
  model.cov.info,
  model.competing.incidence.rates = NULL,
  return.lp = FALSE,
  apply.snp.profile = NULL,
  model.snp.info = NULL,
  model.bin.fh.name = NULL,
  cut.time = NULL,
  apply.cov.profile.2 = NULL,
  model.formula.2 = NULL,
  model.log.RR.2 = NULL,
  model.ref.dataset.2 = NULL,
  model.ref.dataset.weights.2 = NULL,
  model.cov.info.2 = NULL,
  model.bin.fh.name.2 = NULL,
  n.imp = 5,
  return.refs.risk = FALSE,
  ...
)
```

## Arguments

- apply.age.start:

  Age(s) at the start of the risk interval: a single value or one per
  profile. Fractional ages are allowed.

- apply.age.interval.length:

  Length(s) of the risk interval in years: a single value or one per
  profile. Fractional lengths are allowed.

- apply.cov.profile:

  A data frame of covariate profiles (columns matched by name; extra
  columns are ignored); missing values are allowed.

- model.formula:

  A formula such as `Y ~ x1 + as.factor(x2) + x1:x3` describing the
  design matrix of the relative-risk model. The response is ignored.
  Leave `NULL` (with the other covariate arguments) for a SNP-only
  model.

- model.disease.incidence.rates:

  Age-specific annual incidence rates: a two-column table (integer age,
  rate) or a three-column table (start age, end age, annual rate in that
  age band). Must cover consecutive ages including every year of age
  used by the risk intervals.

- model.log.RR:

  Log relative risks of the design-matrix columns (no intercept), as a
  named vector or a one-column matrix with row names; the names must be
  the column names that `model.formula` generates.

- model.ref.dataset:

  A data frame of risk factors for a sample representing the population,
  with no missing values. Columns are matched to the covariates by name.

- model.ref.dataset.weights:

  Optional non-negative sampling weights of the rows of
  `model.ref.dataset` (normalized internally).

- model.cov.info:

  A list with one element per covariate of `model.formula`, each a list
  with fields `name` (the covariate name), `type` (`"continuous"` or
  `"factor"`) and, for factors, `levels` (all levels) and optionally
  `ref` (the referent level).

- model.competing.incidence.rates:

  Age-specific annual rates of competing mortality, in the same layouts.
  `NULL` means no competing mortality.

- return.lp:

  If `TRUE`, also return the linear predictors.

- apply.snp.profile:

  A matrix or data frame of genotypes (risk-allele counts 0, 1, 2 or
  `NA`), one column per SNP of `model.snp.info`, matched by name when
  the columns are named after all the SNPs and by position when they are
  named after none (naming only some is an error). If `NULL`, all SNPs
  are treated as missing.

- model.snp.info:

  A data frame with columns `snp.name`, `snp.odds.ratio` (per-allele
  odds ratio) and `snp.freq` (risk-allele frequency).

- model.bin.fh.name:

  Name of the binary family-history covariate (see
  [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md));
  `NULL` or `NA` when not applicable.

- cut.time:

  Cut age(s) at which the model changes: a single value or one per
  profile.

- apply.cov.profile.2, model.formula.2, model.log.RR.2,
  model.ref.dataset.2, model.ref.dataset.weights.2, model.cov.info.2,
  model.bin.fh.name.2:

  The model after the cut age (see Details).

- n.imp:

  Number of simulated SNP genotype sets of the reference population (1
  to 20).

- return.refs.risk:

  If `TRUE`, also return the risks of each model's reference population
  over the first profile's part of the interval.

- ...:

  Not used. The removed argument `use.c.code` is accepted with a
  deprecation warning and ignored.

## Value

An object of class `"icareSplit"`, a list with

- `risk`: one-column matrix (`Risk_Estimate`) of absolute risks;

- `details`: data frame with the start of the interval (`Int_Start`),
  the age at which the second model starts (`Cut_Age`), the end of the
  interval (`Int_End`), the covariate profile and the risk;

- `lps.1`, `lps.2`: linear predictors of both models (if `return.lp`);

- `refs.risk.1`, `refs.risk.2`: reference risks of both models (if
  `return.refs.risk`).

Without a cut age and second model, the result of
[`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
with class `"icareSplit"`.

## Details

Arguments ending in `.2` describe the model after the cut age; each one
left `NULL` is taken from its counterpart before the cut age (the
reference weights are reused with the reference dataset). The disease
and competing incidence rates, the SNP information and `n.imp` are
shared. If `cut.time` is outside a profile's interval, that profile's
risk comes from one model only. If no `.2` argument and no `cut.time` is
given, this is
[`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
over the whole interval.

For profiles with missing risk factors, \\AR_1\\, \\S_1\\ and \\AR_2\\
are each averaged over the reference population of their model (the two
models may have different reference datasets, so the missing values
cannot be completed jointly). The result can therefore differ slightly
from the risk over the whole interval even when both models are the
same; in the bundled example the relative difference is below 0.1%.

## See also

[`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)

## Examples

``` r
data(bc_data, package = "iCARE")

# Risk from age 40 to 80 with different log relative risks after age 50
results <- computeAbsoluteRiskSplitInterval(
  model.formula = bc_model_formula,
  model.cov.info = bc_model_cov_info,
  model.log.RR = bc_model_log_or,
  model.log.RR.2 = bc_model_log_or_post_50,
  model.ref.dataset = ref_cov_dat,
  model.ref.dataset.2 = ref_cov_dat_post_50,
  model.disease.incidence.rates = bc_inc,
  model.competing.incidence.rates = mort_inc,
  apply.age.start = 40,
  apply.age.interval.length = 40,
  cut.time = 50,
  apply.cov.profile = new_cov_prof
)
results
#> Absolute risks of 3 profiles
#> 
#>   Int_Start Cut_Age Int_End Risk_Estimate
#> 1        40      50      80    0.12393471
#> 2        40      50      80    0.08796538
#> 3        40      50      80    0.13718498
#> Profiles and risks: `$details`; risks: `$risk`.
```
