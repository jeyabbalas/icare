# Build and apply an absolute risk model

Estimates the absolute risk of a disease over an age interval for each
covariate (and/or SNP) profile. The model combines

- a model for the **relative risks** of the risk factors (log relative
  risks of the design-matrix columns of `model.formula`, and/or SNP odds
  ratios),

- a **reference dataset** of risk factors representing the population in
  which the model is applied,

- the population's age-specific **disease incidence rates**, and

- optionally, age-specific **competing mortality rates**.

The baseline hazard is calibrated so that the model reproduces the
population incidence rates, and risks of profiles with missing risk
factors are averaged coherently over the reference population.

## Usage

``` r
computeAbsoluteRisk(
  model.formula = NULL,
  model.cov.info = NULL,
  model.snp.info = NULL,
  model.log.RR = NULL,
  model.ref.dataset = NULL,
  model.ref.dataset.weights = NULL,
  model.disease.incidence.rates,
  model.competing.incidence.rates = NULL,
  model.bin.fh.name = NA,
  n.imp = 5,
  apply.age.start,
  apply.age.interval.length,
  apply.cov.profile = NULL,
  apply.snp.profile = NULL,
  return.lp = FALSE,
  return.refs.risk = FALSE,
  ...
)
```

## Arguments

- model.formula:

  A formula such as `Y ~ x1 + as.factor(x2) + x1:x3` describing the
  design matrix of the relative-risk model. The response is ignored.
  Leave `NULL` (with the other covariate arguments) for a SNP-only
  model.

- model.cov.info:

  A list with one element per covariate of `model.formula`, each a list
  with fields `name` (the covariate name), `type` (`"continuous"` or
  `"factor"`) and, for factors, `levels` (all levels) and optionally
  `ref` (the referent level).

- model.snp.info:

  A data frame with columns `snp.name`, `snp.odds.ratio` (per-allele
  odds ratio) and `snp.freq` (risk-allele frequency).

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

- model.disease.incidence.rates:

  Age-specific annual incidence rates: a two-column table (integer age,
  rate) or a three-column table (start age, end age, annual rate in that
  age band). Must cover consecutive ages including every year of age
  used by the risk intervals.

- model.competing.incidence.rates:

  Age-specific annual rates of competing mortality, in the same layouts.
  `NULL` means no competing mortality.

- model.bin.fh.name:

  Name of a binary (0/1) family-history covariate, used with
  `model.snp.info` to simulate genotypes given family history and to
  attenuate its log relative risk. `NA` when not applicable.

- n.imp:

  Number of simulated SNP genotype sets of the reference population (1
  to 20).

- apply.age.start:

  Age(s) at the start of the risk interval: a single value or one per
  profile. Fractional ages are allowed.

- apply.age.interval.length:

  Length(s) of the risk interval in years: a single value or one per
  profile. Fractional lengths are allowed.

- apply.cov.profile:

  A data frame of covariate profiles (columns matched by name; extra
  columns are ignored); missing values are allowed.

- apply.snp.profile:

  A matrix or data frame of genotypes (risk-allele counts 0, 1, 2 or
  `NA`), one column per SNP of `model.snp.info`, matched by name when
  the columns are named after all the SNPs and by position when they are
  named after none (naming only some is an error). If `NULL`, all SNPs
  are treated as missing.

- return.lp:

  If `TRUE`, also return the linear predictors.

- return.refs.risk:

  If `TRUE`, also return the risks of the reference population over the
  first profile's interval.

- ...:

  Not used. The removed argument `use.c.code` is accepted with a
  deprecation warning and ignored.

## Value

An object of class `"icare"`, a list with

- `risk`: one-column matrix (`Risk_Estimate`) of absolute risks;

- `details`: data frame with the interval (`Int_Start`, `Int_End`), the
  risk (`Risk_Estimate`) and the profile (genotypes, then covariates);

- `beta.used`: one-column matrix of the log relative risks used (SNPs
  first; the family-history coefficient after attenuation);

- `lps`: one-column matrix of linear predictors (if `return.lp`); for
  profiles with missing risk factors, the model-averaged linear
  predictor;

- `refs.risk`: risks of the reference rows (if `return.refs.risk`).

The attribute `"calibration"` holds the calibrated baseline hazard: a
data frame with the `age`, the population `incidence.rate`, the
`baseline.hazard` and the survivor-weighted mean relative risk
`expected.rr` at each age of the incidence table.

## Details

**Model.** The hazard of a person with risk factors \\z\\ at age \\t\\
is \\\lambda_0(t)\exp(\beta' x(z))\\. The absolute risk over \\\[a, a +
\tau)\\ is computed over whole years of age \\t\\, as in earlier
versions of iCARE: \$\$AR = \sum_t e_t \lambda_0(t) RR
\exp\\-\sum\_{u\<t} e_u (\lambda_0(u) RR + m(u))\\,\$\$ where \\RR =
\exp(\beta' x(z))\\, \\m\\ is the competing mortality rate and \\e_t\\
is the time the interval spends in the year of age \\\[t, t + 1)\\: 1
inside the interval for whole-number ages, and the corresponding
fraction of a year at fractional ends.

**Calibration.** The baseline hazard solves \\\lambda(t) = \lambda_0(t)
E\[RR \mid T \ge t\]\\, the expectation being over the (weighted)
reference population still free of the disease at age \\t\\, starting
from the youngest age of the incidence table.

**Missing risk factors.** For a profile with missing risk factors, the
risk is \\E\[AR(z_o, Z_m) \mid z_o\]\\: reference rows are grouped into
100 percentile bins of the linear predictor of the observed risk
factors, and the risk is the weighted mean, over the rows in the
profile's bin, of the risk of the profile's own observed risk factors
completed with each row's values of the missing ones. A profile with no
observed risk factor gets the population-average risk.

**SNPs.** SNP genotypes of the reference population are simulated under
Hardy-Weinberg equilibrium from the allele frequencies (`n.imp` sets),
conditionally on family history when `model.bin.fh.name` is given, in
which case the family-history log relative risk is attenuated by half
the variance of the polygenic score. Missing genotypes in
`apply.snp.profile` are handled like other missing risk factors. Results
depend on the random number generator; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R package to
build, validate and apply absolute risk models. *PLoS ONE* 15(2):
e0228198.
[doi:10.1371/journal.pone.0228198](https://doi.org/10.1371/journal.pone.0228198)

## See also

[`computeAbsoluteRiskSplitInterval()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRiskSplitInterval.md)
for models that change at a cut age,
[`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
to validate a model.

## Examples

``` r
data(bc_data, package = "iCARE")

# Covariates and 72 SNPs, with family-history adjustment
set.seed(50)
results <- computeAbsoluteRisk(
  model.formula = bc_model_formula,
  model.cov.info = bc_model_cov_info,
  model.snp.info = bc_72_snps,
  model.log.RR = bc_model_log_or,
  model.ref.dataset = ref_cov_dat,
  model.disease.incidence.rates = bc_inc,
  model.competing.incidence.rates = mort_inc,
  model.bin.fh.name = "famhist",
  apply.age.start = 50,
  apply.age.interval.length = 30,
  apply.cov.profile = new_cov_prof,
  apply.snp.profile = new_snp_prof,
  return.refs.risk = TRUE
)
results
#> Absolute risks of 3 profiles
#> 
#>   Int_Start Int_End Risk_Estimate
#> 1        50      80    0.10173346
#> 2        50      80    0.09031505
#> 3        50      80    0.16899913
#> Profiles and risks: `$details`; risks: `$risk`.
summary(results)
#>  Risk_Estimate    
#>  Min.   :0.09032  
#>  1st Qu.:0.09602  
#>  Median :0.10173  
#>  Mean   :0.12035  
#>  3rd Qu.:0.13537  
#>  Max.   :0.16900  

# A SNP-only model
snp_only <- computeAbsoluteRisk(
  model.snp.info = bc_72_snps,
  model.disease.incidence.rates = bc_inc,
  model.competing.incidence.rates = mort_inc,
  apply.age.start = 50,
  apply.age.interval.length = 30,
  apply.snp.profile = new_snp_prof
)
snp_only$risk
#>      Risk_Estimate
#> [1,]    0.08822113
#> [2,]    0.10509014
#> [3,]    0.13469787
```
