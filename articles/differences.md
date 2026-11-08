# Differences from the official iCARE

This package is an **unofficial**, refactored distribution of iCARE. The
official package is maintained by its authors and distributed by
[Bioconductor](https://bioconductor.org/packages/iCARE/) (source:
[git.bioconductor.org/packages/iCARE](https://git.bioconductor.org/packages/iCARE),
mirrored at [github.com/bioc/iCARE](https://github.com/bioc/iCARE)),
with a development repository at
[github.com/parichoy/iCARE](https://github.com/parichoy/iCARE).

This article lists every difference from the official version 1.40.0
that a user can notice. The functions, their arguments and the names of
their results are unchanged, so existing scripts run as before; what
changes is that known defects are fixed, a few outputs are added or
corrected, and the figures are redrawn. Several fixes were first made in
iCARE 1.33.1.9002, a patched development version prepared by the
maintainer of this distribution; they are marked *(1.33.1.9002)*.

## At a glance

| Area | Change |
|:---|:---|
| Implementation | Vectorized R throughout; the C code and the `use.c.code` option are removed; faster than the official compiled code on every benchmark |
| Reference weights | Used correctly in the calibration of the baseline hazard *(1.33.1.9002)* |
| Ages | Fractional start ages and interval lengths are used exactly |
| Rates | Only the ages inside the intervals need rates: a missing rate where an interval ends no longer gives `NA` or bin-averaged risks |
| Incidence rates | Rates by age band (three-column tables) are read as yearly rates |
| Missing risk factors | The profile’s own observed risk factors are kept, as in the published method, instead of being replaced by the donors’ |
| Missing risk factors: intervals | Each profile’s risk uses its own interval and does not depend on the other profiles *(1.33.1.9002)* |
| Family history with SNPs | Its log relative risk is attenuated in every coding of the covariate |
| SNP profiles | Columns named after the SNPs are matched by name, not by position |
| SNP-only reference risks | `refs.risk` holds the risks of simulated people, no longer averages of five people |
| Split intervals | Competing mortality before the cut age is accounted for; reference weights are kept after it *(1.33.1.9002)* |
| AUC | Tied risk scores count 1/2 *(1.33.1.9002)*; O(n log n) time and O(n) memory |
| Risk groups | Ties, empty groups and cut-offs are handled, and the tests’ degrees of freedom follow the groups formed *(1.33.1.9002)* |
| Cases diagnosed at exit | Kept when ages have decimals, despite rounding error |
| Hosmer-Lemeshow test | Variance under the null hypothesis by default; the earlier statistic is an option |
| E/O in sampled data | The confidence interval accounts for the variability of E |
| Study incidence | Person-time rates by attained age, with confidence intervals |
| Results | Numeric group results with E/O and counts; new E/O per group with confidence intervals |
| Figures | Redrawn with ggplot2: square panels with their statistics, publication-ready |
| Messages | Clear errors and warnings; nothing printed to the console |

The sections below give the details.

## Implementation and interface

- **No compiled code.** The package is written in vectorized R. The C
  code of the official package and the `use.c.code` argument are
  removed; `computeAbsoluteRisk(use.c.code = 1)` still runs, with a
  deprecation warning, and the setting is ignored. Every workload of the
  benchmark is faster than the official version with its compiled code
  (see
  [`vignette("performance", package = "iCARE")`](https://jeyabbalas.github.io/icare/articles/performance.md)),
  and installing the package needs no compiler.
- **Dependencies.** `plotrix`, `gtools` and `Hmisc` are no longer
  needed; the figures use `ggplot2` (3.5.0 or later) and `patchwork`.
- **Model built once.** When
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
  computes risks for the study and for the reference population, the
  model (including the simulated SNP genotypes) is built once and
  shared, instead of once per call.
- **Methods.** Results have
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  results of
  [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
  (class `"icare"`) and
  [`computeAbsoluteRiskSplitInterval()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRiskSplitInterval.md)
  (class `"icareSplit"`) print a compact table and plot the distribution
  of the risks; results of
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
  (class `"icareValid"`) print a report that fits in 80 columns and plot
  with
  [`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md).
- **Calibration attribute.** The result of
  [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
  carries the calibrated baseline hazard
  (`attr(result, "calibration")`).
- **Messages.** Invalid inputs stop with errors that name the argument
  and the problem (for example, the ages missing from an incidence
  table). The official
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
  printed messages with [`print()`](https://rdrr.io/r/base/print.html);
  this version prints nothing and uses warnings and errors.

## Absolute risks: `computeAbsoluteRisk()`

- **Reference weights in the calibration** *(1.33.1.9002)*. The official
  version filled the subject-by-age weight matrix of the baseline-hazard
  calibration by row, which scrambled non-uniform
  `model.ref.dataset.weights` across subjects and ages. Every risk
  computed with reference weights was therefore miscalibrated (uniform
  weights, the default, were not affected). The baseline hazard now
  solves the calibration equation with each subject’s own weight; with
  integer weights it equals the hazard obtained by replicating each
  reference row as many times as its weight.
- **Fractional ages.** The official version needed whole-number start
  ages (a fractional start such as 60.5 was treated as 61 when mixed
  with whole ages, and stopped with an error otherwise) and counted any
  fraction of a year at the end of an interval as a whole year
  (intervals of 2.1 and 2.5 years gave the 3-year risk). Here, the time
  spent in each year of age is used exactly: a fraction of a year
  contributes that fraction of the year’s hazard.
- **Rates at the end of an interval.** The official version also read
  the rates at the age where an interval *ends*, an age the interval
  does not include, and multiplied them by zero. A missing (`NA`) rate
  there, such as the example mortality rates from age 85 on, therefore
  made the risks `NA` (`NA` times 0 is `NA`). Worse, when the profiles
  had different intervals and any of them ended at such an age, every
  profile with complete risk factors got `NA` and was then silently
  treated as if its risk factors were missing: its risk and linear
  predictor became averages over its percentile bin of the reference
  population. Rates are now needed only for the years of age inside the
  intervals, and a missing one is reported with the ages concerned.
- **Rates by age band.** Incidence and mortality tables with three
  columns (first age, last age, rate) were misread by the official
  version, which divided each band’s rate among the ages of the band:
  each year of a 5-year band got a fifth or a sixth of the yearly rate.
  Each year of age now gets the band’s rate.
- **Missing risk factors use the profile’s own interval**
  *(1.33.1.9002)*. For profiles with missing risk factors, the official
  version computed the reference risks of the percentile bins over the
  interval of the profiles as a whole (`unique(apply.age.start)` over
  all profiles), which was wrong whenever profiles had different
  intervals, and gave every profile with all risk factors missing the
  risk over the first profile’s interval. Each profile now gets the risk
  over its own interval, and a profile’s risk does not depend on which
  other profiles are computed with it.
- **Missing risk factors: the profile’s own observed values.** For a
  profile with missing risk factors, iCARE estimates the average risk
  over the values the missing risk factors take given the observed ones
  (Pal Choudhury et al. 2020, equation 3). The donors are the reference
  rows in the profile’s percentile bin of the linear predictor of the
  observed risk factors, and the paper averages “the risks obtained from
  multiple filled-in values” of the profile’s risk score. The official
  code instead averaged the donors’ *own* complete-data risks, which
  replaces the profile’s observed risk factors by the donors’. All
  profiles in a bin then got the same risk and linear predictor. The
  estimate missed its target wherever the bins are wide, as in the tails
  of the risk distribution, and differed from the complete-data risk
  even when the missing risk factors had no effect. The risk now
  combines the profile’s own observed risk factors with each donor’s
  values of the missing ones, including in interactions between an
  observed and a missing risk factor. The linear predictor is the
  profile’s observed part plus the donors’ average missing part. In the
  example model, we set 15% of each risk factor to missing, at random,
  for 2,000 women of the reference dataset. For the 1,751 women left
  with a missing risk factor, the 30-year risks changed by 0.3%
  (median), and by more than 10% for 8.5% of them. A woman with
  high-risk levels of all her risk factors and an unknown BMI now gets a
  30-year risk of 54%, against 58% with her BMI known; the official
  version gave her 25%. Risks of profiles with complete risk factors are
  unchanged.
- **Family history with SNPs.** The SNPs explain part of the familial
  risk, so the log relative risk of family history is attenuated when
  SNPs are in the model. The attenuation applies to the effect of family
  history in any coding of the covariate, including a factor whose
  reference level is “1” (design column `famhist0`), for which the
  official version stops with an error.
- **SNP columns matched by name.** The official version took the columns
  of `apply.snp.profile` in the order of `model.snp.info`, whatever
  their names. Columns named after the SNPs are now matched by name, so
  their order does not matter; unnamed columns are still taken in order,
  and naming only some of the SNPs is an error.
- **SNP-only reference risks.** For a model with SNPs only, the official
  `refs.risk` averaged the risks of the `n.imp` sets of simulated people
  position by position. Since the people of different sets are
  unrelated, each value was an average of five people’s risks, and the
  distribution was far too narrow (in the example, a standard deviation
  of 0.013 instead of 0.030 for the 30-year risk). `refs.risk` now holds
  the risks of the first set of 10,000 simulated people, a sample from
  the population distribution. For models with covariates and SNPs, each
  reference row’s risk is averaged over its own simulated genotype sets,
  as before.
- **Single profiles.** A single profile with a model of continuous
  covariates only failed in the official version (“non-conformable
  arguments”); it works here.
- **Design checks.** Formulas whose transformations give infinite or
  undefined values (for example `log(0)`) stop with an error naming the
  rows, instead of producing wrong risks.

## Split intervals: `computeAbsoluteRiskSplitInterval()`

- **Competing mortality before the cut age.** The official version
  combined the risks before and after the cut age as
  $`AR_1 + (1 - AR_1) AR_2`$, which assumes that everyone free of the
  disease at the cut age is alive. This version uses
  $`AR_1 + S_1 AR_2`$, where $`S_1`$ is the probability of reaching the
  cut age alive and free of the disease. The official formula overstates
  the risk, more so at older ages. With the same model on both sides and
  a whole-number cut age, the split risk now equals the risk over the
  whole interval.
- **Reference weights after the cut** *(1.33.1.9002)*. When no second
  reference dataset is given, the second model reuses the reference
  dataset *and its weights*; the official version dropped the weights.
- **`return.refs.risk`** *(1.33.1.9002)* is passed on when there is a
  single model.

## Model validation: `ModelValidation()`

### Corrections

- **Ties in the AUC** *(1.33.1.9002)*. The official AUC counted a case
  and a non-case with the same risk score as discordant (0 instead of
  1/2), which biases the AUC down when scores are tied, as they often
  are for questionnaire-based scores.
- **Exact risks for the study subjects.** Because of the defect with
  rates at the end of an interval (above), when
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
  computed the risks of the example studies, which follow women up to
  age 85, the official version replaced each subject’s risk and risk
  score by the average over her percentile bin: the 5,285 women of the
  nested case-control study had only 101 distinct risk scores. The risks
  and scores are now exact; on that study the AUC is 0.592 instead of
  0.587 (with the ties of the averaged scores also counted as 0).
- **Settings of the model object** *(1.33.1.9002)*.
  `model.ref.dataset.weights` and `n.imp` given in `iCARE.model.object`
  are used (the official version ignored them), settings left out take
  the defaults of
  [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md),
  and objects without `return.lp` or `return.refs.risk` no longer fail.
- **Supplied risks are kept** *(1.33.1.9002)*. When `predicted.risk` or
  `linear.predictor` is supplied, only the missing one is computed; the
  official version recomputed both, silently overwriting a supplied
  risk.
- **Onsets after exit** *(1.33.1.9002)*. A case whose onset is recorded
  after exit is censored at exit, now consistently in all statistics,
  and with a warning.
- **Cases diagnosed at exit.** Onsets are compared with the follow-up up
  to rounding error. When follow-up ends at diagnosis and ages are
  recorded as decimals, the follow-up computed as exit minus entry can
  fall just below the recorded onset time (in floating point, 70.94 -
  61.71 \< 9.23). The official version then turned the case into a
  non-case, silently: in simulated data with ages to two decimals, it
  dropped nearly half of such cases.
- **Risk-score groups** *(1.33.1.9002)*. Groups are formed at (weighted)
  quantiles of the risk score, with the repairs needed for tied scores:
  values tied at a quantile form their own group instead of stopping
  with “breaks are not unique”, and no subject is left without a group
  (in the official weighted branch, rounding could drop the subjects
  with the lowest or highest score). Empty groups are dropped, and the
  degrees of freedom of the Hosmer-Lemeshow test ($`G`$) and of the
  relative-risk test ($`G - 1`$) follow the number of groups $`G`$
  actually formed; the official version always used
  `number.of.percentiles`. Duplicated `linear.predictor.cutoffs` are
  merged, and cut-offs outside the range of the scores are ignored.
- **Reference population.** `reference.entry.age` and
  `reference.exit.age` may give one window per reference row (the
  official version used the first row’s window for all), the reference
  weights are used, and with SNPs the reference risk score is each row’s
  own linear predictor averaged over its simulated genotype sets (the
  official version computed it as if every SNP were missing).
- **Study incidence.** The official study incidence at age $`a`$ was the
  proportion of subjects, among those who entered at least a year before
  $`a`$ and were followed past $`a`$, with onset in $`[a, a + 1)`$. It
  ignored the first year of follow-up and follow-up ending within the
  year. The study incidence is now the person-time rate: onsets at each
  year of attained age divided by the person-years lived at that age, on
  the scale of the population incidence rates, with 95% confidence
  intervals (exact Poisson intervals in a full cohort, Wald intervals
  with a sampling term in a sampled study) and the numbers of cases and
  person-years. `Study_Incidence_Rate` is a data frame whose first two
  columns, `ages` and `study_incidence`, are those of the official
  matrix, followed by `lower`, `upper`, `cases` and `person_years`.
- **p-values** are computed in the upper tail, so very small p-values
  are reported instead of 0, and the `htest` objects describe their data
  (“observed and expected risks in 10 risk-score groups”).

### Changed statistics

- **Hosmer-Lemeshow variance.** The official statistic divided each
  group’s squared difference by the variance of the *observed* risk,
  $`O_g(1 - O_g)/N_g`$. With few cases per group that test rejects
  well-calibrated models too often, and a group without cases makes it
  infinite. In simulations with perfectly calibrated risks and 10 groups
  of about 10 cases, it rejected 21.8% of the time at the 5% level. The
  default is now the variance *under the null hypothesis* of
  calibration, $`E_g(1 - E_g)/N_g`$, as in the classical Hosmer-Lemeshow
  test; it rejected 5.6% of the time in the same simulations (see “How
  well the intervals and tests work” in
  [`vignette("model-validation")`](https://jeyabbalas.github.io/icare/articles/model-validation.md)).
  `hosmer.lemeshow.variance = "observed"` gives the official statistic.
- **E/O confidence interval in sampled data.** In a nested case-control
  or other sampled study, $`E`$ is a weighted sample mean, so it varies
  with the sample and covaries with $`O`$. The official interval treated
  $`E`$ as known. The variance of $`\log(E/O)`$ now includes both,
  through the linearization $`u_i = (p_i/E - y_i/O)/N`$; depending on
  the data, the interval can be wider or narrower. Intervals in a full
  cohort are unchanged.

### New results and checks

- `Category_Results` has numeric columns (the official ones were
  character strings) and adds each group’s E/O ratio with its 95%
  confidence interval and its numbers of subjects and cases. The E/O
  ratios are also returned as
  `Category_Specific_Expected_to_Observed_Ratio`.
- `Risk_Prediction_Interval` is a single label, e.g. `"5 years"`,
  `"Observed Followup"` or, for subject-specific intervals,
  `"Subject-specific (2 to 5 years)"`.
- `dataset` and `model.name` are stored with the results and used by
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md).
- Input checks: an error when no interval is given (neither
  `predicted.risk.interval` nor `total.followup.validation = TRUE`);
  cases must have a finite time of onset, while non-cases may have `Inf`
  or `NA`; unknown names in `iCARE.model.object` are reported with a
  warning; the settings `use.c.code`, `return.lp` and `return.refs.risk`
  of earlier model objects are accepted and ignored.

## Figures: `plotModelValidation()`

The official figure drew five base-graphics panels, with flattened
calibration plots and the statistics in a separate text panel, and
returned `NULL`. The new figure

- is drawn with ggplot2 and patchwork and returned as an object that can
  be modified and saved with
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html);
- has **square** panels, and the calibration panels have the same limits
  on both axes, so the line of perfect calibration runs at 45 degrees;
- shows each panel’s statistics with it: E/O and the Hosmer-Lemeshow
  test with the absolute-risk calibration, the relative-risk test with
  the relative-risk calibration, the AUC with the discrimination panel;
- compares the incidence of the study and of the population **by
  attained age**, with the study’s confidence intervals, on a log scale
  (the official panel compared a loess-smoothed population curve with
  the raw study proportions);
- offers an optional panel of the E/O ratio of each risk-score group
  with its confidence interval (`panels = "eo"`), and a choice of
  panels;
- weights the risk-score densities of a sampled study by the inverse
  inclusion probabilities;
- has a title, a subtitle describing the study, and a caption explaining
  the points, bars and tests.

The arguments of the official function keep their meaning, except
`risk.score.plot.percent.smooth`, which is accepted and ignored.
`panels`, `base.size` and `draw` are new.

## What stays the same

- The functions, their arguments, and the names and meaning of their
  results, except where listed above.
- The model: yearly hazards, the calibration of the baseline hazard, the
  simulation of SNP genotypes (with family history), and the donors for
  missing risk factors: the reference rows in the profile’s percentile
  bin (100 unweighted bins of the observed linear predictor; the
  averages over donors are weighted).
- Results that were not affected by a defect agree with the official
  version to rounding error: the test suite compares risks, linear
  predictors and validation statistics with recorded outputs of the
  official version 1.40.0 and of 1.33.1.9002.
- When
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
  computes the risks itself, each subject’s risk covers her follow-up
  within the interval, as before; supply `predicted.risk` to validate
  fixed-horizon risks.
- The risk score of the AUC and of the risk groups is the linear
  predictor, which excludes age.
