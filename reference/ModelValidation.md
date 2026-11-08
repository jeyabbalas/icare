# Validate an absolute risk model in a cohort or a sample from a cohort

Assesses the calibration and discrimination of an absolute risk model in
independent data from a cohort study (full cohort) or from a sample
drawn from a cohort with known inclusion probabilities, such as a nested
case-control study. Risks are validated over a fixed interval after
study entry (e.g. 5-year risks) or over each subject's observed
follow-up.

The function reports

- **calibration in the large**: the ratio of the expected to the
  observed number of cases (E/O) with its confidence interval;

- **calibration by risk group**: subjects are grouped by their risk
  score (deciles by default), and the observed and expected absolute
  risks, relative risks and E/O ratios of the groups are compared, with
  a Hosmer-Lemeshow test of the absolute risks and a test of the
  relative risks;

- **discrimination**: the area under the ROC curve (AUC) of the risk
  score;

- the study's **incidence rates by attained age**, to compare with the
  population incidence rates of the model.

## Usage

``` r
ModelValidation(
  study.data,
  total.followup.validation = FALSE,
  predicted.risk = NULL,
  predicted.risk.interval = NULL,
  linear.predictor = NULL,
  iCARE.model.object = list(model.formula = NULL, model.cov.info = NULL, model.snp.info =
    NULL, model.log.RR = NULL, model.ref.dataset = NULL, model.ref.dataset.weights =
    NULL, model.disease.incidence.rates = NULL, model.competing.incidence.rates = NULL,
    model.bin.fh.name = NA, apply.cov.profile = NULL, apply.snp.profile = NULL, n.imp =
    5),
  number.of.percentiles = 10,
  reference.entry.age = NULL,
  reference.exit.age = NULL,
  predicted.risk.ref = NULL,
  linear.predictor.ref = NULL,
  linear.predictor.cutoffs = NULL,
  dataset = "Example Dataset",
  model.name = "Example Risk Prediction Model",
  hosmer.lemeshow.variance = c("expected", "observed")
)
```

## Arguments

- study.data:

  Data frame with one row per subject and the columns

  - `observed.outcome`: 1 if the disease developed during follow-up, 0
    otherwise;

  - `study.entry.age`: age at entry into the study;

  - `study.exit.age`: age at the last follow-up;

  - `time.of.onset`: years from entry to disease onset (`Inf` or `NA`
    for subjects without the disease);

  - `sampling.weights` (optional): for a sample from the cohort, such as
    a nested case-control study, each subject's probability of inclusion
    in the sample.

- total.followup.validation:

  If `TRUE`, validate risks over each subject's observed follow-up; if
  `FALSE` (default), over `predicted.risk.interval`.

- predicted.risk:

  Optional vector of predicted absolute risks, one per subject. If
  `NULL`, the risks are computed from `iCARE.model.object`.

- predicted.risk.interval:

  Years after entry over which risks are validated (e.g. 5): one value,
  or one per subject. Not used when `total.followup.validation = TRUE`.

- linear.predictor:

  Optional vector of risk scores (the linear predictor \\x'\beta\\), one
  per subject. If `NULL`, computed from `iCARE.model.object`. A supplied
  `predicted.risk` or `linear.predictor` is always used as given; only a
  missing one is computed.

- iCARE.model.object:

  Named list of arguments of
  [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)
  describing the model and the subjects: `model.formula`,
  `model.cov.info`, `model.snp.info`, `model.log.RR`,
  `model.ref.dataset`, `model.ref.dataset.weights`,
  `model.disease.incidence.rates`, `model.competing.incidence.rates`,
  `model.bin.fh.name`, `apply.cov.profile` and `apply.snp.profile` (the
  subjects' risk factors, one row per row of `study.data`) and `n.imp`.
  Settings that earlier versions needed here (`use.c.code`, `return.lp`,
  `return.refs.risk`) are accepted and ignored.

- number.of.percentiles:

  Number of risk-score groups formed at quantiles of the risk score
  (default 10, deciles).

- reference.entry.age, reference.exit.age:

  Optional ages (one value, or one per row of `model.ref.dataset`) over
  which to compute the absolute risks of the reference population,
  returned as `Reference_Absolute_Risk`.

- predicted.risk.ref, linear.predictor.ref:

  Optional absolute risks and risk scores of the reference population,
  returned as given.

- linear.predictor.cutoffs:

  Optional cut points of the risk score defining the risk groups,
  instead of quantiles.

- dataset, model.name:

  Names of the dataset and the model, shown by
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md).

- hosmer.lemeshow.variance:

  Variance in the Hosmer-Lemeshow statistic: `"expected"` (default) or
  `"observed"` (earlier versions); see Details.

## Value

An object of class `"icareValid"`: a list with

- `Subject_Specific_Observed_Outcome`, `Adjusted_Followup`: outcome and
  follow-up over the risk prediction interval;

- `Risk_Prediction_Interval`: the interval, e.g. `"5 years"` or
  `"Observed Followup"`;

- `Subject_Specific_Predicted_Absolute_Risk`,
  `Subject_Specific_Risk_Score`: each subject's predicted risk and risk
  score;

- `Reference_Absolute_Risk`, `Reference_Risk_Score`: those of the
  reference population (`NULL` unless supplied or requested);

- `Population_Incidence_Rate`: the model's disease incidence rates;

- `Study_Incidence_Rate`: data frame of the study's incidence rates by
  attained age (`ages`, `study_incidence`), their confidence limits
  (`lower`, `upper`), and the (weighted) `cases` and `person_years`;

- `Category_Results`: data frame with one row per risk-score group: the
  group (`Categories`); its observed and predicted absolute risks and
  relative risks with confidence intervals for the observed ones; its
  E/O ratio with confidence interval (`Expected_to_Observed_Ratio`,
  `CI_Expected_to_Observed_Ratio_Lower`,
  `CI_Expected_to_Observed_Ratio_Upper`); and its numbers of subjects
  and cases in the data;

- `Category_Specific_Observed_Absolute_Risk`,
  `Category_Specific_Predicted_Absolute_Risk`,
  `Category_Specific_Observed_Relative_Risk`,
  `Category_Specific_Predicted_Relative_Risk`,
  `Category_Specific_Expected_to_Observed_Ratio`: the same, as vectors;

- `Variance_Matrix_Absolute_Risk`, `Variance_Matrix_LogRelative_Risk`:
  covariance matrices of the groups' observed absolute risks and log
  relative risks;

- `Hosmer_Lemeshow_Results`, `HL_pvalue`: the Hosmer-Lemeshow test (an
  `htest` object) and its p-value;

- `RR_test_result`, `RR_test_pvalue`: the relative-risk test and its
  p-value;

- `AUC`, `Variance_AUC`, `CI_AUC`: the AUC, its variance and 95%
  confidence interval;

- `Overall_Expected_to_Observed_Ratio`,
  `CI_Overall_Expected_to_Observed_Ratio`: E/O and its 95% confidence
  interval;

- `input.args`: the model formula, the study data, the names and the
  settings, used by [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

**Follow-up.** Each subject is followed from `study.entry.age` for the
shorter of the risk prediction interval and the observed follow-up, and
counts as a case only if the disease develops within that time
(`Adjusted_Followup`, `Subject_Specific_Observed_Outcome`). When the
risks are computed here, each subject's risk covers that same time from
`study.entry.age`, so that a subject followed for 3 years of a 5-year
interval contributes a 3-year risk. To validate fixed-horizon risks,
compute them yourself and supply `predicted.risk` and
`linear.predictor`. Cases whose onset is recorded after exit are treated
as non-cases censored at exit, with a warning.

**Sampled data.** When `study.data` has a `sampling.weights` column
(each subject's probability of inclusion in the sample), all estimates
are weighted by the inverse probabilities, and the variances add the
variability due to sampling.

**Risk-score groups.** Subjects are grouped at quantiles of
`linear.predictor` (weighted quantiles in sampled data), or at
`linear.predictor.cutoffs`. Values tied at a quantile form their own
group, and groups left empty are dropped; the tests use the number of
groups actually formed, \\G\\.

**Tests.** The Hosmer-Lemeshow statistic is \\\sum_g (O_g - E_g)^2 /
V_g\\ on \\G\\ degrees of freedom, where \\O_g\\ and \\E_g\\ are the
observed and expected absolute risks of group \\g\\. With
`hosmer.lemeshow.variance = "expected"` (the default), \\V_g\\ is the
variance under the null hypothesis of calibration, \\E_g(1 - E_g)/N_g\\
(plus a sampling term in sampled data), as in the classical
Hosmer-Lemeshow test. With `"observed"`, \\V_g\\ is the variance of
\\O_g\\, as in earlier versions of iCARE; that test rejects too often
when groups have few cases and is infinite when a group has none. The
relative-risk test compares the observed and expected relative risks of
the groups (each group's risk divided by the mean over groups) by a Wald
test of the log relative risks on \\G - 1\\ degrees of freedom; it needs
a case in every group. Like the `"observed"` Hosmer-Lemeshow test, it
rejects a calibrated model too often when groups have few cases. In
simulations at the 5% level, it rejected 11.5% of the time with 10 cases
per group on average, 7.7% with 30 and 6.7% with 50 (see
[`vignette("model-validation")`](https://jeyabbalas.github.io/icare/articles/model-validation.md)).
With few cases, form fewer groups.

**Confidence intervals** are 95% Wald intervals: on the log scale for
E/O ratios and relative risks, and on the risk scale for absolute risks
and the AUC. The E/O interval in a full cohort treats the expected risks
as known; in sampled data it also accounts for the sampling variability
of the expected risk and its covariance with the observed risk.

**AUC.** The probability that a case has a higher risk score than a
non-case, counting ties as 1/2, weighted by the inverse inclusion
probabilities in sampled data, with the variance of DeLong et al.
(1988). The risk score is the linear predictor, which excludes age, so
the AUC measures discrimination by the risk factors within the age
distribution of the study.

**Study incidence.** Onsets divided by person-years at each year of
attained age over the observed follow-up, with exact Poisson intervals
in a full cohort and weighted estimates with Wald intervals in sampled
data.

## References

Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R package to
build, validate and apply absolute risk models. *PLOS ONE* 15(2):
e0228198.
[doi:10.1371/journal.pone.0228198](https://doi.org/10.1371/journal.pone.0228198)

DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas
under two or more correlated receiver operating characteristic curves: a
nonparametric approach. *Biometrics* 44(3): 837-845.

Hosmer DW, Lemeshow S (1980). Goodness of fit tests for the multiple
logistic regression model. *Communications in Statistics - Theory and
Methods* 9(10): 1043-1069.

## See also

[`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md),
[`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md)

## Examples

``` r
data(bc_data, package = "iCARE")

# 5-year risks in a random sample of the full validation cohort
set.seed(50)
cohort <- validation.cohort.data[sample(nrow(validation.cohort.data), 5000), ]
covariates <- all.vars(bc_model_formula)[-1]
validation <- ModelValidation(
  study.data = cohort,
  predicted.risk.interval = 5,
  iCARE.model.object = list(
    model.formula = bc_model_formula,
    model.cov.info = bc_model_cov_info,
    model.log.RR = bc_model_log_or,
    model.ref.dataset = ref_cov_dat,
    model.disease.incidence.rates = bc_inc,
    model.competing.incidence.rates = mort_inc,
    apply.cov.profile = cohort[, covariates]
  ),
  dataset = "Simulated cohort",
  model.name = "Breast cancer risk model"
)
validation
#> iCARE model validation
#> 
#> Dataset:                   Simulated cohort
#> Model:                     Breast cancer risk model
#> Model formula:             diagnosis ~ famhist + as.factor(menarche_dec) +
#>                            as.factor(parity) + as.factor(birth_dec) +
#>                            as.factor(agemeno_dec) + as.factor(height_dec) +
#>                            as.factor(bmi_dec) + as.factor(rd_menohrt) +
#>                            rd2_everhrt_e + rd2_everhrt_c + rd2_currhrt +
#>                            as.factor(alcoholdweek_dec) + as.factor(ever_smoke) +
#>                            as.factor(rd_menohrt) * as.factor(bmi_dec)
#> Risk prediction interval:  5 years
#> Subjects (cases):          5,000 (96)
#> Follow-up (years):         mean 5.00 (range 5-5)
#> Age at entry (years):      mean 60.9 (range 50-72)
#> 
#> Calibration
#>   E/O:                     0.954 (95% CI 0.783-1.163)
#>   Hosmer-Lemeshow test:    chi-square 16.41, df 10, p = 0.09
#>                            (variance from the expected risks)
#>   Relative-risk test:      chi-square 21.55, df 9, p = 0.01
#> 
#> Discrimination
#>   AUC:                     0.579 (95% CI 0.521-0.638)
#> 
#> Risk-score groups
#>             Group Subjects Cases Exp. %      Obs. % (95% CI)      E/O (95% CI)
#>   [-0.961,-0.505]      500    11   1.06  2.20 (0.91 to 3.49)  0.48 (0.27-0.86)
#>   (-0.505,-0.385]      500     2   1.25 0.40 (-0.15 to 0.95) 3.14 (0.79-12.50)
#>   (-0.385,-0.301]      500     4   1.36  0.80 (0.02 to 1.58)  1.70 (0.64-4.51)
#>    (-0.301,-0.22]      500     7   1.51  1.40 (0.37 to 2.43)  1.08 (0.52-2.25)
#>    (-0.22,-0.137]      500    14   1.63  2.80 (1.35 to 4.25)  0.58 (0.35-0.98)
#>  (-0.137,-0.0512]      500     7   1.79  1.40 (0.37 to 2.43)  1.28 (0.61-2.67)
#>  (-0.0512,0.0386]      500    13   1.94  2.60 (1.21 to 3.99)  0.75 (0.44-1.27)
#>    (0.0386,0.146]      500    11   2.12  2.20 (0.91 to 3.49)  0.96 (0.54-1.72)
#>     (0.146,0.304]      500    13   2.50  2.60 (1.21 to 3.99)  0.96 (0.56-1.64)
#>      (0.304,1.07]      500    14   3.16  2.80 (1.35 to 4.25)  1.13 (0.67-1.89)
validation$Category_Results
#>          Categories Observed_Absolute_Risk Predicted_Absolute_Risk
#> 1   [-0.961,-0.505]                  0.022              0.01060187
#> 2   (-0.505,-0.385]                  0.004              0.01254380
#> 3   (-0.385,-0.301]                  0.008              0.01357997
#> 4    (-0.301,-0.22]                  0.014              0.01507666
#> 5    (-0.22,-0.137]                  0.028              0.01633529
#> 6  (-0.137,-0.0512]                  0.014              0.01793208
#> 7  (-0.0512,0.0386]                  0.026              0.01938218
#> 8    (0.0386,0.146]                  0.022              0.02115120
#> 9     (0.146,0.304]                  0.026              0.02497459
#> 10     (0.304,1.07]                  0.028              0.03158403
#>    CI_Absolute_Risk_Lower CI_Absolute_Risk_Upper Observed_Relative_Risk
#> 1             0.009142639            0.034857361              1.1458333
#> 2            -0.001532619            0.009532619              0.2083333
#> 3             0.000191423            0.015808577              0.4166667
#> 4             0.003701510            0.024298490              0.7291667
#> 5             0.013539503            0.042460497              1.4583333
#> 6             0.003701510            0.024298490              0.7291667
#> 7             0.012051188            0.039948812              1.3541667
#> 8             0.009142639            0.034857361              1.1458333
#> 9             0.012051188            0.039948812              1.3541667
#> 10            0.013539503            0.042460497              1.4583333
#>    Predicted_Relative_Risk CI_Relative_Risk_Lower CI_Relative_Risk_Upper
#> 1                0.5788258             0.66115533              1.9858178
#> 2                0.6848483             0.05302386              0.8185519
#> 3                0.7414199             0.16030426              1.0830100
#> 4                0.8231339             0.35919611              1.4802054
#> 5                0.8918507             0.90462652              2.3509549
#> 6                0.9790305             0.35919611              1.4802054
#> 7                1.0582005             0.82218267              2.2303649
#> 8                1.1547832             0.66115533              1.9858178
#> 9                1.3635270             0.82218267              2.2303649
#> 10               1.7243802             0.90462652              2.3509549
#>    Expected_to_Observed_Ratio CI_Expected_to_Observed_Ratio_Lower
#> 1                   0.4819032                           0.2686254
#> 2                   3.1359488                           0.7864526
#> 3                   1.6974963                           0.6395965
#> 4                   1.0769040                           0.5160683
#> 5                   0.5834031                           0.3480799
#> 6                   1.2808632                           0.6138085
#> 7                   0.7454683                           0.4359466
#> 8                   0.9614182                           0.5359195
#> 9                   0.9605610                           0.5617319
#> 10                  1.1280012                           0.6730073
#>    CI_Expected_to_Observed_Ratio_Upper Number_of_Subjects Number_of_Cases
#> 1                            0.8645151                500              11
#> 2                           12.5044730                500               2
#> 3                            4.5051739                500               4
#> 4                            2.2472263                500               7
#> 5                            0.9778190                500              14
#> 6                            2.6728375                500               7
#> 7                            1.2747501                500              13
#> 8                            1.7247459                500              11
#> 9                            1.6425583                500              13
#> 10                           1.8905986                500              14
```
