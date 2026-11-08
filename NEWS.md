# iCARE 2.0.0 (unofficial distribution)

The first release of this unofficial, refactored distribution of iCARE. It is
rewritten from scratch on the basis of Bioconductor iCARE 1.40.0 and of
1.33.1.9002, a patched development version prepared by the maintainer of this
distribution. The functions, their arguments and the names of their results
are unchanged. `vignette("differences")` explains each change below.

## Changes that can alter results

* The reference weights (`model.ref.dataset.weights`) are used correctly in
  the calibration of the baseline hazard; the official version scrambled
  non-uniform weights across subjects and ages.
* Fractional start ages and interval lengths are used exactly, instead of
  rounding a fractional start up and counting a fraction of a year as a whole
  year.
* Rates are needed only for the years of age inside the intervals. In the
  official version, a missing rate at the age where an interval ends made the
  risks `NA` or, when intervals differed, silently replaced the risk of every
  profile by the average over its percentile bin; this also affected the
  risks that `ModelValidation()` computed for the example studies.
* For a profile with missing risk factors, the risk combines the profile's
  own observed risk factors with the donors' values of the missing ones, as
  in the published method (Pal Choudhury et al. 2020, equation 3), instead of
  averaging the donors' own risks. Risks of profiles in the tails of the
  risk distribution can change by more than 10%.
* Each profile with missing risk factors gets the risk over its own interval,
  independently of the other profiles.
* Rates given by age band are read as yearly rates for every age of the band.
* For SNP-only models, `refs.risk` holds the risks of simulated people
  instead of averages over the `n.imp` sets of simulated people.
* `computeAbsoluteRiskSplitInterval()` combines the two intervals as
  `AR1 + S1 * AR2`, accounting for death from other causes before the cut
  age, and keeps the reference weights after the cut age.
* `ModelValidation()`:
  * the AUC counts tied risk scores as 1/2;
  * the Hosmer-Lemeshow statistic uses the variance under the null
    hypothesis by default (`hosmer.lemeshow.variance = "observed"` restores
    the earlier statistic);
  * in sampled studies, the E/O confidence interval accounts for the
    variability of E;
  * risk-score groups handle tied scores and empty groups, and the degrees of
    freedom of the tests follow the number of groups formed;
  * the study incidence is the person-time rate by attained age, with
    confidence intervals;
  * a supplied `predicted.risk` or `linear.predictor` is never recomputed,
    the reference weights and `n.imp` of the model object are used, and
    onsets after exit are censored consistently;
  * cases diagnosed at exit are kept when ages have decimals (rounding error
    in exit minus entry turned many of them into non-cases);
  * the reference population's risks use each row's own window.

## New features

* `Category_Results` has numeric columns and adds each group's E/O ratio with
  its confidence interval and the numbers of subjects and cases;
  `Category_Specific_Expected_to_Observed_Ratio` is new.
* `plotModelValidation()` is redrawn with ggplot2 and patchwork: square
  panels with their statistics, incidence by attained age, an optional panel
  of E/O by risk group (`panels`), `base.size`, and the figure is returned
  for saving with `ggplot2::ggsave()`.
* `print()`, `summary()` and `plot()` methods for the results of
  `computeAbsoluteRisk()`, `computeAbsoluteRiskSplitInterval()` and
  `ModelValidation()`.
* The calibrated baseline hazard is returned as
  `attr(result, "calibration")`.
* With SNPs, family history may be coded in any way, for example as a factor
  with reference level "1" (the official version stopped with an error);
  the attenuation of its log relative risk applies to the effect of family
  history.
* Columns of `apply.snp.profile` named after the SNPs are matched by name.
* Clear error messages for invalid inputs, including formulas that give
  infinite or undefined design values.

## Implementation

* Pure, vectorized R: the C code and the `use.c.code` option are removed (the
  option is accepted with a deprecation warning and ignored). Every
  benchmarked workload is faster than the official version with its compiled
  code; see `vignette("performance")`.
* `plotrix`, `gtools` and `Hmisc` are no longer dependencies; `ggplot2` and
  `patchwork` are.
* The code is organized in small documented modules (ages and rates, the
  covariate and SNP models, the risk engine, missing risk factors, and the
  validation statistics), with a test suite that compares the results with
  brute-force computations and with recorded outputs of earlier versions.
* A simulation study (`tools/simulations/validation-simulation.R`) checks
  the coverage of the confidence intervals and the size of the tests of
  `ModelValidation()`. The relative-risk test, unchanged from earlier
  versions, rejects calibrated models too often when risk groups have few
  cases, and its documentation now says so.
