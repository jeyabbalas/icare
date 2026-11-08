# iCARE: Individualized Coherent Absolute Risk Estimation

iCARE builds, applies and validates models for the absolute risk of a
disease. An absolute risk model combines

1.  a model for the **relative risks** of the risk factors, typically
    estimated in a cohort or case-control study,

2.  a **reference dataset** of risk factors representing the population
    in which the model is applied,

3.  the population's age-specific **disease incidence rates**, and

4.  optionally, age-specific rates of **competing mortality**.

The main functions are

- [`computeAbsoluteRisk()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRisk.md):
  estimate absolute risks for covariate and/or SNP profiles, including
  profiles with missing risk factors;

- [`computeAbsoluteRiskSplitInterval()`](https://jeyabbalas.github.io/icare/reference/computeAbsoluteRiskSplitInterval.md):
  the same when the relative-risk model or the risk-factor distribution
  changes at a cut age;

- [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md):
  validate a model in a cohort or nested case-control study (calibration
  and discrimination);

- [`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md):
  publication-ready figures of a validation.

## Unofficial distribution

This is an **unofficial**, refactored distribution of iCARE. The
official package is distributed by Bioconductor
(<https://bioconductor.org/packages/iCARE/>) and on GitHub
(<https://github.com/parichoy/iCARE>). This version keeps the
user-facing interface, fixes known defects and has no compiled code; see
[`vignette("differences", package = "iCARE")`](https://jeyabbalas.github.io/icare/articles/differences.md)
for every difference.

## References

Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R package to
build, validate and apply absolute risk models. *PLoS ONE* 15(2):
e0228198.
[doi:10.1371/journal.pone.0228198](https://doi.org/10.1371/journal.pone.0228198)

## See also

Useful links:

- <https://jeyabbalas.github.io/icare/>

- <https://github.com/jeyabbalas/icare>

- Report bugs at <https://github.com/jeyabbalas/icare/issues>

## Author

**Maintainer**: Jeya Balaji Balasubramanian <jeyabbalas@gmail.com>
(Refactored, unofficial distribution)

Authors:

- Jeya Balaji Balasubramanian <jeyabbalas@gmail.com> (Refactored,
  unofficial distribution)

- Parichoy Pal Choudhury (Author of the original iCARE package)

- Paige Maas (Author of the original iCARE package)

- William Wheeler (Author of the original iCARE package)

- Nilanjan Chatterjee (Author of the original iCARE package)
