#' Example data: a breast cancer risk model and validation studies
#'
#' @description
#' Example inputs for [computeAbsoluteRisk()],
#' [computeAbsoluteRiskSplitInterval()], [ModelValidation()] and
#' [plotModelValidation()]: a breast cancer relative-risk model with 13 risk
#' factors and 72 SNPs, a reference dataset, incidence and mortality rates,
#' profiles to apply the model to, and simulated studies to validate it. Load
#' them all with `data(bc_data, package = "iCARE")`.
#'
#' @format
#' **The model**
#' * `bc_model_formula`: formula of the relative-risk model, with the 13 risk
#'   factors below and an interaction of HRT use with body mass index.
#' * `bc_model_cov_info`: list describing each risk factor: its name, its type
#'   (`"continuous"` or `"factor"`), and for factors the levels and the
#'   reference level.
#' * `bc_model_log_or`: named vector of the 77 log odds ratios of the model,
#'   one per column of its design matrix.
#' * `bc_model_log_or_post_50`: the log odds ratios for women aged 50 years
#'   or older (an example for [computeAbsoluteRiskSplitInterval()]).
#' * `bc_72_snps`: data frame of 72 breast cancer SNPs with columns
#'   `snp.name`, `snp.odds.ratio` (per-allele odds ratio) and `snp.freq`
#'   (risk-allele frequency), from Michailidou et al. (2017).
#'
#' **The population**
#' * `ref_cov_dat`: data frame of the 13 risk factors in 14,137 women,
#'   representative of the US population (imputed from NHIS 2010 and NHANES;
#'   references 4 and 5).
#' * `ref_cov_dat_post_50`: the reference dataset for women aged 50 years or
#'   older.
#' * `bc_inc`: matrix of the annual breast cancer incidence rates by single
#'   year of age (columns: age, rate) from SEER (reference 3).
#' * `mort_inc`: matrix of the annual mortality rates from causes other than
#'   breast cancer by single year of age (columns: age, rate) from CDC WONDER
#'   (reference 1).
#'
#' **Profiles**
#' * `new_cov_prof`: data frame of the risk factors of three women.
#' * `new_snp_prof`: data frame of the genotypes (risk-allele counts) of the
#'   72 SNPs of the same three women; `NA` marks a missing genotype.
#'
#' **Validation studies**
#' * `validation.cohort.data`: simulated full cohort of 50,000 women with the
#'   columns
#'   * `id`: subject identifier;
#'   * the 13 risk factors: `famhist` (family history of breast cancer in a
#'     first-degree relative), `menarche_dec`, `parity`, `birth_dec`,
#'     `agemeno_dec`, `height_dec`, `bmi_dec` (categories of age at menarche,
#'     number of births, age at first birth, age at menopause, height and body
#'     mass index), `rd_menohrt` (menopausal status and HRT use),
#'     `rd2_everhrt_e`, `rd2_everhrt_c`, `rd2_currhrt` (ever use of
#'     estrogen-only or combined HRT, current HRT use), `alcoholdweek_dec`
#'     (categories of drinks per week) and `ever_smoke`;
#'   * `study.entry.age`, `study.exit.age`: ages at entry and exit;
#'   * `observed.outcome`: 1 if breast cancer developed during follow-up;
#'   * `time.of.onset`: years from entry to onset (`Inf` without onset);
#'   * `observed.followup`: years of follow-up;
#'   * `inclusion`: 0 (a placeholder for the indicator of inclusion in the
#'     nested case-control study; see the examples).
#' * `validation.nested.case.control.data`: a simulated case-control study of
#'   5,285 women nested in the cohort, with the same columns and the 72 SNP
#'   genotypes (columns named by the SNPs' rs identifiers).
#' * `output`: a result of [computeAbsoluteRisk()] from an earlier version of
#'   iCARE.
#'
#' @references
#' 1. Centers for Disease Control and Prevention (CDC), National Center for
#'    Health Statistics (NCHS). Underlying Cause of Death 1999-2011 on CDC
#'    WONDER Online Database, released 2014.
#' 2. Michailidou K, Beesley J, Lindstrom S, et al. (2015). Genome-wide
#'    association analysis of more than 120,000 individuals identifies 15 new
#'    susceptibility loci for breast cancer. *Nature Genetics* 47: 373-380.
#' 3. Surveillance, Epidemiology, and End Results (SEER) Program SEER*Stat
#'    Database: Incidence - SEER 18 Regs Research Data, Nov 2011 Sub.
#'    National Cancer Institute.
#' 4. 2010 National Health Interview Survey (NHIS) Public Use Data Release.
#'    National Center for Health Statistics, 2011.
#' 5. National Health and Nutrition Examination Survey Questionnaire.
#'    National Center for Health Statistics, 2010.
#' 6. Michailidou K, Lindstrom S, Dennis J, et al. (2017). Association
#'    analysis identifies 65 new breast cancer risk loci. *Nature* 551: 92-94.
#'
#' @examples
#' data(bc_data, package = "iCARE")
#' bc_model_formula
#' head(bc_inc)
#' head(bc_72_snps)
#'
#' # Inclusion probabilities of the nested case-control study, estimated from
#' # the cohort by logistic regression (as in the vignette)
#' cohort <- validation.cohort.data
#' ncc <- validation.nested.case.control.data
#' cohort$inclusion <- as.numeric(cohort$id %in% ncc$id)
#' selection <- glm(
#'   inclusion ~ observed.outcome * (study.entry.age + observed.followup),
#'   data = cohort, family = binomial
#' )
#' probability <- fitted(selection)[match(ncc$id, cohort$id)]
#' summary(probability)
#' @name bc_data
#' @aliases bc_72_snps bc_inc bc_model_formula bc_model_cov_info bc_model_log_or bc_model_log_or_post_50 mort_inc new_cov_prof new_snp_prof ref_cov_dat ref_cov_dat_post_50 validation.cohort.data validation.nested.case.control.data output
#' @docType data
#' @keywords datasets
NULL
