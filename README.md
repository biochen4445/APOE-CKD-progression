# APOE and CKD progression in an EMR-linked biobank

Analysis code for an EMR adaptation of Hsu et al., *JAMA* 2005;293:2892-2899 ("Apolipoprotein E and progression of chronic kidney disease"), applied to the China Medical University Hospital (CMUH) iHi cohort.

The code follows the statistical analysis plan **SAP v0.1-HL (2026-09-30)**. Every analytic choice the SAP fixes, or leaves open, is a setting in `config/config.yml`; the scripts contain no hard-coded thresholds. Points in the SAP that still need a decision are listed in `docs/SAP_review_v0.1_HL.md`.

> No patient data are stored in this repository. `data/raw/`, `data/derived/` and `output/` are git-ignored. The pipeline runs end to end on synthetic data from `sim/simulate_data.R`.

## Design in brief

| Item | Definition (SAP section) |
| --- | --- |
| Index date | First outpatient serum creatinine at age 45-65 years (Sec. 2) |
| Exposure | APOE from rs429358 and rs7412; e2 and e4 allele counts, reference e3 (Sec. 3) |
| Primary outcome | First of: creatinine rise >= 0.4 mg/dL confirmed >= 90 days later; CKD hospitalisation (principal diagnosis not AKI); death with CKD (Sec. 4) |
| Censoring | Death, last CMUH encounter, or EMR end date (Sec. 5) |
| Models | Cox; Model 1 age + sex; Model 2 + BMI, diabetes, SBP, DBP, antihypertensive use, CHD, index eGFR; Model 3 + HDL, LDL, log TG (Sec. 5-6) |
| Primary test | 2-df likelihood-ratio test of e2 and e4 counts in Model 3 (proposed; review item A1) |

## Repository layout

```
config/
  config.yml               all analytic parameters
  codelists/codelists.csv  ICD-9-CM, ICD-10-CM, procedure and ATC code lists ("status" = SAP or draft-verify)
R/
  00_setup.R               packages, config, helper functions (sourced by every stage)
  01_import_validate.R     schema checks, date and code standardisation
  02_apoe_genotype_qc.R    genotype calls, call rate, HWE, allele frequencies
  03_cohort.R              index date, inclusion/exclusion, Figure 0 flow
  04_outcomes.R            primary and secondary outcomes, censoring
  05_covariates.R          index covariates, analysis data set, missingness
  06_table1.R              Table 1 by genotype and by outcome status, with SMD
  07_incidence.R           Figure 1: incidence per 1,000 person-years by genotype
  08_cox_primary.R         Table 2 (Model 1-3), genotype models, PH tests, secondary outcomes
  09_bias_checks.R         follow-up and death by genotype, Fine-Gray, testing intensity
  10_sensitivity.R         SAP Sec. 6.5 sensitivity analyses plus review-proposed ones
  11_power.R               Schoenfeld events required and detectable HR
  functions/               reusable functions (eGFR, APOE calling, cohort, outcomes, covariates, models)
sim/simulate_data.R        synthetic data matching the input schema
tests/testthat/            unit tests (eGFR, APOE calls, HWE, confirmed-event logic, code matching)
docs/                      SAP review notes
data/README.md             input schema for the EMR extracts
run_all.R                  runs stages 01-11 and writes output/logs/run_log.md
```

## Running

Requirements: R >= 4.2 with `data.table`, `survival`, `ggplot2`, `yaml` (and `testthat` for tests).

```r
install.packages(c("data.table", "survival", "ggplot2", "yaml", "testthat", "renv"))
```

From the repository root:

```bash
Rscript run_all.R --simulate                 # synthetic data -> full pipeline (about 1 min)
Rscript -e 'testthat::test_dir("tests/testthat")'
```

For the real analysis, place the extracts described in `data/README.md` in `data/raw/`, set `data_source: "emr"` and `emr_end_date` in `config/config.yml`, then run `Rscript run_all.R`. Each stage can also be run on its own (for example `Rscript R/08_cox_primary.R`) once the earlier stages have produced their files in `data/derived/`.

Outputs go to `output/tables/` (CSV), `output/figures/` (PDF and PNG) and `output/logs/`. `run_log.md` records the run date, git commit, md5 of the config and code lists, package versions and the key estimates, so any reported number can be traced to a run.

## Before the first real run

1. Resolve the open items in `docs/SAP_review_v0.1_HL.md` and the SAP checklist (EMR date range, index setting, death data source).
2. Verify every code list row marked `draft-verify`, and add NHI dialysis order codes (not included).
3. Confirm the counted allele and strand for rs429358 (C = e4) and rs7412 (T = e2) against the genotyping array manifest.
4. Run `renv::init()` and commit `renv.lock` so the package versions are fixed.

## Publishing to GitHub

```bash
git init
git add .
git commit -m "Analysis code for APOE x CKD progression (SAP v0.1-HL)"
git remote add origin https://github.com/<account>/<repo>.git
git push -u origin main
```

Check `git status` before each commit: only code, config, code lists, documentation and tests should be tracked. Add a licence file of your choice before making the repository public.

## Reference

Hsu CC, Kao WHL, Coresh J, et al. Apolipoprotein E and progression of chronic kidney disease. *JAMA*. 2005;293(23):2892-2899. doi:10.1001/jama.293.23.2892
