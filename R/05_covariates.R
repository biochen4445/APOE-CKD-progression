# Purpose: Derive index covariates (SAP Sec. 5) and assemble the analysis data set.
# Input(s): data/derived/{inputs,cohort,outcomes,apoe}.rds
# Output(s): data/derived/analysis.rds; output/tables/s05_missingness.csv
# Key methods/packages: data.table
# Notes: Missing data handled by complete-case analysis (SAP silent; review A4). This table
#        reports the share missing for each model covariate.

source("R/00_setup.R")
log_msg("== 05 covariates and analysis data ==")
d <- load_derived("inputs", cfg)
cohort <- load_derived("cohort", cfg)
out <- load_derived("outcomes", cfg)
apoe <- load_derived("apoe", cfg)

covs <- derive_covariates(cfg, d, cohort, cl)
analysis <- assemble_analysis(cfg, out, covs, apoe)
save_derived(analysis, "analysis", cfg)

mvars <- model_covs(cfg, "m3")
miss <- data.table::data.table(
  variable = mvars,
  n_missing = vapply(mvars, function(v) sum(is.na(analysis[[v]])), numeric(1)))
miss[, pct_missing := round(100 * n_missing / nrow(analysis), 2)]
miss <- rbind(miss, data.table::data.table(
  variable = "Model 3 complete cases",
  n_missing = nrow(analysis) - nrow(complete_sample(analysis, mvars)),
  pct_missing = round(100 * (1 - nrow(complete_sample(analysis, mvars)) / nrow(analysis)), 2)))
write_table(miss, "s05_missingness", cfg)
log_msg("05 done: analysis n = ", nrow(analysis), "; Model 3 complete cases = ",
        nrow(complete_sample(analysis, mvars)))
