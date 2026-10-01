# Purpose: Cox model helpers -- formulas, complete-case samples, HR extraction, joint tests.
# Key packages: survival

model_covs <- function(cfg, set, drop = NULL, add = NULL) {
  covs <- unlist(cfg$models[[set]])
  if (isTRUE(cfg$models$adjust_pcs)) covs <- c(covs, paste0("pc", seq_len(cfg$models$n_pcs)))
  setdiff(unique(c(covs, add)), drop)
}

surv_lhs <- function(outcome, timescale = c("followup", "age"), status = FALSE) {
  timescale <- match.arg(timescale)
  ev <- if (status) sprintf("%s_status", outcome) else sprintf("%s_event", outcome)
  if (timescale == "followup") sprintf("Surv(tstart, %s_tstop, %s)", outcome, ev)
  else sprintf("Surv(age_entry, %s_age_exit, %s)", outcome, ev)
}

surv_vars <- function(outcome, timescale = "followup") {
  if (timescale == "followup") c("tstart", paste0(outcome, "_tstop"), paste0(outcome, "_event"))
  else c("age_entry", paste0(outcome, "_age_exit"), paste0(outcome, "_event"))
}

cox_formula <- function(outcome, terms, timescale = "followup", cluster = NULL) {
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  if (!is.null(cluster)) rhs <- paste0(rhs, " + cluster(", cluster, ")")
  stats::as.formula(paste(surv_lhs(outcome, timescale), "~", rhs))
}

complete_sample <- function(data, vars) data[stats::complete.cases(data[, ..vars])]

tidy_hr <- function(fit, terms) {
  b  <- stats::coef(fit)[terms]
  se <- sqrt(diag(stats::vcov(fit)))[terms]   # robust variance when cluster() is used
  data.table::data.table(term = terms, hr = exp(b), lcl = exp(b - 1.959964 * se),
                         ucl = exp(b + 1.959964 * se), p = 2 * stats::pnorm(-abs(b / se)))
}

# Joint test of `terms`: likelihood-ratio when variance is model-based, Wald when robust
joint_test <- function(full, reduced, terms) {
  if (!is.null(full$naive.var)) {
    b <- stats::coef(full)[terms]; V <- stats::vcov(full)[terms, terms]
    stat <- as.numeric(t(b) %*% solve(V) %*% b); type <- "Wald (robust)"
  } else {
    stat <- 2 * (full$loglik[2] - reduced$loglik[2]); type <- "LRT"
  }
  list(stat = stat, df = length(terms), p = stats::pchisq(stat, length(terms), lower.tail = FALSE),
       type = type)
}

# Fit the APOE model for one outcome / covariate set on a complete-case sample.
# exposure: character vector of terms; joint test is over all exposure columns in the fit.
fit_apoe <- function(cfg, data, outcome, set, exposure = c("e2_count", "e4_count"),
                     drop = NULL, add = NULL, timescale = "followup", cluster = NULL,
                     sample_vars = NULL, label = NULL) {
  covs <- model_covs(cfg, set, drop, add)
  if (timescale == "age") covs <- setdiff(covs, "age")
  vars <- unique(c(surv_vars(outcome, timescale), exposure, covs, cluster, sample_vars))
  dat <- complete_sample(data, vars)
  full <- survival::coxph(cox_formula(outcome, c(exposure, covs), timescale, cluster),
                          data = dat, ties = "efron", model = TRUE)
  red  <- survival::coxph(cox_formula(outcome, covs, timescale, cluster), data = dat, ties = "efron")
  exp_terms <- grep(paste0("^(", paste(exposure, collapse = "|"), ")"), names(stats::coef(full)), value = TRUE)
  jt <- joint_test(full, red, exp_terms)
  hr <- tidy_hr(full, exp_terms)
  hr[, `:=`(outcome = outcome, model = set, label = if (is.null(label)) set else label,
            n = nrow(dat), events = sum(dat[[paste0(outcome, "_event")]]),
            joint_stat = jt$stat, joint_df = jt$df, joint_p = jt$p, joint_type = jt$type)]
  data.table::setcolorder(hr, c("label", "outcome", "model", "n", "events", "term"))
  list(table = hr[], fit = full, data = dat)
}

term_label <- function(term) {
  lab <- c(e2_count = "Per e2 allele (vs e3)", e4_count = "Per e4 allele (vs e3)", e4_carrier = "E4 carrier vs non-E4")
  out <- lab[term]
  geno <- grepl("^apoe_genotype", term)
  out[geno] <- sub("^apoe_genotype", "", term[geno])
  unname(out)
}
