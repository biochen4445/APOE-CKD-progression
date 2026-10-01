# Purpose: Table 1 -- index characteristics by APOE genotype and by CKD progression status (SAP Sec. 6.1).
# Input(s): data/derived/analysis.rds
# Output(s): output/tables/table1a_by_genotype.csv, table1b_by_outcome.csv
# Key methods/packages: mean (SD) or median [IQR]; n (%); standardized mean difference (SMD)
# Notes: Table 1a SMD = largest absolute SMD of any genotype vs e3/e3. Table 1b SMD = event vs no event.

source("R/00_setup.R")
log_msg("== 06 Table 1 ==")
a <- load_derived("analysis", cfg)
a[, male := as.integer(sex == "M")]

vars <- data.table::data.table(
  var   = c("age", "male", "bmi", "diabetes", "sbp", "dbp", "antihtn", "chd",
            "index_cr", "egfr_index", "hdl", "ldl", "tg", "fu_years"),
  label = c("Age, y", "Male", "BMI, kg/m2", "Diabetes", "Systolic BP, mm Hg", "Diastolic BP, mm Hg",
            "Antihypertensive use", "Coronary heart disease", "Serum creatinine, mg/dL",
            "eGFR (CKD-EPI 2021), mL/min/1.73 m2", "HDL cholesterol, mg/dL", "LDL cholesterol, mg/dL",
            "Triglycerides, mg/dL", "Follow-up, y"),
  type  = c("mean", "bin", "mean", "bin", "mean", "mean", "bin", "bin",
            "mean", "mean", "mean", "mean", "median", "median"))

summ_var <- function(x, type) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  switch(type,
    mean   = sprintf("%.1f (%.1f)", mean(x), stats::sd(x)),
    median = sprintf("%.1f [%.1f-%.1f]", stats::median(x), stats::quantile(x, .25), stats::quantile(x, .75)),
    bin    = sprintf("%s (%.1f)", format(sum(x), big.mark = ","), 100 * mean(x)))
}
smd <- function(x1, x2, type) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  if (type == "bin") {
    p1 <- mean(x1); p2 <- mean(x2); den <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
    if (den == 0) return(0) else return((p1 - p2) / den)
  }
  (mean(x1) - mean(x2)) / sqrt((stats::var(x1) + stats::var(x2)) / 2)
}
make_table <- function(data, group, levels, smd_fun) {
  rows <- lapply(seq_len(nrow(vars)), function(i) {
    v <- vars$var[i]; t <- vars$type[i]
    cells <- c(list(Overall = summ_var(data[[v]], t)),
               stats::setNames(lapply(levels, function(l) summ_var(data[get(group) == l][[v]], t)), levels))
    data.table::as.data.table(c(list(Characteristic = vars$label[i]), cells,
                                list(n_missing = sum(is.na(data[[v]])), SMD = round(smd_fun(v, t), 3))))
  })
  head <- data.table::as.data.table(c(list(Characteristic = "N"),
    list(Overall = format(nrow(data), big.mark = ",")),
    stats::setNames(lapply(levels, function(l) format(data[get(group) == l, .N], big.mark = ",")), levels),
    list(n_missing = NA, SMD = NA)))
  rbind(head, data.table::rbindlist(rows))
}

glev <- levels(a$apoe_genotype)
glev <- unlist(cfg$apoe$genotype_order)[unlist(cfg$apoe$genotype_order) %in% glev]
t1a <- make_table(a, "apoe_genotype", glev, function(v, t) {
  ref <- a[apoe_genotype == cfg$apoe$reference_genotype][[v]]
  s <- vapply(setdiff(glev, cfg$apoe$reference_genotype),
              function(l) smd(a[apoe_genotype == l][[v]], ref, t), numeric(1))
  s[which.max(abs(s))]
})
data.table::setnames(t1a, "SMD", "max_abs_SMD_vs_e3e3")

a[, progression := data.table::fifelse(primary_event == 1, "CKD progression", "No progression")]
t1b <- make_table(a, "progression", c("CKD progression", "No progression"), function(v, t) {
  smd(a[primary_event == 1][[v]], a[primary_event == 0][[v]], t)
})

write_table(t1a, "table1a_by_genotype", cfg)
write_table(t1b, "table1b_by_outcome", cfg)
log_msg("06 done")
