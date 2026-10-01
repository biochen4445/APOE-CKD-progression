# Purpose: Read raw EMR / genotype extracts, check them against the expected schema, standardise
#          dates and codes, and save one list object for the later stages.
# Input(s): <raw dir>/{person,apoe_genotype,labs,diagnoses,procedures,medications,vitals}.csv
#           optional: catastrophic.csv (id, date, category)
#           Schema: data/README.md
# Output(s): data/derived/inputs.rds; output/tables/s01_input_summary.csv, s01_data_checks.csv
# Key methods/packages: data.table
# Notes: Creatinine must be in mg/dL; the script stops if values look like umol/L.

source("R/00_setup.R")
log_msg("== 01 import and validate ==")

schema <- list(
  person        = c("id", "sex", "birth_date", "last_visit_date", "death_date"),
  apoe_genotype = c("id", "rs429358", "rs7412"),
  labs          = c("id", "date", "test", "value", "setting"),
  diagnoses     = c("id", "encounter_id", "date", "setting", "code_system", "code", "dx_rank"),
  procedures    = c("id", "encounter_id", "date", "setting", "code_system", "code"),
  medications   = c("id", "date", "atc"),
  vitals        = c("id", "date", "setting", "sbp", "dbp", "height_cm", "weight_kg")
)
date_cols <- c("birth_date", "last_visit_date", "death_date", "consent_date", "date", "discharge_date")

d <- list()
for (nm in names(schema)) {
  x <- read_raw(cfg, nm)
  miss <- setdiff(schema[[nm]], names(x))
  stop_if(length(miss) > 0, nm, ": missing required column(s): ", paste(miss, collapse = ", "))
  for (dc in intersect(date_cols, names(x))) x[, (dc) := data.table::as.IDate(get(dc))]
  if ("code" %in% names(x)) x[, code := norm_code(code)]
  if ("atc" %in% names(x)) x[, atc := norm_code(atc)]
  if ("death_cause_code" %in% names(x)) x[, death_cause_code := norm_code(death_cause_code)]
  if ("setting" %in% names(x)) x[, setting := toupper(setting)]
  if ("test" %in% names(x)) x[, test := tolower(test)]
  d[[nm]] <- x
}
cat_ill <- read_raw(cfg, "catastrophic", required = FALSE)
if (!is.null(cat_ill)) {
  cat_ill[, date := data.table::as.IDate(date)]
  d$catastrophic <- cat_ill
}

# ---- Checks -----------------------------------------------------------------------
chk <- list()
add_chk <- function(table, check, n_bad, action = "") {
  chk[[length(chk) + 1]] <<- data.table::data.table(table = table, check = check, n = n_bad, action = action)
}
add_chk("person", "duplicated id", sum(duplicated(d$person$id)), "stop if > 0")
stop_if(any(duplicated(d$person$id)), "person.csv has duplicated id")
add_chk("person", "sex not M/F", sum(!d$person$sex %in% c("M", "F")), "excluded at cohort step E3")
add_chk("person", "death_date before birth_date", d$person[death_date < birth_date, .N], "review")
add_chk("person", "last_visit_date after death_date + 1 d", d$person[last_visit_date > death_date + 1L, .N], "review")
add_chk("apoe_genotype", "id not in person", sum(!d$apoe_genotype$id %in% d$person$id), "dropped at cohort step I1")
cr <- d$labs[test == "creatinine"]
stop_if(nrow(cr) == 0, "labs has no rows with test == 'creatinine'")
stop_if(stats::median(cr$value, na.rm = TRUE) > 20,
        "Median creatinine > 20: values look like umol/L. Convert to mg/dL (divide by 88.4).")
add_chk("labs", "creatinine <= 0 or > 20 mg/dL", cr[value <= 0 | value > 20, .N], "<= 0 dropped; review > 20")
add_chk("labs", "setting not OPD/IPD/ER", sum(!d$labs$setting %in% c("OPD", "IPD", "ER")), "review")
known_tests <- c("creatinine", "hdl", "ldl", "tg", "hba1c", "glucose_fasting", "glucose_random")
add_chk("labs", "unrecognised test names (ignored)", sum(!d$labs$test %in% known_tests), "review")
add_chk("diagnoses", "code_system not ICD9CM/ICD10CM",
        sum(!d$diagnoses$code_system %in% c("ICD9CM", "ICD10CM")), "review")
add_chk("diagnoses", "ICD9CM code dated >= 2016-01-01",
        d$diagnoses[code_system == "ICD9CM" & date >= as.IDate("2016-01-01"), .N], "review (Taiwan switched 2016)")
add_chk("diagnoses", "IPD rows without encounter_id", d$diagnoses[setting == "IPD" & is.na(encounter_id), .N],
        "cannot be grouped into admissions")
add_chk("diagnoses", "IPD rows without discharge_date",
        if ("discharge_date" %in% names(d$diagnoses)) d$diagnoses[setting == "IPD" & is.na(discharge_date), .N] else nrow(d$diagnoses[setting == "IPD"]),
        "admission date used as event date")
checks <- data.table::rbindlist(chk)

summ <- data.table::rbindlist(lapply(names(d), function(nm) {
  x <- d[[nm]]
  dc <- intersect(c("date", "birth_date"), names(x))[1]
  rng <- if (is.na(dc)) c(NA, NA) else as.character(range(x[[dc]], na.rm = TRUE))
  data.table::data.table(table = nm, rows = nrow(x), ids = data.table::uniqueN(x$id),
                         date_min = rng[1], date_max = rng[2])
}))

d$labs <- d$labs[!(test == "creatinine" & (is.na(value) | value <= 0))]
save_derived(d, "inputs", cfg)
write_table(summ, "s01_input_summary", cfg)
write_table(checks, "s01_data_checks", cfg)
log_msg("01 done: ", nrow(d$person), " persons")
