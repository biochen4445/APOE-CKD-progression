# Purpose: Generate synthetic EMR + APOE data that match the input schema (data/README.md), so the
#          pipeline can be run, tested and shared without patient data.
# Input(s): none (seed fixed)
# Output(s): data/simulated/{person,apoe_genotype,labs,diagnoses,procedures,medications,vitals}.csv
# Key methods/packages: data.table
# Notes: True data-generating effects: HR 1.15 per e2 allele and 0.85 per e4 allele on latent CKD
#        progression. Numbers from simulated data are NOT study results.
#        Includes deliberate edge cases: e1-implying genotypes, unclear dosages, transient inpatient
#        AKI spikes, AKI-principal admissions, baseline dialysis, unconfirmed creatinine rises.

suppressPackageStartupMessages(library(data.table))
set.seed(20260914)
N <- as.integer(Sys.getenv("SIM_N", 20000))
out_dir <- "data/simulated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
end_date <- as.IDate("2025-12-31")
icd10_start <- as.IDate("2016-01-01")
dx_sys <- function(date) ifelse(date < icd10_start, "ICD9CM", "ICD10CM")
px_sys <- function(date) ifelse(date < icd10_start, "ICD9PROC", "ICD10PCS")
pick <- function(date, icd9, icd10) ifelse(date < icd10_start, icd9, icd10)

# ---- Persons ------------------------------------------------------------------------
p <- data.table(id = sprintf("SIM%06d", seq_len(N)),
                sex = sample(c("M", "F"), N, replace = TRUE))
p[, birth_date := as.IDate("1940-01-01") + sample(0:(35 * 365), N, replace = TRUE)]
p[, first_visit := as.IDate("2005-01-01") + sample(0:(13 * 365), N, replace = TRUE)]
p[, last_visit_date := pmin(first_visit + as.integer(round(365.25 * (0.5 + rexp(N, 1 / 9)))), end_date)]
p[, age0 := as.numeric(first_visit - birth_date) / 365.25]
p[, family_id := sprintf("F%06d", sample(seq_len(round(0.8 * N)), N, replace = TRUE))]
p[, consent_date := as.IDate("2010-01-01") + sample(0:(10 * 365), N, replace = TRUE)]

# APOE haplotypes (East Asian-like frequencies, plus rare e1)
draw <- function(n) sample(c("e2", "e3", "e4", "e1"), n, replace = TRUE, prob = c(0.08, 0.8295, 0.09, 0.0005))
h1 <- draw(N); h2 <- draw(N)
p[, e2c := (h1 == "e2") + (h2 == "e2")]
p[, e4c := (h1 == "e4") + (h2 == "e4")]
p[, n4 := (h1 %in% c("e4", "e1")) + (h2 %in% c("e4", "e1"))]
p[, n2 := (h1 %in% c("e2", "e1")) + (h2 %in% c("e2", "e1"))]

# Risk factors
p[, dm := rbinom(N, 1, 0.12 + 0.04 * (sex == "M"))]
p[, htn := rbinom(N, 1, 0.30)]
p[, chd := rbinom(N, 1, 0.05)]
p[, base_cr := ifelse(sex == "M", 0.95, 0.72) * exp(rnorm(N, 0, 0.15))]
p[, height := ifelse(sex == "M", 170, 158) + rnorm(N, 0, 6)]
p[, weight := ifelse(sex == "M", 70, 58) + rnorm(N, 0, 10)]

# Latent CKD progression (years after first visit) and death
p[, lp := log(1.15) * e2c + log(0.85) * e4c + 0.7 * dm + 0.4 * htn + 0.3 * chd +
      0.2 * (sex == "M") + 0.04 * (age0 - 55)]
p[, t_prog := rexp(N, 0.03 * exp(lp))]
p[, t_death := rexp(N, 0.006 * exp(0.07 * (age0 - 55) + 0.4 * dm + 0.3 * (t_prog < 10)))]
p[, death_date := first_visit + as.integer(round(365.25 * t_death))]
p[death_date >= last_visit_date, death_date := NA]
p[!is.na(death_date), last_visit_date := death_date]
p[, fu_y := as.numeric(last_visit_date - first_visit) / 365.25]
p[, eskd := rbinom(N, 1, 0.10) == 1 & t_prog < fu_y - 1]
p[, t_eskd := t_prog + runif(N, 0.5, 3)]

# ---- Outpatient visits ---------------------------------------------------------------
nv <- rpois(N, 3 * p$fu_y) + 1L
v <- data.table(id = rep(p$id, nv))
v <- merge(v, p[, .(id, first_visit, last_visit_date, fu_y, base_cr, t_prog, sex, dm, htn, height, weight)], by = "id")
v[, date := first_visit + as.integer(round(runif(.N) * as.numeric(last_visit_date - first_visit)))]
v[, t := as.numeric(date - first_visit) / 365.25]

# Creatinine: slow drift; step increase after latent progression
v[, cr := base_cr * (1 + 0.004 * t) +
       ifelse(t > t_prog, 0.45 + 0.08 * (t - t_prog), 0)]
v[, cr := round(cr * exp(rnorm(.N, 0, 0.05)), 2)]
labs_cr <- v[runif(.N) < 0.55, .(id, date, test = "creatinine", value = cr, setting = "OPD")]
# a few isolated outpatient rises that are not confirmed later (tests the confirmation rule)
spk <- labs_cr[sample(.N, round(0.01 * .N))]
spk[, value := value + 0.6]
labs_cr <- rbind(labs_cr, spk[, date := date + 1L])

labs_other <- v[runif(.N) < 0.40]
labs_other <- rbind(
  labs_other[, .(id, date, test = "hdl", value = round(rnorm(.N, ifelse(sex == "M", 45, 55), 10)))],
  labs_other[, .(id, date, test = "ldl", value = round(rnorm(.N, 115, 30)))],
  labs_other[, .(id, date, test = "tg", value = round(exp(rnorm(.N, log(120) + 0.3 * dm, 0.45))))],
  labs_other[, .(id, date, test = "hba1c", value = round(ifelse(dm == 1, rnorm(.N, 7.4, 1), rnorm(.N, 5.6, 0.35)), 1))],
  labs_other[, .(id, date, test = "glucose_fasting", value = round(ifelse(dm == 1, rnorm(.N, 140, 25), rnorm(.N, 95, 10))))])
labs_other[, setting := "OPD"]

vit <- v[runif(.N) < 0.6, .(id, date, setting = "OPD",
                            sbp = round(125 + 12 * htn + rnorm(.N, 0, 12)),
                            dbp = round(78 + 6 * htn + rnorm(.N, 0, 8)),
                            height_cm = round(height + rnorm(.N, 0, 0.5), 1),
                            weight_kg = round(weight + rnorm(.N, 0, 1.5), 1))]
vit[runif(.N) < 0.15, c("height_cm", "weight_kg") := NA]

# Medications at visits
meds <- rbind(
  v[htn == 1 & runif(.N) < 0.85, .(id, date, atc = "C09AA05", days_supply = 28L)],
  v[dm == 1 & runif(.N) < 0.85, .(id, date, atc = "A10BA02", days_supply = 28L)],
  v[runif(.N) < 0.05, .(id, date, atc = "N02BE01", days_supply = 7L)])

# ---- Diagnoses -------------------------------------------------------------------------
dx_dm <- p[dm == 1, .(id, date = first_visit + as.integer(runif(.N, -3, 1) * 365))]
dx_dm[, `:=`(code = pick(date, "25000", "E119"))]
dx_chd <- p[chd == 1, .(id, date = first_visit - as.integer(runif(.N, 0, 5) * 365))]
dx_chd[, code := pick(date, "41401", "I2510")]
dx_ckd_prev <- p[runif(N) < 0.03, .(id, date = first_visit - as.integer(runif(.N, 0, 2) * 365))]
dx_ckd_prev[, code := pick(date, "5853", "N183")]
opd_dx <- rbind(dx_dm, dx_chd, dx_ckd_prev)
opd_dx[, `:=`(encounter_id = NA_character_, setting = "OPD", dx_rank = 1L, discharge_date = as.IDate(NA))]

# Inpatient admissions: background, CKD-related after progression, AKI-principal
adm_n <- rpois(N, 0.06 * p$fu_y)
adm <- data.table(id = rep(p$id, adm_n))
adm <- merge(adm, p[, .(id, first_visit, last_visit_date, base_cr, t_prog)], by = "id")
adm[, date := first_visit + as.integer(round(runif(.N) * as.numeric(last_visit_date - first_visit)))]
adm[, type := ifelse(runif(.N) < 0.15, "aki", "other")]
ckd_adm <- p[t_prog < fu_y & runif(N) < 0.35, .(id, first_visit, last_visit_date, base_cr, t_prog)]
ckd_adm[, date := pmin(first_visit + as.integer(365.25 * (t_prog + runif(.N, 0.2, 2))), last_visit_date)]
ckd_adm[, type := "ckd"]
adm <- rbind(adm, ckd_adm)
adm[, encounter_id := sprintf("A%07d", seq_len(.N))]
adm[, discharge_date := pmin(date + sample(2:10, .N, replace = TRUE), last_visit_date)]
adm[discharge_date < date, discharge_date := date]
ipd_dx <- rbind(
  adm[, .(id, encounter_id, date, discharge_date, dx_rank = 1L,
          code = data.table::fcase(type == "aki", pick(date, "5849", "N179"),
                                   type == "ckd", pick(date, "4019", "I10"),
                                   rep(TRUE, .N), pick(date, "486", "J189")))],
  adm[type %in% c("aki", "ckd") & runif(.N) < 0.8,
      .(id, encounter_id, date, discharge_date, dx_rank = 2L, code = pick(date, "5853", "N183"))])
ipd_dx[, setting := "IPD"]
# inpatient creatinine, with AKI spikes (must not count: outcome uses OPD values)
ipd_cr <- adm[, .(id, date = date + 1L, test = "creatinine",
                  value = round(base_cr * ifelse(type == "aki", 2.2, 1.05), 2), setting = "IPD")]

diagnoses <- rbind(opd_dx, ipd_dx, use.names = TRUE)

# ---- Dialysis: incident ESKD and a few prevalent (excluded at baseline) --------------------
esk <- p[eskd == TRUE]
dial <- esk[, .(date = first_visit + as.integer(365.25 * (t_eskd + seq(0, max(0, fu_y - t_eskd), by = 1 / 12)))),
            by = id]
dial <- merge(dial, p[, .(id, last_visit_date)], by = "id")[date <= last_visit_date]
prev <- p[runif(N) < 0.004, .(id, date = first_visit - 200L)]
dial <- rbind(dial[, .(id, date)], prev)
procedures <- dial[, .(id, encounter_id = NA_character_, date, setting = "OPD",
                       code_system = px_sys(date), code = pick(date, "3995", "5A1D70Z"))]
dial_dx <- dial[, .(id, encounter_id = NA_character_, date, discharge_date = as.IDate(NA), setting = "OPD",
                    dx_rank = 1L, code = pick(date, "V451", "Z992"))]
diagnoses <- rbind(diagnoses, dial_dx, use.names = TRUE)
diagnoses[, code_system := dx_sys(date)]
setcolorder(diagnoses, c("id", "encounter_id", "date", "discharge_date", "setting", "code_system", "code", "dx_rank"))

# ---- Genotype file -----------------------------------------------------------------------
g <- p[, .(id, rs429358 = n4 + rnorm(N, 0, 0.02), rs7412 = n2 + rnorm(N, 0, 0.02))]
bad <- sample(N, round(0.01 * N)); g[bad, rs429358 := rs429358 + rnorm(length(bad), 0, 0.4)]
g[sample(N, round(0.005 * N)), rs7412 := NA]
g[, `:=`(info_rs429358 = 0.97, info_rs7412 = 0.99)]
for (k in 1:10) g[, (paste0("pc", k)) := round(rnorm(N, 0, 0.01), 5)]
g[, `:=`(rs429358 = round(rs429358, 3), rs7412 = round(rs7412, 3))]
g <- g[-sample(N, round(0.02 * N))]   # 2% of EMR persons not genotyped

# ---- Write ---------------------------------------------------------------------------------
person <- p[, .(id, sex, birth_date, last_visit_date, death_date, consent_date, family_id)]
labs <- rbind(labs_cr, labs_other, ipd_cr)
fwrite(person, file.path(out_dir, "person.csv"))
fwrite(g, file.path(out_dir, "apoe_genotype.csv"))
fwrite(labs, file.path(out_dir, "labs.csv"))
fwrite(diagnoses, file.path(out_dir, "diagnoses.csv"))
fwrite(procedures, file.path(out_dir, "procedures.csv"))
fwrite(meds, file.path(out_dir, "medications.csv"))
fwrite(vit, file.path(out_dir, "vitals.csv"))
message(sprintf("Simulated data written to %s: %d persons, %d lab rows, %d diagnosis rows",
                out_dir, N, nrow(labs), nrow(diagnoses)))
