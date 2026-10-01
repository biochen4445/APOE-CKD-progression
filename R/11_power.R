# Purpose: Power (SAP Sec. 6.6): events needed per allele HR (Schoenfeld) and detectable HR given
#          the observed number of events.
# Input(s): config power block; data/derived/analysis.rds (observed events and allele frequencies, if present)
# Output(s): output/tables/s11_power.csv
# Key methods/packages: Schoenfeld (1983): D = (z_{1-a/2} + z_{power})^2 / (var(X) * log(HR)^2),
#            var(X) = 2p(1-p) for an additive allele count under HWE.
# Notes: Ignores adjustment for the other allele and covariates (small effect for APOE).

source("R/00_setup.R")
log_msg("== 11 power ==")
pw <- cfg$power
z <- stats::qnorm(1 - pw$alpha / 2) + stats::qnorm(pw$power)
events_needed <- function(hr, p) z^2 / (2 * p * (1 - p) * log(hr)^2)
detectable_hr <- function(D, p) exp(c(-1, 1) * z / sqrt(D * 2 * p * (1 - p)))

grid <- data.table::CJ(allele = c("e2", "e4"), hr = unlist(pw$hr_grid))
grid[, freq := ifelse(allele == "e2", pw$freq_e2, pw$freq_e4)]
grid[, events_required := ceiling(events_needed(hr, freq))]
grid[, basis := "SAP assumed allele frequency"]

obs <- NULL
f <- file.path(cfg$paths$derived, "analysis.rds")
if (file.exists(f)) {
  a <- readRDS(f)
  D <- sum(a$primary_event, na.rm = TRUE)
  pe2 <- mean(a$e2_count) / 2; pe4 <- mean(a$e4_count) / 2
  obs <- data.table::data.table(
    allele = c("e2", "e4"), freq = c(pe2, pe4), observed_events = D,
    detectable_hr_lower = c(detectable_hr(D, pe2)[1], detectable_hr(D, pe4)[1]),
    detectable_hr_upper = c(detectable_hr(D, pe2)[2], detectable_hr(D, pe4)[2]),
    basis = "Observed cohort (data_source in config)")
}
write_table(rbind(grid, obs, fill = TRUE), "s11_power", cfg)
log_msg("11 done: events for e4 HR 0.85 = ", grid[allele == "e4" & hr == 0.85, events_required])
