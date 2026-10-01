# Purpose: Call APOE genotypes from rs429358 / rs7412 and report genotype QC (SAP Sec. 3).
# Input(s): data/derived/inputs.rds
# Output(s): data/derived/apoe.rds; output/tables/s02_apoe_qc.csv, s02_snp_hwe.csv,
#            s02_apoe_genotype_hwe.csv
# Key methods/packages: hard calls from dosages (tolerance in config), HWE exact test per SNP,
#            3-allele chi-square HWE across the six genotypes
# Notes: Allele coding (which allele is counted) must be verified against the array manifest.
#        Double heterozygotes are called e2/e4 by convention; e1-implying combinations are excluded.

source("R/00_setup.R")
log_msg("== 02 APOE genotype QC ==")
d <- load_derived("inputs", cfg)
g <- data.table::copy(d$apoe_genotype)
ac <- cfg$apoe

g[, n4 := hardcall(rs429358, ac$hardcall_tolerance)]
g[, n2 := hardcall(rs7412, ac$hardcall_tolerance)]
fail_info <- c(rs429358 = 0L, rs7412 = 0L)
if ("info_rs429358" %in% names(g)) {
  fail_info["rs429358"] <- g[info_rs429358 < ac$min_info, .N]
  g[info_rs429358 < ac$min_info, n4 := NA]
}
if ("info_rs7412" %in% names(g)) {
  fail_info["rs7412"] <- g[info_rs7412 < ac$min_info, .N]
  g[info_rs7412 < ac$min_info, n2 := NA]
}
g <- cbind(g, call_apoe(g$n4, g$n2))
g[, apoe_qc_pass := !is.na(apoe_genotype)]

qc <- data.table::data.table(
  item = c("Genotyped participants",
           "rs429358 not hard-callable (dosage or INFO)", "rs7412 not hard-callable (dosage or INFO)",
           "rs429358 INFO below threshold", "rs7412 INFO below threshold",
           "Combination implies e1 (excluded)", "Double heterozygote called e2/e4",
           "APOE genotype passing QC"),
  n = c(nrow(g), g[is.na(n4), .N], g[is.na(n2), .N], fail_info["rs429358"], fail_info["rs7412"],
        g[e1_implied == TRUE, .N], g[n4 == 1 & n2 == 1, .N], g[apoe_qc_pass == TRUE, .N]))
qc[, pct := round(100 * n / nrow(g), 2)]

snp_hwe <- data.table::rbindlist(lapply(c("n4", "n2"), function(v) {
  x <- g[!is.na(get(v)), get(v)]
  h <- c(sum(x == 0), sum(x == 1), sum(x == 2))
  data.table::data.table(snp = ifelse(v == "n4", "rs429358", "rs7412"),
                         counted_allele = ifelse(v == "n4", ac$rs429358_counted_allele, ac$rs7412_counted_allele),
                         n = length(x), n_0 = h[1], n_1 = h[2], n_2 = h[3],
                         counted_allele_freq = (h[2] + 2 * h[3]) / (2 * length(x)),
                         hwe_exact_p = hwe_exact(h[1], h[2], h[3]))
}))

ok <- g[apoe_qc_pass == TRUE]
hw <- hwe_apoe_chisq(ok$apoe_genotype)
geno_tab <- data.table::data.table(
  genotype = c("e2/e2", "e2/e3", "e3/e3", "e2/e4", "e3/e4", "e4/e4"),
  observed = hw$observed, expected = round(hw$expected, 1))
geno_tab[, pct := round(100 * observed / sum(observed), 2)]
geno_tab <- rbind(geno_tab, data.table::data.table(
  genotype = c("allele freq e2", "allele freq e3", "allele freq e4", "HWE chi-square (3 df) P"),
  observed = c(round(as.numeric(hw$allele_freq), 4), signif(hw$p, 3)), expected = NA, pct = NA))

pcs <- grep("^pc[0-9]+$", names(g), value = TRUE)
apoe <- g[, c("id", "n4", "n2", "apoe_genotype", "e1_implied", "e2_count", "e4_count", "apoe_qc_pass", pcs), with = FALSE]
save_derived(apoe, "apoe", cfg)
write_table(qc, "s02_apoe_qc", cfg)
write_table(snp_hwe, "s02_snp_hwe", cfg)
write_table(geno_tab, "s02_apoe_genotype_hwe", cfg)
log_msg("02 done: ", ok[, .N], " pass QC; e2 freq ", round(hw$allele_freq["e2"], 3),
        ", e4 freq ", round(hw$allele_freq["e4"], 3))
