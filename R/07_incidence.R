# Purpose: Genotype-specific incidence rates of CKD progression (SAP Sec. 6.2; Figure 1).
# Input(s): data/derived/analysis.rds
# Output(s): output/tables/fig1_incidence_by_genotype.csv; output/figures/fig1_incidence_by_genotype.pdf/.png
# Key methods/packages: events / person-years x 1,000; exact (Garwood) Poisson 95% CI; ggplot2
# Notes: Person-time counted from entry (index, or consent date with delayed entry) to exit.

source("R/00_setup.R")
log_msg("== 07 incidence ==")
a <- load_derived("analysis", cfg)
a <- a[!is.na(primary_tstop)]
a[, py := primary_tstop - tstart]

glev <- unlist(cfg$apoe$genotype_order)
rate_row <- function(dat, lab) {
  ev <- sum(dat$primary_event); pt <- sum(dat$py)
  ci <- pois_exact_ci(ev, pt)
  data.table::data.table(genotype = lab, n = nrow(dat), events = ev, person_years = round(pt, 1),
                         rate_per_1000py = 1000 * ci$rate, lcl = 1000 * ci$lcl, ucl = 1000 * ci$ucl)
}
inc <- rbind(
  data.table::rbindlist(lapply(glev, function(g) rate_row(a[apoe_genotype == g], g))),
  rate_row(a, "Overall"))
inc[, `rate (95% CI)` := fmt_ci(rate_per_1000py, lcl, ucl, 1)]
write_table(inc, "fig1_incidence_by_genotype", cfg)

pd <- inc[genotype != "Overall"]
pd[, genotype := factor(genotype, levels = glev)]
p <- ggplot(pd, aes(x = genotype, y = rate_per_1000py)) +
  geom_hline(yintercept = inc[genotype == "Overall", rate_per_1000py], linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.15, colour = "grey30") +
  geom_point(size = 2.6, colour = "#1f3b6f") +
  geom_text(aes(y = ucl, label = paste0("n = ", prettyNum(n, big.mark = ","), "\n", prettyNum(events, big.mark = ","), " events")),
            vjust = -0.3, size = 2.8, colour = "grey30", lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = "APOE genotype", y = "CKD progression per 1,000 person-years (95% CI)",
       caption = "Dashed line: overall rate. Genotypes ordered from most e2 to most e4.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
save_figure(p, "fig1_incidence_by_genotype", cfg)
log_msg("07 done")
