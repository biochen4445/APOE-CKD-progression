# Purpose: Call APOE genotypes from rs429358 and rs7412, and Hardy-Weinberg tests.
#
# Input coding: n4 = count (0/1/2) of the rs429358 allele defining e4 (C);
#               n2 = count (0/1/2) of the rs7412 allele defining e2 (T).
# Haplotypes: e2 = (rs429358 T, rs7412 T); e3 = (T, C); e4 = (C, C); e1 = (C, T).
# Unphased double heterozygote (n4 = 1, n2 = 1) is e2/e4 or e1/e3; by convention it is
# called e2/e4 because e1 is extremely rare. Combinations that require e1 on either
# haplotype (n4 + n2 > 2, or n4 = 2 with n2 >= 1, etc.) are flagged and set to NA.

hardcall <- function(dosage, tol = 0.1) {
  r <- round(dosage)
  r[abs(dosage - r) > tol | r < 0 | r > 2] <- NA
  as.integer(r)
}

call_apoe <- function(n4, n2) {
  key <- paste(n4, n2)
  map <- c("0 0" = "e3/e3", "1 0" = "e3/e4", "2 0" = "e4/e4",
           "0 1" = "e2/e3", "0 2" = "e2/e2", "1 1" = "e2/e4")
  geno <- unname(map[key])
  both_called <- !is.na(n4) & !is.na(n2)
  data.table::data.table(
    apoe_genotype = geno,
    e1_implied    = both_called & is.na(geno),
    e2_count      = ifelse(is.na(geno), NA_integer_, as.integer(n2)),
    e4_count      = ifelse(is.na(geno), NA_integer_, as.integer(n4))
  )
}

# Hardy-Weinberg exact test for one bi-allelic SNP (Wigginton, Cutler & Abecasis 2005)
hwe_exact <- function(n_hom1, n_het, n_hom2) {
  obs_homr <- min(n_hom1, n_hom2); obs_homc <- max(n_hom1, n_hom2)
  rare <- 2 * obs_homr + n_het
  n <- n_het + obs_homc + obs_homr
  if (n == 0 || rare == 0) return(1)
  probs <- numeric(rare + 1)
  mid <- floor(rare * (2 * n - rare) / (2 * n))
  if ((rare %% 2) != (mid %% 2)) mid <- mid + 1
  probs[mid + 1] <- 1
  hets <- mid; homr <- (rare - mid) / 2; homc <- n - hets - homr
  while (hets > 1) {
    probs[hets - 1] <- probs[hets + 1] * hets * (hets - 1) / (4 * (homr + 1) * (homc + 1))
    hets <- hets - 2; homr <- homr + 1; homc <- homc + 1
  }
  hets <- mid; homr <- (rare - mid) / 2; homc <- n - hets - homr
  while (hets <= rare - 2) {
    probs[hets + 3] <- probs[hets + 1] * 4 * homr * homc / ((hets + 2) * (hets + 1))
    hets <- hets + 2; homr <- homr - 1; homc <- homc - 1
  }
  probs <- probs / sum(probs)
  min(1, sum(probs[probs <= probs[n_het + 1] * (1 + 1e-7)]))
}

# Chi-square HWE test for the 3-allele APOE system (6 genotypes, 3 df)
hwe_apoe_chisq <- function(genotypes) {
  g <- genotypes[!is.na(genotypes)]
  n <- length(g)
  alle <- unlist(strsplit(g, "/"))
  p <- table(factor(alle, levels = c("e2", "e3", "e4"))) / (2 * n)
  levs <- c("e2/e2", "e2/e3", "e3/e3", "e2/e4", "e3/e4", "e4/e4")
  expct <- n * c(p["e2"]^2, 2 * p["e2"] * p["e3"], p["e3"]^2,
                 2 * p["e2"] * p["e4"], 2 * p["e3"] * p["e4"], p["e4"]^2)
  obs <- as.numeric(table(factor(g, levels = levs)))
  stat <- sum((obs - expct)^2 / expct)
  list(allele_freq = p, observed = obs, expected = as.numeric(expct),
       chisq = stat, df = 3, p = stats::pchisq(stat, 3, lower.tail = FALSE))
}
