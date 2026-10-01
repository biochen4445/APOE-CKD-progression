test_that("APOE genotypes are called from rs429358 (e4 allele) and rs7412 (e2 allele) counts", {
  n4 <- c(0, 1, 2, 0, 0, 1, 2, 1, 2, NA)
  n2 <- c(0, 0, 0, 1, 2, 1, 1, 2, 2, 0)
  r <- call_apoe(n4, n2)
  expect_equal(r$apoe_genotype, c("e3/e3", "e3/e4", "e4/e4", "e2/e3", "e2/e2", "e2/e4", NA, NA, NA, NA))
  expect_equal(r$e1_implied, c(rep(FALSE, 6), TRUE, TRUE, TRUE, FALSE))
  expect_equal(r$e2_count[1:6], c(0L, 0L, 0L, 1L, 2L, 1L))
  expect_equal(r$e4_count[1:6], c(0L, 1L, 2L, 0L, 0L, 1L))
})

test_that("hard calls reject ambiguous dosages", {
  expect_equal(hardcall(c(0.02, 0.95, 1.5, 2.08, -0.3), 0.1), c(0L, 1L, NA, 2L, NA))
})

test_that("HWE exact test behaves at the extremes", {
  expect_gt(hwe_exact(25, 50, 25), 0.5)        # perfect HWE
  expect_lt(hwe_exact(50, 0, 50), 1e-10)       # no heterozygotes
  expect_equal(hwe_exact(100, 0, 0), 1)        # monomorphic
})
