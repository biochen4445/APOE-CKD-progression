test_that("CKD-EPI 2021 reproduces reference values", {
  # At Scr = kappa the creatinine terms equal 1
  expect_equal(egfr_ckdepi2021(0.7, 50, "F"), 142 * 0.9938^50 * 1.012)
  expect_equal(egfr_ckdepi2021(0.9, 50, "M"), 142 * 0.9938^50)
  # 60-year-old man, Scr 1.0 mg/dL -> about 86 mL/min/1.73 m2
  expect_equal(round(egfr_ckdepi2021(1.0, 60, "M")), 86)
  expect_true(all(diff(egfr_ckdepi2021(c(0.6, 0.8, 1.2, 2), 55, "F")) < 0))
})

test_that("MDRD uses the configured constant and female factor", {
  expect_equal(egfr_mdrd(1, 1, "M", 175), 175)
  expect_equal(egfr_mdrd(1, 1, "F", 186.3), 186.3 * 0.742)
})
