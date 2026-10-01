d <- function(x) as.IDate(x)
x <- data.table(
  id   = c("a", "a", "a", "b", "b", "b", "c", "c", "d"),
  date = d(c("2020-01-01", "2020-02-01", "2020-06-01",   # a: rise confirmed 121 d later
             "2020-01-01", "2020-03-01", "2020-07-01",   # b: rise, next value >= 90 d is normal
             "2020-01-01", "2020-01-20",                 # c: rise, no value >= 90 d later
             "2020-05-01")),                             # d: single rise
  flag = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))

test_that("confirmation by the next value >= 90 days later", {
  ev <- detect_sustained(x, 90, "next", TRUE)
  expect_equal(ev$id, "a")
  expect_equal(ev$event_date, d("2020-02-01"))
})

test_that("rule 'any' accepts a qualifying value later than the next one", {
  y <- rbind(x, data.table(id = "b", date = d("2020-12-01"), flag = TRUE))
  ev <- detect_sustained(y, 90, "any", TRUE)
  expect_setequal(ev$id, c("a", "b"))
  expect_equal(ev[id == "b", event_date], d("2020-01-01"))
})

test_that("without confirmation the first flagged value is the event", {
  ev <- detect_sustained(x, 90, "next", FALSE)
  expect_setequal(ev$id, c("a", "b", "c", "d"))
  expect_equal(ev[id == "a", event_date], d("2020-02-01"))
})

test_that("code lists match on prefixes within the right coding system", {
  cl <- data.table(list = "l", code_system = c("ICD10CM", "ICD9CM"), code_prefix = c("N18", "585"))
  expect_equal(match_codelist(c("N183", "N183", "5853", "I10"), c("ICD10CM", "ICD9CM", "ICD9CM", "ICD10CM"), "l", cl),
               c(TRUE, FALSE, TRUE, FALSE))
})
