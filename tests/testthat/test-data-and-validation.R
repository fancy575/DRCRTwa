test_that("example data have the required event-history structure", {
  irt <- wa_irt_data()
  expect_true(all(c("id", "treatment", "time", "status", "z1", "z2", "z3") %in%
                    names(irt)))
  expect_setequal(unique(irt$status), 0:3)
  expect_true(all(vapply(split(irt$status, irt$id), function(z) {
    sum(z %in% c(0L, 3L)) == 1L
  }, logical(1L))))

  crt <- wa_crt_data()
  expect_true("cluster" %in% names(crt))
  key <- interaction(crt$cluster, crt$id, drop = TRUE)
  expect_true(all(vapply(split(crt$status, key), function(z) {
    sum(z %in% c(0L, 3L)) == 1L
  }, logical(1L))))
  expect_lt(length(unique(crt$id)), length(unique(key)))
})

test_that("event-weight names and dimensions are checked", {
  dat <- wa_irt_data()
  expect_error(
    DRCRT_WA(
      survival::Surv(time, status) ~ z1,
      data = dat, id = id, treatment = treatment,
      censoring_formula = ~ z2,
      terminal_event = 3,
      event_weights = c("1" = 1),
      times = 1
    ),
    "exactly match"
  )
  expect_error(
    DRCRT_WA(
      survival::Surv(time, status) ~ z1,
      data = dat, id = id, treatment = treatment,
      censoring_formula = ~ z2,
      terminal_event = 3,
      event_weights = c(1, 1, 1),
      times = 1
    ),
    "must have length"
  )
})
