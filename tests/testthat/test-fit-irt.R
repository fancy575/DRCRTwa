test_that("DRCRT_WA fits an individually randomized trial", {
  dat <- wa_irt_data()
  fit <- DRCRT_WA(
    survival::Surv(time, status) ~ z1 + z3,
    data = dat,
    id = id,
    treatment = treatment,
    censoring_formula = ~ z1 + z2,
    terminal_formula = ~ z1 + z2 + z3,
    terminal_event = 3,
    event_weights = c("1" = 1, "2" = 2),
    times = c(0.5, 1),
    control = DRCRTwa_control(
      grid_per_unit = 20, min_grid = 20, max_grid = 100,
      keep_nuisance = TRUE, retain_data = TRUE
    )
  )
  expect_s3_class(fit, "DRCRTwa")
  expect_identical(fit$analysis$design, "IRT")
  expect_identical(fit$analysis$estimands, "individual-average")
  expect_equal(nrow(fit$estimates), 2L * 3L)
  expect_equal(nrow(fit$arm_estimates), 2L * 2L * 3L)
  expect_true(all(is.finite(
    fit$estimates$estimate[fit$estimates$method == "dr"]
  )))
  expect_s3_class(summary(fit, times = 0.75), "summary.DRCRTwa")
  expect_s3_class(plot(fit), "ggplot")
})

test_that("an IRT rejects the cluster-average estimand", {
  dat <- wa_irt_data()
  expect_error(
    DRCRT_WA(
      survival::Surv(time, status) ~ z1,
      data = dat, id = id, treatment = treatment,
      censoring_formula = ~ z2,
      terminal_event = 3,
      event_weights = c("1" = 1, "2" = 1),
      times = 1,
      estimand = "cluster-average"
    ),
    "only the individual-average"
  )
})
