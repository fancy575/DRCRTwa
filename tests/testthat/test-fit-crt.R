test_that("DRCRT_WA fits both CRT estimands", {
  dat <- wa_crt_data()
  fit <- DRCRT_WA(
    survival::Surv(time, status) ~ z1 + z3,
    data = dat,
    id = id,
    treatment = treatment,
    cluster = cluster,
    censoring_formula = ~ z1 + z2,
    terminal_formula = ~ z1 + z2 + z3,
    terminal_event = 3,
    event_weights = c("1" = 1, "2" = 2),
    times = 1,
    control = DRCRTwa_control(
      grid_per_unit = 20, min_grid = 20, max_grid = 100,
      keep_nuisance = FALSE, retain_data = TRUE
    )
  )
  expect_s3_class(fit, "DRCRTwa")
  expect_identical(fit$analysis$design, "CRT")
  expect_setequal(
    fit$analysis$estimands,
    c("individual-average", "cluster-average")
  )
  expect_equal(nrow(fit$estimates), 2L * 3L)
  expect_true(all(is.finite(
    fit$estimates$estimate[fit$estimates$method == "dr"]
  )))
  expect_equal(fit$analysis$n_clusters, 18L)
})
