library(DRCRTwa)

dat <- wa_crt_data()
fit <- DRCRT_WA(
  survival::Surv(time, status) ~ z1 + z3,
  data = dat,
  id = id,
  treatment = treatment,
  cluster = cluster,
  censoring_formula = survival::Surv(time, status) ~ z1 + z2,
  terminal_event = 3,
  event_weights = c(`1` = 1, `2` = 2),
  times = c(0.5, 1, 2)
)
summary(fit, times = c(0.5, 1, 2))
plot(fit)
