.add_inference <- function(x, design, n_units, conf.level) {
  alpha <- 1 - conf.level
  if (identical(design, "CRT")) {
    df <- n_units - 2L
    if (df < 1L) {
      .drcrtwa_stop("CRT inference requires at least three randomized clusters.")
    }
    critical <- stats::qt(1 - alpha / 2, df = df)
    stat <- x$estimate / x$std_error
    p_value <- 2 * stats::pt(-abs(stat), df = df)
    distribution <- paste0("t(", df, ")")
  } else {
    df <- Inf
    critical <- stats::qnorm(1 - alpha / 2)
    stat <- x$estimate / x$std_error
    p_value <- 2 * stats::pnorm(-abs(stat))
    distribution <- "normal"
  }
  x$lower <- x$estimate - critical * x$std_error
  x$upper <- x$estimate + critical * x$std_error
  x$statistic <- stat
  x$p_value <- p_value
  x$conf.level <- conf.level
  attr(x, "df") <- df
  attr(x, "reference_distribution") <- distribution
  x
}

.name_one_nuisance_fit <- function(x, coefficient_names, horizon) {
  if (is.null(x)) return(NULL)
  if (length(x$beta)) names(x$beta) <- coefficient_names
  if (length(x$base_increment)) {
    x$grid_time <- seq_len(length(x$base_increment)) *
      horizon / length(x$base_increment)
  }
  x
}

.name_nuisance_fits <- function(x, prep, horizon) {
  if (is.null(x)) return(NULL)
  for (arm in c("control", "treatment")) {
    x$censoring[[arm]] <- .name_one_nuisance_fit(
      x$censoring[[arm]], colnames(prep$X_censoring), horizon)
    x$terminal[[arm]] <- .name_one_nuisance_fit(
      x$terminal[[arm]], colnames(prep$X_terminal), horizon)
    for (code in names(x$recurrent[[arm]])) {
      x$recurrent[[arm]][[code]] <- .name_one_nuisance_fit(
        x$recurrent[[arm]][[code]], colnames(prep$X_recurrent), horizon)
    }
  }
  x
}

.grid_for_times <- function(prep, times) {
  raw <- pmax(prep$control$min_grid,
              ceiling(prep$control$grid_per_unit * times))
  out <- as.integer(pmin(prep$control$max_grid, raw))
  if (any(raw > prep$control$max_grid)) {
    warning("The numerical integration grid was capped at `control$max_grid` for ",
            "one or more requested times.", call. = FALSE)
  }
  out
}

.fit_prepared_drcrtwa <- function(prep, times, conf.level,
                                  keep_influence = prep$control$keep_influence,
                                  keep_nuisance = prep$control$keep_nuisance) {
  times <- sort(unique(as.numeric(times)))
  n_grid <- .grid_for_times(prep, times)
  raw <- drcrtwa_fit_cpp(
    arm = as.integer(prep$arm),
    unit = as.integer(prep$unit),
    cluster_size = as.integer(prep$cluster_size),
    X_censoring = prep$X_censoring,
    X_terminal = prep$X_terminal,
    X_recurrent = prep$X_recurrent,
    final_time = as.numeric(prep$final_time),
    final_status = as.integer(prep$final_status),
    recurrent_times = unname(prep$recurrent_times),
    recurrent_codes = as.integer(prep$recurrent_codes),
    event_weights = as.numeric(prep$event_weights),
    terminal_code = as.integer(prep$terminal_code),
    horizons = times,
    n_grid = n_grid,
    prob0 = as.numeric(prep$prob0),
    prob1 = as.numeric(prep$prob1),
    target_weights = unname(prep$target_weights),
    estimand_names = prep$estimands,
    keep_influence = keep_influence,
    keep_nuisance = keep_nuisance
  )

  contrast_parts <- list()
  arm_parts <- list()
  diagnostic_parts <- list()
  influence <- if (keep_influence) vector("list", length(raw$horizons)) else NULL
  for (h in seq_along(raw$horizons)) {
    one <- raw$horizons[[h]]
    contrast_parts[[h]] <- .bind_rows_base(one$contrasts)
    arm_parts[[h]] <- .bind_rows_base(one$arms)
    diagnostic_parts[[h]] <- one$diagnostics
    if (keep_influence) {
      influence[[h]] <- one$influence
      names(influence)[[h]] <- format(times[[h]], digits = 15)
    }
  }
  contrasts <- .bind_rows_base(contrast_parts)
  arms <- .bind_rows_base(arm_parts)
  diagnostics <- .bind_rows_base(diagnostic_parts)

  contrasts <- .add_inference(contrasts, prep$design, prep$n_units, conf.level)
  arms <- .add_inference(arms, prep$design, prep$n_units, conf.level)
  arms$arm_label <- ifelse(
    arms$arm == 1L,
    unname(prep$treatment_labels[["treatment"]]),
    unname(prep$treatment_labels[["control"]])
  )
  diagnostics$arm_label <- ifelse(
    diagnostics$arm == 1L,
    unname(prep$treatment_labels[["treatment"]]),
    unname(prep$treatment_labels[["control"]])
  )
  diagnostics$event_code[
    diagnostics$model == "terminal" & is.na(diagnostics$event_code)
  ] <- prep$terminal_code

  ord_c <- order(contrasts$estimand, contrasts$method, contrasts$time)
  ord_a <- order(arms$estimand, arms$method, arms$arm, arms$time)
  contrasts <- contrasts[ord_c, , drop = FALSE]
  arms <- arms[ord_a, , drop = FALSE]
  rownames(contrasts) <- NULL
  rownames(arms) <- NULL

  list(
    estimates = contrasts,
    arm_estimates = arms,
    diagnostics = diagnostics,
    influence = influence,
    nuisance = .name_nuisance_fits(raw$nuisance, prep, max(times)),
    times = times,
    n_grid = n_grid
  )
}

#' Doubly robust while-alive estimation for randomized trials
#'
#' Fits the exposure-weighted while-alive estimator for an individually
#' randomized trial (IRT) or cluster-randomized trial (CRT). The recurrent-event
#' and terminal-event outcome working models are fitted arm specifically, as is
#' the censoring working model. The compiled engine uses Cox working models for
#' terminal events and censoring, LWYY marginal proportional-rate working models
#' for each recurrent-event type, censoring-martingale augmentation, and
#' independent-unit influence-function inference.
#'
#' The input data must be in long event-history form. Each participant has one
#' row for every recurrent event and exactly one closing row. Status zero denotes
#' censoring or administrative end of observation; positive nonterminal codes
#' denote recurrent-event types; and `terminal_event` denotes the terminal event.
#'
#' @param formula A two-sided formula such as
#'   `survival::Surv(time, status) ~ z1 + z3`. Its right side specifies the
#'   recurrent-event working model and, by default, the terminal-event working
#'   model.
#' @param data Long-format data frame.
#' @param id Unquoted participant-ID column or its name as a character string.
#' @param treatment Unquoted binary treatment column or its name.
#' @param cluster For a CRT, the unquoted cluster-ID column or its name. Leave
#'   `NULL` for an IRT.
#' @param censoring_formula Formula specifying the censoring working model, for
#'   example `survival::Surv(time, status) ~ z1 + z2`. Only its right side is
#'   used; the event coding is taken from `formula`.
#' @param terminal_formula Formula specifying the terminal-event working model.
#'   It defaults to `formula`.
#' @param terminal_event Positive integer status code for the terminal event.
#'   By default, the largest observed positive status code is used.
#' @param event_weights Numeric weights for recurrent-event types. A named vector
#'   such as `c(`1` = 1, `2` = 2)` is recommended. Names must exactly match the
#'   recurrent-event codes. Unnamed weights are matched to sorted codes.
#' @param times Positive analysis times. The default evaluates every observed
#'   positive event time and the largest follow-up time.
#' @param estimand For a CRT, `"individual-average"`,
#'   `"cluster-average"`, or `"both"`; the default is both. An IRT has only
#'   the individual-average estimand.
#' @param treatment_levels Optional two-element vector ordered as control then
#'   treatment.
#' @param treatment_probability Known probability Pr(A=1): one scalar, one value
#'   per participant, one per independent unit, one per long-format row, or the
#'   name of a data column. If omitted, the observed independent-unit allocation
#'   fraction is used.
#' @param conf.level Confidence level for pointwise Wald intervals.
#' @param control A list returned by [DRCRTwa_control()].
#'
#' @return An object of class `DRCRTwa` containing treatment contrasts, arm
#'   trajectories, pointwise confidence intervals, diagnostics, nuisance fits,
#'   and optionally independent-unit influence values.
#'
#' @examples
#' dat <- wa_irt_data()
#' fit <- DRCRT_WA(
#'   survival::Surv(time, status) ~ z1 + z3,
#'   data = dat,
#'   id = id,
#'   treatment = treatment,
#'   censoring_formula = survival::Surv(time, status) ~ z1 + z2,
#'   terminal_event = 3,
#'   event_weights = c(`1` = 1, `2` = 2),
#'   times = c(0.5, 1, 2)
#' )
#' summary(fit)
#'
#' @export
DRCRT_WA <- function(formula,
                     data,
                     id,
                     treatment,
                     cluster = NULL,
                     censoring_formula = formula,
                     terminal_formula = formula,
                     terminal_event = NULL,
                     event_weights = NULL,
                     times = NULL,
                     estimand = NULL,
                     treatment_levels = NULL,
                     treatment_probability = NULL,
                     conf.level = 0.95,
                     control = DRCRTwa_control()) {
  call <- match.call()
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    .drcrtwa_stop("`conf.level` must lie strictly between zero and one.")
  }
  if (!is.list(control)) {
    .drcrtwa_stop("`control` must be created by `DRCRTwa_control()`.")
  }
  defaults <- DRCRTwa_control()
  missing_control <- setdiff(names(defaults), names(control))
  if (length(missing_control)) control[missing_control] <- defaults[missing_control]

  id_name <- .resolve_column_name(substitute(id), data, "id", parent.frame())
  treatment_name <- .resolve_column_name(
    substitute(treatment), data, "treatment", parent.frame())
  cluster_expr <- substitute(cluster)
  cluster_name <- if (identical(cluster_expr, quote(NULL))) {
    NULL
  } else {
    .resolve_column_name(cluster_expr, data, "cluster", parent.frame())
  }

  prep <- .prepare_drcrtwa_data(
    formula = formula,
    censoring_formula = censoring_formula,
    terminal_formula = terminal_formula,
    data = data,
    id_name = id_name,
    treatment_name = treatment_name,
    cluster_name = cluster_name,
    terminal_event = terminal_event,
    event_weights = event_weights,
    treatment_levels = treatment_levels,
    treatment_probability = treatment_probability,
    times = times,
    estimand = estimand,
    control = control
  )
  fitted <- .fit_prepared_drcrtwa(
    prep = prep,
    times = prep$times,
    conf.level = conf.level,
    keep_influence = control$keep_influence,
    keep_nuisance = control$keep_nuisance
  )

  unit_arm <- prep$arm[match(seq_len(prep$n_units) - 1L, prep$unit)]
  analysis <- list(
    design = prep$design,
    n_long_rows = prep$n_long_rows,
    n_participants = prep$n_subjects,
    n_clusters = if (identical(prep$design, "CRT")) prep$n_units else NA_integer_,
    independent_units = prep$n_units,
    units_control = sum(unit_arm == 0L),
    units_treatment = sum(unit_arm == 1L),
    followup_range = range(prep$final_time),
    terminal_code = prep$terminal_code,
    recurrent_codes = prep$recurrent_codes,
    event_weights = prep$event_weights,
    recurrent_event_counts = prep$event_counts,
    terminal_events = prep$terminal_events,
    censored_participants = prep$censored_participants,
    estimands = prep$estimands,
    treatment_labels = prep$treatment_labels,
    randomization_probability_source = prep$probability_source,
    probability_treatment_range = range(prep$prob1)
  )

  object <- list(
    call = call,
    formula = formula,
    censoring_formula = censoring_formula,
    terminal_formula = terminal_formula,
    analysis = analysis,
    estimates = fitted$estimates,
    arm_estimates = fitted$arm_estimates,
    diagnostics = fitted$diagnostics,
    nuisance = fitted$nuisance,
    influence = fitted$influence,
    times = fitted$times,
    n_grid = fitted$n_grid,
    conf.level = conf.level,
    control = control,
    internal = if (isTRUE(control$retain_data)) prep else NULL
  )
  class(object) <- "DRCRTwa"
  object
}
