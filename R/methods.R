.evaluate_drcrtwa_times <- function(object, times, keep_all = FALSE) {
  if (is.null(times)) {
    return(list(
      estimates = object$estimates,
      arm_estimates = object$arm_estimates,
      times = object$times
    ))
  }
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      any(times <= 0)) {
    .drcrtwa_stop("`times` must contain at least one positive finite value.")
  }
  requested <- unique(as.numeric(times))
  if (any(requested > object$analysis$followup_range[[2L]] +
          1e-10 * max(1, object$analysis$followup_range[[2L]]))) {
    .drcrtwa_stop("Requested times cannot exceed the largest observed follow-up time.")
  }
  match_existing <- .match_times(object$times, requested)
  missing <- requested[is.na(match_existing)]
  extra <- NULL
  if (length(missing)) {
    if (is.null(object$internal)) {
      .drcrtwa_stop("The fitted object does not retain processed data. Refit with ",
                    "`DRCRTwa_control(retain_data = TRUE)` to evaluate new times.")
    }
    extra <- .fit_prepared_drcrtwa(
      prep = object$internal,
      times = missing,
      conf.level = object$conf.level,
      keep_influence = FALSE,
      keep_nuisance = FALSE
    )
  }
  estimates <- object$estimates
  arm_estimates <- object$arm_estimates
  available <- object$times
  if (!is.null(extra)) {
    estimates <- rbind(estimates, extra$estimates)
    arm_estimates <- rbind(arm_estimates, extra$arm_estimates)
    available <- c(available, extra$times)
  }
  select_time_rows <- function(x) {
    pieces <- lapply(requested, function(z) {
      hit <- abs(x$time - z) <= 1e-8 * pmax(1, abs(z))
      x[hit, , drop = FALSE]
    })
    .bind_rows_base(pieces)
  }
  list(
    estimates = if (keep_all) estimates else select_time_rows(estimates),
    arm_estimates = if (keep_all) arm_estimates else select_time_rows(arm_estimates),
    times = requested
  )
}

#' Extract a while-alive trajectory
#'
#' @param object A fitted `DRCRTwa` object.
#' @param times Optional analysis times. New times are evaluated from the
#'   retained subject histories when necessary.
#' @param type `"arms"` for treatment-specific while-alive rates or
#'   `"contrast"` for the treatment-minus-control contrast.
#' @param estimand Optional estimand names.
#' @param method One or more of `"dr"`, `"ipcw"`, and `"or"`.
#'
#' @return A data frame containing estimates, standard errors, and pointwise
#'   confidence intervals.
#' @export
wa_trajectory <- function(object, times = NULL,
                          type = c("arms", "contrast"),
                          estimand = NULL,
                          method = "dr") {
  if (!inherits(object, "DRCRTwa")) {
    .drcrtwa_stop("`object` must inherit from class `DRCRTwa`.")
  }
  type <- match.arg(type)
  method <- match.arg(method, c("dr", "ipcw", "or"), several.ok = TRUE)
  evaluated <- .evaluate_drcrtwa_times(object, times)
  out <- if (identical(type, "arms")) {
    evaluated$arm_estimates
  } else {
    evaluated$estimates
  }
  if (is.null(estimand)) estimand <- object$analysis$estimands
  bad <- setdiff(estimand, object$analysis$estimands)
  if (length(bad)) {
    .drcrtwa_stop("Unknown estimand: ", paste(bad, collapse = ", "), ".")
  }
  out <- out[out$estimand %in% estimand & out$method %in% method, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Print, summarize, and plot a DRCRTwa fit
#'
#' Methods for compact printing, time-specific summaries, and trajectory plots
#' from a fitted `DRCRTwa` object.
#'
#' @param x,object A fitted `DRCRTwa` object or, for
#'   `print.summary.DRCRTwa()`, an object returned by `summary()`.
#' @param digits Number of digits used when printing.
#' @param times Optional times to summarize or plot. The summary default is the
#'   largest fitted time; the plot default uses the full fitted trajectory.
#' @param estimand Optional estimand names.
#' @param method One or more methods for `summary()`, or one method for
#'   `plot()`: `"dr"`, `"ipcw"`, or `"or"`.
#' @param type `"arms"` for treatment-specific trajectories or `"contrast"`
#'   for treatment-minus-control trajectories.
#' @param conf.int Logical; include pointwise confidence ribbons.
#' @param ... For plotting, additional arguments passed to [ggplot2::theme()];
#'   otherwise unused.
#'
#' @return `print()` returns its input invisibly. `summary()` returns an object
#'   of class `summary.DRCRTwa`. `plot()` and `autoplot()` return a `ggplot`
#'   object.
#'
#' @name DRCRTwa-methods
NULL

#' @rdname DRCRTwa-methods
#' @export
print.DRCRTwa <- function(x, digits = 3L, ...) {
  cat("Doubly robust while-alive fit (DRCRTwa)\n")
  cat("Design: ", x$analysis$design, "\n", sep = "")
  if (identical(x$analysis$design, "CRT")) {
    cat("Participants: ", x$analysis$n_participants,
        "; clusters: ", x$analysis$n_clusters, "\n", sep = "")
  } else {
    cat("Participants: ", x$analysis$n_participants, "\n", sep = "")
  }
  cat("Treatment: ", x$analysis$treatment_labels[["treatment"]],
      " versus ", x$analysis$treatment_labels[["control"]], "\n", sep = "")
  cat("Recurrent-event weights: ",
      paste0(names(x$analysis$event_weights), "=",
             format(x$analysis$event_weights, trim = TRUE), collapse = ", "),
      "; terminal code: ", x$analysis$terminal_code, "\n", sep = "")
  cat("Fitted times: ", length(x$times), " (",
      format(min(x$times), digits = digits), " to ",
      format(max(x$times), digits = digits), ")\n", sep = "")

  endpoint <- x$estimates[
    x$estimates$method == "dr" &
      abs(x$estimates$time - max(x$times)) <= 1e-10 * max(1, max(x$times)),
    c("time", "estimand", "estimate", "std_error", "lower", "upper"),
    drop = FALSE
  ]
  names(endpoint) <- c("Time", "Estimand", "Estimate", "SE", "Lower", "Upper")
  for (j in c("Time", "Estimate", "SE", "Lower", "Upper"))
    endpoint[[j]] <- round(endpoint[[j]], digits)
  cat("\nDR treatment contrast at the largest fitted time:\n")
  print(endpoint, row.names = FALSE)
  invisible(x)
}

#' @rdname DRCRTwa-methods
#' @export
summary.DRCRTwa <- function(object, times = NULL, estimand = NULL,
                            method = "dr", ...) {
  if (is.null(times)) times <- max(object$times)
  method <- match.arg(method, c("dr", "ipcw", "or"), several.ok = TRUE)
  if (is.null(estimand)) estimand <- object$analysis$estimands
  bad <- setdiff(estimand, object$analysis$estimands)
  if (length(bad)) {
    .drcrtwa_stop("Unknown estimand: ", paste(bad, collapse = ", "), ".")
  }
  evaluated <- .evaluate_drcrtwa_times(object, times)
  estimates <- evaluated$estimates
  arms <- evaluated$arm_estimates
  estimates <- estimates[
    estimates$estimand %in% estimand & estimates$method %in% method,
    , drop = FALSE
  ]
  arms <- arms[arms$estimand %in% estimand & arms$method %in% method,
               , drop = FALSE]
  out <- list(
    call = object$call,
    analysis = object$analysis,
    formula = object$formula,
    censoring_formula = object$censoring_formula,
    terminal_formula = object$terminal_formula,
    times = evaluated$times,
    methods = method,
    estimands = estimand,
    estimates = estimates,
    arm_estimates = arms,
    conf.level = object$conf.level
  )
  class(out) <- "summary.DRCRTwa"
  out
}

#' @rdname DRCRTwa-methods
#' @export
print.summary.DRCRTwa <- function(x, digits = 3L, ...) {
  a <- x$analysis
  cat("DRCRTwa analysis summary\n")
  cat("Design: ", a$design, "\n", sep = "")
  if (identical(a$design, "CRT")) {
    cat("Participants: ", a$n_participants,
        "; randomized clusters: ", a$n_clusters,
        " (control ", a$units_control,
        ", treatment ", a$units_treatment, ")\n", sep = "")
  } else {
    cat("Participants: ", a$n_participants,
        " (control ", a$units_control,
        ", treatment ", a$units_treatment, ")\n", sep = "")
  }
  cat("Long-format rows: ", a$n_long_rows,
      "; terminal events: ", a$terminal_events,
      "; censored/administratively closed: ", a$censored_participants,
      "\n", sep = "")
  cat("Recurrent event types: ",
      paste0(names(a$recurrent_event_counts), " (n=",
             a$recurrent_event_counts, ", w=",
             a$event_weights[names(a$recurrent_event_counts)], ")",
             collapse = "; "), "\n", sep = "")
  cat("Treatment probability: ", a$randomization_probability_source,
      "; Pr(A=1) range ",
      paste(round(a$probability_treatment_range, digits), collapse = " to "),
      "\n", sep = "")
  cat("Pointwise confidence level: ",
      format(100 * x$conf.level, trim = TRUE), "%\n", sep = "")

  arm <- x$arm_estimates[, c(
    "time", "estimand", "arm_label", "method", "estimate", "std_error",
    "lower", "upper", "weighted_burden", "rmst"
  ), drop = FALSE]
  names(arm) <- c("Time", "Estimand", "Treatment", "Method", "Rate", "SE",
                  "Lower", "Upper", "Burden", "RMST")
  numeric_arm <- c("Time", "Rate", "SE", "Lower", "Upper", "Burden", "RMST")
  arm[numeric_arm] <- lapply(arm[numeric_arm], round, digits = digits)
  cat("\nArm-specific exposure-weighted while-alive rates:\n")
  print(arm, row.names = FALSE)

  contrast <- x$estimates[, c(
    "time", "estimand", "method", "estimate", "std_error", "lower", "upper",
    "p_value"
  ), drop = FALSE]
  names(contrast) <- c("Time", "Estimand", "Method", "Difference", "SE",
                       "Lower", "Upper", "P value")
  numeric_contrast <- c("Time", "Difference", "SE", "Lower", "Upper", "P value")
  contrast[numeric_contrast] <- lapply(
    contrast[numeric_contrast], round, digits = digits)
  cat("\nTreatment-minus-control contrasts:\n")
  print(contrast, row.names = FALSE)
  invisible(x)
}

#' @rdname DRCRTwa-methods
#' @export
plot.DRCRTwa <- function(x, type = c("arms", "contrast"), times = NULL,
                         estimand = NULL, method = "dr",
                         conf.int = TRUE, ...) {
  type <- match.arg(type)
  method <- match.arg(method, c("dr", "ipcw", "or"))
  dat <- wa_trajectory(x, times = times, type = type,
                       estimand = estimand, method = method)
  if (!nrow(dat)) .drcrtwa_stop("No trajectory rows match the plotting request.")

  if (identical(type, "arms")) {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = .data$time, y = .data$estimate, color = .data$arm_label,
        fill = .data$arm_label, linetype = .data$arm_label,
        group = .data$arm_label
      )
    )
    if (isTRUE(conf.int)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.16, color = NA
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::labs(
        x = "Time",
        y = "Exposure-weighted while-alive rate",
        color = "Treatment", fill = "Treatment", linetype = "Treatment"
      )
  } else {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = .data$time, y = .data$estimate,
                   group = .data$estimand)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5)
    if (isTRUE(conf.int)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.16
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::labs(
        x = "Time",
        y = "Treatment-minus-control while-alive rate difference"
      )
  }
  if (length(unique(dat$estimand)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$estimand),
                              scales = "free_y")
  }
  p + ggplot2::theme_classic() + ggplot2::theme(...)
}

#' @rdname DRCRTwa-methods
#' @export
autoplot.DRCRTwa <- function(object, ...) {
  plot.DRCRTwa(object, ...)
}
