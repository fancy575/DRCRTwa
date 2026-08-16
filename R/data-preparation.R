.prepare_event_weights <- function(event_weights, recurrent_codes) {
  k <- length(recurrent_codes)
  if (!k) {
    .drcrtwa_stop("At least one recurrent-event type must be observed.")
  }
  if (is.null(event_weights)) {
    out <- rep(1, k)
    names(out) <- as.character(recurrent_codes)
    return(out)
  }
  if (!is.numeric(event_weights) || any(!is.finite(event_weights))) {
    .drcrtwa_stop("`event_weights` must be a finite numeric vector.")
  }
  if (any(event_weights < 0)) {
    .drcrtwa_stop("Event weights must be nonnegative.")
  }
  if (!any(event_weights > 0)) {
    .drcrtwa_stop("At least one recurrent-event weight must be positive.")
  }
  code_names <- as.character(recurrent_codes)
  weight_names <- names(event_weights)
  unnamed <- is.null(weight_names) ||
    all(is.na(weight_names) | !nzchar(weight_names))
  if (unnamed) {
    if (length(event_weights) != k) {
      .drcrtwa_stop("An unnamed `event_weights` vector must have length ", k,
                    ", one value for each recurrent-event type.")
    }
    out <- as.numeric(event_weights)
    names(out) <- code_names
    return(out)
  }
  if (anyNA(weight_names) || any(!nzchar(weight_names)) ||
      anyDuplicated(weight_names)) {
    .drcrtwa_stop("Named event weights must have unique, nonempty names.")
  }
  missing <- setdiff(code_names, names(event_weights))
  extra <- setdiff(names(event_weights), code_names)
  if (length(missing) || length(extra)) {
    msg <- c(
      if (length(missing)) paste0("missing codes: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("unknown codes: ", paste(extra, collapse = ", "))
    )
    .drcrtwa_stop("The names of `event_weights` must exactly match the observed ",
                  "recurrent-event codes (", paste(msg, collapse = "; "), ").")
  }
  out <- as.numeric(event_weights[code_names])
  names(out) <- code_names
  out
}

.map_treatment <- function(x, treatment_levels = NULL) {
  if (anyNA(x)) .drcrtwa_stop("Treatment contains missing values.")
  x_chr <- as.character(x)
  observed <- unique(x_chr)
  if (is.null(treatment_levels)) {
    if (is.numeric(x) && setequal(sort(unique(as.numeric(x))), c(0, 1))) {
      lev <- c("0", "1")
    } else if (is.factor(x)) {
      lev <- levels(droplevels(x))
    } else {
      lev <- sort(observed)
    }
  } else {
    if (length(treatment_levels) != 2L || anyNA(treatment_levels) ||
        anyDuplicated(as.character(treatment_levels))) {
      .drcrtwa_stop("`treatment_levels` must contain exactly two distinct values, ",
                    "ordered as control then treatment.")
    }
    lev <- as.character(treatment_levels)
  }
  if (length(lev) != 2L || !setequal(observed, lev)) {
    .drcrtwa_stop("Treatment must have exactly two observed values. Supply ",
                  "`treatment_levels = c(control, treatment)` when their order ",
                  "is not unambiguous.")
  }
  arm <- match(x_chr, lev) - 1L
  list(
    arm = as.integer(arm),
    levels = lev,
    labels = c(control = lev[[1L]], treatment = lev[[2L]])
  )
}

.prepare_probabilities <- function(treatment_probability, data, subject_data,
                                   first_index, subject_key, unit, arm,
                                   design, treatment_labels) {
  n <- nrow(subject_data)
  n_units <- max(unit) + 1L
  if (is.null(treatment_probability)) {
    unit_arm <- arm[match(seq_len(n_units) - 1L, unit)]
    p1 <- rep(mean(unit_arm == 1L), n)
    source <- "empirical randomization fraction"
  } else {
    z <- treatment_probability
    if (is.character(z) && length(z) == 1L && z %in% names(data)) {
      z <- data[[z]]
    }
    if (!is.numeric(z)) {
      .drcrtwa_stop("`treatment_probability` must be NULL, a numeric value/vector, ",
                    "or the name of a numeric data column containing Pr(A=1).")
    }
    if (length(z) == 1L) {
      p1 <- rep(as.numeric(z), n)
    } else if (length(z) == 2L && length(z) != n) {
      if (!is.null(names(z))) {
        zn <- names(z)
        trt_candidates <- c("1", "treatment", treatment_labels[["treatment"]])
        hit <- match(trt_candidates, zn, nomatch = 0L)
        hit <- hit[hit > 0L]
        if (length(hit)) {
          p1 <- rep(as.numeric(z[[hit[[1L]]]]), n)
        } else if (abs(sum(z) - 1) < 1e-8) {
          p1 <- rep(as.numeric(z[[2L]]), n)
        } else {
          .drcrtwa_stop("Could not identify Pr(A=1) from the named ",
                        "`treatment_probability` vector.")
        }
      } else if (abs(sum(z) - 1) < 1e-8) {
        p1 <- rep(as.numeric(z[[2L]]), n)
      } else {
        .drcrtwa_stop("A length-two `treatment_probability` vector must contain ",
                      "control and treatment probabilities summing to one.")
      }
    } else if (length(z) == nrow(data)) {
      ok <- .is_constant_within(z, subject_key)
      if (any(!ok)) {
        .drcrtwa_stop("Row-level treatment probabilities must be constant within ",
                      "participant.")
      }
      p1 <- as.numeric(z[first_index])
    } else if (length(z) == n) {
      p1 <- as.numeric(z)
    } else if (length(z) == n_units) {
      p1 <- as.numeric(z[unit + 1L])
    } else {
      .drcrtwa_stop("`treatment_probability` has an unsupported length. Use one ",
                    "value, one value per participant, one per independent unit, ",
                    "or one per long-format row.")
    }
    source <- "user supplied"
  }
  if (any(!is.finite(p1)) || any(p1 <= 0 | p1 >= 1)) {
    .drcrtwa_stop("All treatment probabilities must lie strictly between zero and one.")
  }
  if (identical(design, "CRT")) {
    ok <- .is_constant_within(p1, unit)
    if (any(!ok)) {
      .drcrtwa_stop("Treatment probabilities must be constant within cluster for a CRT.")
    }
  }
  list(prob0 = 1 - p1, prob1 = p1, source = source)
}

.prepare_drcrtwa_data <- function(formula, censoring_formula, terminal_formula,
                                  data, id_name, treatment_name,
                                  cluster_name = NULL, terminal_event = NULL,
                                  event_weights = NULL,
                                  treatment_levels = NULL,
                                  treatment_probability = NULL,
                                  times = NULL, estimand = NULL,
                                  control = DRCRTwa_control()) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (!nrow(data)) .drcrtwa_stop("`data` has no rows.")

  response <- .extract_surv_response(formula)
  required <- unique(c(id_name, treatment_name, cluster_name,
                       response$start, response$time, response$status))
  missing_required <- setdiff(required, names(data))
  if (length(missing_required)) {
    .drcrtwa_stop("Required columns not found in `data`: ",
                  paste(missing_required, collapse = ", "), ".")
  }
  if (anyNA(data[required])) {
    .drcrtwa_stop("ID, treatment, cluster, time, start, and status fields may not ",
                  "contain missing values.")
  }

  time <- data[[response$time]]
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    .drcrtwa_stop("The follow-up time in the Surv response must be finite, numeric, ",
                  "and nonnegative.")
  }
  if (!is.null(response$start)) {
    start <- data[[response$start]]
    if (!is.numeric(start) || any(!is.finite(start)) || any(start < 0) ||
        any(start >= time)) {
      .drcrtwa_stop("For `Surv(start, stop, status)`, all rows must satisfy ",
                    "0 <= start < stop.")
    }
  }
  status_raw <- data[[response$status]]
  if (!is.numeric(status_raw) || any(!is.finite(status_raw)) ||
      any(status_raw < 0) || any(abs(status_raw - round(status_raw)) > 1e-8)) {
    .drcrtwa_stop("Status codes must be finite nonnegative integers.")
  }
  status <- as.integer(round(status_raw))
  positive_codes <- sort(unique(status[status > 0L]))
  if (!length(positive_codes)) {
    .drcrtwa_stop("No positive event code is present in `data`.")
  }
  if (is.null(terminal_event)) {
    terminal_code <- max(positive_codes)
  } else {
    if (!is.numeric(terminal_event) || length(terminal_event) != 1L ||
        !is.finite(terminal_event) || terminal_event <= 0 ||
        abs(terminal_event - round(terminal_event)) > 1e-8) {
      .drcrtwa_stop("`terminal_event` must be one positive integer status code.")
    }
    terminal_code <- as.integer(round(terminal_event))
  }
  recurrent_codes <- setdiff(positive_codes, terminal_code)
  weights <- .prepare_event_weights(event_weights, recurrent_codes)

  design <- if (is.null(cluster_name)) "IRT" else "CRT"
  if (anyNA(data[[id_name]]) ||
      (!is.null(cluster_name) && anyNA(data[[cluster_name]]))) {
    .drcrtwa_stop("Participant and cluster identifiers may not be missing.")
  }
  if (is.null(cluster_name)) {
    subject_key <- as.character(data[[id_name]])
  } else {
    subject_key <- paste(as.character(data[[cluster_name]]),
                         as.character(data[[id_name]]), sep = "\r")
  }

  reserved <- unique(c(id_name, treatment_name, cluster_name,
                       response$start, response$time, response$status))
  formula_vars <- unique(unlist(lapply(
    list(formula, censoring_formula, terminal_formula), .rhs_variables),
    use.names = FALSE))
  forbidden <- intersect(formula_vars, reserved)
  if (length(forbidden)) {
    .drcrtwa_stop("ID, treatment, cluster, follow-up-time, and status variables ",
                  "cannot be nuisance covariates because the models are fitted ",
                  "arm specifically. Remove: ", paste(forbidden, collapse = ", "), ".")
  }
  .validate_baseline_variables(
    data, subject_key,
    list(formula, censoring_formula, terminal_formula),
    reserved = character()
  )

  # Every participant must end with exactly one closing row: either status 0
  # (censoring/administrative end) or the terminal-event code.
  key_levels <- unique(subject_key)
  row_by_subject <- split(seq_len(nrow(data)), factor(subject_key,
                                                     levels = key_levels))
  first_index <- integer(length(key_levels))
  closing_index <- integer(length(key_levels))
  for (i in seq_along(row_by_subject)) {
    ind <- row_by_subject[[i]]
    first_index[[i]] <- ind[[which.min(time[ind])]]
    close <- ind[status[ind] %in% c(0L, terminal_code)]
    if (length(close) != 1L) {
      .drcrtwa_stop("Each participant must have exactly one closing row with ",
                    "status 0 or terminal status ", terminal_code,
                    ". The condition fails for participant key `",
                    key_levels[[i]], "`.")
    }
    closing_index[[i]] <- close
    max_time <- max(time[ind])
    if (abs(time[close] - max_time) > 1e-10 * max(1, abs(max_time))) {
      .drcrtwa_stop("The closing row must occur at the participant's largest ",
                    "follow-up time. The condition fails for participant key `",
                    key_levels[[i]], "`.")
    }
    if (sum(status[ind] == terminal_code) > 1L) {
      .drcrtwa_stop("The terminal event may occur at most once per participant.")
    }
  }

  # Treatment and all baseline model variables must be constant within subject.
  if (any(!.is_constant_within(data[[treatment_name]], subject_key))) {
    .drcrtwa_stop("Treatment must be constant within participant.")
  }
  subject_data <- data[first_index, , drop = FALSE]
  subject_data$.DRCRTwa_subject_key <- key_levels
  n_subjects <- nrow(subject_data)

  treatment_map <- .map_treatment(subject_data[[treatment_name]],
                                  treatment_levels = treatment_levels)
  arm <- treatment_map$arm

  if (identical(design, "CRT")) {
    cluster_subject <- as.character(subject_data[[cluster_name]])
    cluster_levels <- unique(cluster_subject)
    unit <- match(cluster_subject, cluster_levels) - 1L
    cluster_size <- tabulate(unit + 1L, nbins = length(cluster_levels))[unit + 1L]
    unit_treatment <- split(arm, unit)
    if (any(vapply(unit_treatment, function(z) length(unique(z)) != 1L,
                   logical(1L)))) {
      .drcrtwa_stop("Treatment must be constant within cluster for a CRT.")
    }
  } else {
    cluster_levels <- NULL
    unit <- seq_len(n_subjects) - 1L
    cluster_size <- rep.int(1L, n_subjects)
  }

  if (length(unique(arm)) != 2L) {
    .drcrtwa_stop("Both treatment arms must be represented among independent units.")
  }

  final_time <- as.numeric(time[closing_index])
  final_status <- as.integer(status[closing_index])
  if (any(final_time <= 0)) {
    .drcrtwa_stop("Each participant must have positive final follow-up time.")
  }

  recurrent_times <- lapply(recurrent_codes, function(code) {
    rows <- which(status == code)
    if (length(rows)) {
      by_key <- split(as.numeric(time[rows]),
                      factor(subject_key[rows], levels = key_levels))
      lapply(by_key, function(z) sort(as.numeric(z)))
    } else {
      rep(list(numeric()), n_subjects)
    }
  })
  names(recurrent_times) <- as.character(recurrent_codes)

  X_recurrent <- .model_matrix_rhs(formula, subject_data,
                                   "recurrent-event outcome")
  X_censoring <- .model_matrix_rhs(censoring_formula, subject_data,
                                   "censoring")
  X_terminal <- .model_matrix_rhs(terminal_formula, subject_data,
                                  "terminal-event outcome")
  if (nrow(X_recurrent) != n_subjects || nrow(X_censoring) != n_subjects ||
      nrow(X_terminal) != n_subjects) {
    .drcrtwa_stop("Internal model-matrix dimension check failed.")
  }

  probability <- .prepare_probabilities(
    treatment_probability = treatment_probability,
    data = data,
    subject_data = subject_data,
    first_index = first_index,
    subject_key = subject_key,
    unit = unit,
    arm = arm,
    design = design,
    treatment_labels = treatment_map$labels
  )

  if (is.null(estimand)) {
    estimands <- if (identical(design, "IRT")) {
      "individual-average"
    } else {
      c("individual-average", "cluster-average")
    }
  } else {
    estimand <- unique(as.character(estimand))
    if (!length(estimand) || anyNA(estimand) || any(!nzchar(estimand))) {
      .drcrtwa_stop("`estimand` must contain at least one nonempty estimand name.")
    }
    if ("both" %in% estimand) {
      estimand <- c("individual-average", "cluster-average")
    }
    allowed <- c("individual-average", "cluster-average")
    bad <- setdiff(estimand, allowed)
    if (length(bad)) {
      .drcrtwa_stop("Unknown estimand specification: ", paste(bad, collapse = ", "), ".")
    }
    if (identical(design, "IRT") && "cluster-average" %in% estimand) {
      .drcrtwa_stop("An IRT has only the individual-average estimand.")
    }
    estimands <- estimand
  }
  target_weights <- lapply(estimands, function(x) {
    if (identical(x, "cluster-average")) 1 / cluster_size else rep(1, n_subjects)
  })
  names(target_weights) <- estimands

  max_followup <- max(final_time)
  if (is.null(times)) {
    event_times <- sort(unique(as.numeric(time[status > 0L & time > 0])))
    analysis_times <- sort(unique(c(event_times, max_followup)))
  } else {
    if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
        any(times <= 0)) {
      .drcrtwa_stop("`times` must contain at least one positive finite analysis time.")
    }
    analysis_times <- sort(unique(as.numeric(times)))
  }
  if (any(analysis_times > max_followup + 1e-10 * max(1, max_followup))) {
    .drcrtwa_stop("Analysis times cannot exceed the largest observed follow-up time ",
                  "(", format(max_followup), ").")
  }
  list(
    formula = formula,
    censoring_formula = censoring_formula,
    terminal_formula = terminal_formula,
    design = design,
    response = response,
    id_name = id_name,
    treatment_name = treatment_name,
    cluster_name = cluster_name,
    subject_data = subject_data,
    subject_key = key_levels,
    arm = arm,
    treatment_labels = treatment_map$labels,
    unit = as.integer(unit),
    cluster_size = as.integer(cluster_size),
    cluster_levels = cluster_levels,
    n_subjects = n_subjects,
    n_units = max(unit) + 1L,
    final_time = final_time,
    final_status = final_status,
    terminal_code = terminal_code,
    recurrent_codes = as.integer(recurrent_codes),
    recurrent_times = recurrent_times,
    event_weights = weights,
    X_censoring = X_censoring,
    X_terminal = X_terminal,
    X_recurrent = X_recurrent,
    prob0 = probability$prob0,
    prob1 = probability$prob1,
    probability_source = probability$source,
    estimands = estimands,
    target_weights = target_weights,
    times = analysis_times,
    max_followup = max_followup,
    control = control,
    n_long_rows = nrow(data),
    event_counts = stats::setNames(
      vapply(recurrent_codes, function(code) sum(status == code), integer(1L)),
      as.character(recurrent_codes)
    ),
    terminal_events = sum(final_status == terminal_code),
    censored_participants = sum(final_status == 0L)
  )
}
