.drcrtwa_stop <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

.bind_rows_base <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.resolve_column_name <- function(expr, data, argument, envir) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
  } else if (is.character(expr) && length(expr) == 1L) {
    nm <- expr
  } else {
    value <- try(eval(expr, envir = envir), silent = TRUE)
    if (inherits(value, "try-error") || !is.character(value) ||
        length(value) != 1L) {
      .drcrtwa_stop("`", argument,
                    "` must be an unquoted column name or a single character string.")
    }
    nm <- value
  }
  if (!nm %in% names(data)) {
    .drcrtwa_stop("Column `", nm, "` supplied through `", argument,
                  "` is not present in `data`.")
  }
  nm
}

.extract_surv_response <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) < 3L) {
    .drcrtwa_stop("`formula` must be a two-sided formula with a Surv response.")
  }
  lhs <- formula[[2L]]
  if (!is.call(lhs) || !grepl("Surv", paste(deparse(lhs[[1L]]), collapse = ""),
                              fixed = TRUE)) {
    .drcrtwa_stop("The left side of `formula` must be `Surv(time, status)` or ",
                  "`Surv(start, stop, status)`.")
  }
  vars <- all.vars(lhs)
  if (length(vars) == 2L) {
    list(start = NULL, time = vars[[1L]], status = vars[[2L]])
  } else if (length(vars) == 3L) {
    list(start = vars[[1L]], time = vars[[2L]], status = vars[[3L]])
  } else {
    .drcrtwa_stop("The Surv response must contain two or three simple variables.")
  }
}

.rhs_variables <- function(formula) {
  if (!inherits(formula, "formula")) {
    .drcrtwa_stop("All nuisance-model specifications must be formulas.")
  }
  if (length(formula) == 2L) {
    all.vars(formula[[2L]])
  } else {
    all.vars(formula[[3L]])
  }
}

.rhs_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    .drcrtwa_stop("All nuisance-model specifications must be formulas.")
  }
  if (length(formula) == 2L) formula else stats::delete.response(stats::terms(formula))
}

.model_matrix_rhs <- function(formula, data, label) {
  rhs <- .rhs_formula(formula)
  mf <- try(stats::model.frame(rhs, data = data, na.action = stats::na.fail,
                               drop.unused.levels = TRUE), silent = TRUE)
  if (inherits(mf, "try-error")) {
    .drcrtwa_stop("Could not construct the ", label,
                  " model frame. Check missing values and formula variables.\n",
                  as.character(mf))
  }
  mm <- stats::model.matrix(rhs, data = mf)
  if (ncol(mm) && "(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  storage.mode(mm) <- "double"
  mm
}

.is_constant_within <- function(x, group) {
  split_x <- split(x, group, drop = TRUE)
  vapply(split_x, function(z) {
    z <- z[!is.na(z)]
    length(unique(z)) <= 1L
  }, logical(1L))
}

.validate_baseline_variables <- function(data, subject_key, formulas,
                                         reserved = character()) {
  vars <- unique(unlist(lapply(formulas, .rhs_variables), use.names = FALSE))
  vars <- setdiff(vars, reserved)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    .drcrtwa_stop("Formula variables not found in `data`: ",
                  paste(missing_vars, collapse = ", "), ".")
  }
  bad <- character()
  for (v in vars) {
    ok <- .is_constant_within(data[[v]], subject_key)
    if (any(!ok)) bad <- c(bad, v)
  }
  if (length(bad)) {
    .drcrtwa_stop(
      "The current estimator uses baseline nuisance covariates. These variables ",
      "vary within participant: ", paste(unique(bad), collapse = ", "), "."
    )
  }
  invisible(vars)
}

.match_times <- function(available, requested, tolerance = 1e-8) {
  vapply(requested, function(z) {
    hit <- which(abs(available - z) <= tolerance * pmax(1, abs(z)))
    if (length(hit)) hit[[1L]] else NA_integer_
  }, integer(1L))
}

.format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}
