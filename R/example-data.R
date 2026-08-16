.read_drcrtwa_example <- function(file) {
  path <- system.file("extdata", file, package = "DRCRTwa")
  if (!nzchar(path)) {
    .drcrtwa_stop("Installed example data file was not found: ", file, ".")
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Simulated IRT data in long event-history format
#'
#' Returns a reproducible example with one row per recurrent event and one final
#' row per participant. Status 0 denotes censoring or administrative closure,
#' status 1 and 2 denote two recurrent-event types, and status 3 denotes death.
#'
#' @return A data frame with columns `id`, `treatment`, `time`, `status`, `z1`,
#'   `z2`, and `z3`.
#' @export
wa_irt_data <- function() {
  .read_drcrtwa_example("wa_irt.csv")
}

#' Simulated CRT data in long event-history format
#'
#' Returns a reproducible cluster-randomized example. Participant identifiers
#' intentionally repeat across clusters, so `cluster` and `id` jointly identify
#' a participant. Status 0 denotes censoring or administrative closure, status 1
#' and 2 denote recurrent-event types, and status 3 denotes death.
#'
#' @return A data frame with columns `cluster`, `id`, `treatment`, `time`,
#'   `status`, `z1`, `z2`, and `z3`.
#' @export
wa_crt_data <- function() {
  .read_drcrtwa_example("wa_crt.csv")
}
