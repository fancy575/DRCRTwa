#' Control options for DRCRTwa
#'
#' @param grid_per_unit Number of equally spaced numerical-integration bins per
#'   unit of follow-up time.
#' @param min_grid Minimum number of integration bins at each analysis time.
#' @param max_grid Maximum number of integration bins at each analysis time.
#' @param keep_influence Logical; retain independent-unit influence values.
#' @param keep_nuisance Logical; retain nuisance fits at the largest fitted time.
#' @param retain_data Logical; retain processed subject histories and design
#'   matrices so that `summary()` can evaluate additional times.
#'
#' @return A list of control settings.
#' @export
DRCRTwa_control <- function(grid_per_unit = 100,
                            min_grid = 50L,
                            max_grid = 5000L,
                            keep_influence = FALSE,
                            keep_nuisance = TRUE,
                            retain_data = TRUE) {
  if (!is.numeric(grid_per_unit) || length(grid_per_unit) != 1L ||
      !is.finite(grid_per_unit) || grid_per_unit <= 0) {
    .drcrtwa_stop("`grid_per_unit` must be one positive finite number.")
  }
  min_grid <- as.integer(min_grid)
  max_grid <- as.integer(max_grid)
  if (length(min_grid) != 1L || is.na(min_grid) || min_grid < 2L) {
    .drcrtwa_stop("`min_grid` must be an integer of at least 2.")
  }
  if (length(max_grid) != 1L || is.na(max_grid) || max_grid < min_grid) {
    .drcrtwa_stop("`max_grid` must be at least `min_grid`.")
  }
  list(
    grid_per_unit = as.numeric(grid_per_unit),
    min_grid = min_grid,
    max_grid = max_grid,
    keep_influence = isTRUE(keep_influence),
    keep_nuisance = isTRUE(keep_nuisance),
    retain_data = isTRUE(retain_data)
  )
}
