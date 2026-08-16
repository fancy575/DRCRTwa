#' DRCRTwa: Doubly Robust While-Alive Estimation for Randomized Trials
#'
#' The `DRCRTwa` package implements compiled doubly robust estimation of
#' exposure-weighted while-alive recurrent-event rates for individually
#' randomized and cluster-randomized trials in the presence of a terminal event
#' and covariate-dependent censoring. The main fitting function is
#' [DRCRT_WA()].
#'
#' @keywords internal
#' @useDynLib DRCRTwa, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
"_PACKAGE"
