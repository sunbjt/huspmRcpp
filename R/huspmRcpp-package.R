#' huspmRcpp: High-Utility Sequential Pattern Mining
#'
#' Mines high-utility sequential patterns from event data with strictly positive
#' utilities. See [huspm()] for the mining entry point and [prepare_huspm_data()]
#' for the input contract.
#'
#' The `importFrom(Rcpp, evalCpp)` below is load-bearing, not decorative.
#' `DESCRIPTION` declaring `Imports: Rcpp` does not by itself load Rcpp's
#' namespace or its DLL, so Rcpp's C-callables (`Rcpp_precious_remove`,
#' `enterRNGScope`, ...) stay unregistered and every `.Call` into this package
#' fails with "no function 'Rcpp' provided by package 'Rcpp_precious_remove'".
#' `Rcpp::compileAttributes()` does not emit this line either.
#'
#' @useDynLib huspmRcpp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
