#' @useDynLib cstab
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  assign('lkup',lookup(n = 10000, root = 200),envir = parent.env(environment()))
  }
