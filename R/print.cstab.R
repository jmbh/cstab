#' Print method for cstab objects
#'
#' \code{print.cstab} prints key variables of cstab objects.
#'
#' @param x a cstab object (output of functions \code{cStability}).
#' @param \dots additional arguments passed to print.
#'
#' @author
#'    Jonas M. B. Haslbeck <jonas.haslbeck@gmail.com>
#'    Dirk U. Wulff <dirk.wulff@gmail.com>
#'
#' @export

print.cstab <- function(x, ...) {
  if(x$call$norm==FALSE) {
    l <- list('kOptimal'=x$k_instab, 'instabilityPath'=x$instab_path)
    print(l, ...)
    }
  if(x$call$norm==TRUE) {
    l <- list('kOptimal'=x$k_instab_norm, 'instabilityPath'=x$instab_path_norm)
    print(l, ...)
    }
  }


