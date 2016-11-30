#' Plot method for cstab objects
#'
#' \code{plot.cstab} plots \emph{instability} path.
#'
#' @param x a cstab object (output of functions \code{cStability}).
#' @param \dots additional arguments passed to print.
#'
#' @author
#'    Jonas M. B. Haslbeck <jonas.haslbeck@gmail.com>
#'    Dirk U. Wulff <dirk.wulff@gmail.com>
#'
#' @method plot cstab
#'
#' @export
plot.cstab <- function(x, ...) {
  if(x$call$norm==FALSE) {
    graphics::plot(x$call$kseq,x$instab_path,xlab='k',ylab='Instability',las=1,type='b',...)
    }
  if(x$call$norm==TRUE) {
    graphics::plot(x$call$kseq,x$instab_path_norm,xlab='k',ylab='Instability',las=1,type='b',...)
    }
}
