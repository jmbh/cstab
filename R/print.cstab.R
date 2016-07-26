
print.cstab <- function(x, ...) {

  if(x$call$norm==TRUE) {
    ls <- list('kOptimal'=x$kopt_instab_norm, 'InstabilityPath'=x$Instab_path)
    print(ls)
  }

  if(x$call$norm==FALSE) {
    ls <- list('kOptimal'=x$kopt_instab, 'InstabilityPath'=x$Instab_path_norm)
    print(ls)
  }

}
