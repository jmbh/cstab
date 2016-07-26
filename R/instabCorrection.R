

#require(Rcpp)
#sourceCpp('/Users/jmb/Dropbox/MyData/_PhD/__software/mta/src//instabCorrection.cpp')

instabLookup = function(x,y,root = 200, lookupSize = 10000){
  if(!exists('lkup',envir = .GlobalEnv)) assign('lkup',lookup(lookupSize,root),envir=.GlobalEnv)
  a = stabExp(x,get('lkup',.GlobalEnv))
  b = stabExp(y,get('lkup',.GlobalEnv))
  return(a * (1-b) + (1-a) * b)
  }




