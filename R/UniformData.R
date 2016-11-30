
# Creates uniformly distributed data of same dimensionality as input data

UniformData <- function(data) {

  # get dimensions of data
  n         <- nrow(data)
  unifdims  <- t(apply(data, 2, range))

  # sample data and combine to data matrix
  data_synt  <- matrix(NA,nrow=nrow(data),ncol=ncol(data))
  for(i in 1:ncol(data)) data_synt[,i] <- stats::runif(n, unifdims[i,1], unifdims[i,2])
  return(data_synt)
  }
