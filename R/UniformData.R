
# Creates uniformly distributed data of same dimensionality as input data

UniformData <- function(data) {

  # get dimensions of data
  n <- nrow(data)
  unifdims <- t(apply(data, 2, range))
  diffs <- abs(unifdims[,1] - unifdims[,2])
  unifdims2 <- cbind(0,diffs)
  dims <- ncol(data)

  # sample data and combine to data matrix
  unif_dims <- list()
  for(i in 1:dims) unif_dims[[i]] <- runif(n, unifdims2[i,1], unifdims2[i,2])
  data_synt <- do.call(cbind, unif_dims)

  return(data_synt)
}
