

GetMeasures <- function(data,
                        k,
                        method='kmeans',
                        linkage = 'complete',
                        kmIter = 10,
                        type='data')

{

  # clustering
  if(k==1) {
    cl <- rep(1,nrow(data))
  } else {
    check_k <- FALSE
    counter <- 0

      if(method=='kmeans') {
        # k-means clustering
        while(check_k == FALSE) {

        l_km <- list() # re-start k means algorithm a couple of times
        for(km in 1:kmIter) {
          l_km[[km]] <- kcca(data, k=k, kccaFamily("kmeans")) #save whole model
        }

        WCD <- unlist(lapply(l_km, function(x) mean(x@clusinfo$av_dist)))
        km_model <- l_km[[which.min(WCD)]] # pick k-means clustering with smallest WCD

        cl <- km_model@cluster
        if(length(unique(cl))==k) check_k <- TRUE
        counter <- counter + 1
        if(counter>100) stop(paste0('k means solution with ', k, 'centers always converges to solution with empty cluster.'))
      }
    } else {
      # hierarhical clustering
      hcobj <- hclust(dist(data), method = linkage)
      cl <- cutree(hcobj, k)
    }
  }

  # calc distance matric
  dmat <- as.matrix(dist(data))

  # calc within cluster dissimilarity
  l_diss <- list()
  for(i in 1:k) l_diss[[i]] <- mean(dmat[cl==i, cl==i])

  tb <- table(cl)
  WCD <- sum(unlist(l_diss) * tb/sum(tb))

  if(k>1) {
    if(type=='data') {
      Sil <- mean(silhouette(cl, dmat)[,3])
    } else {
      Sil <- NULL
    }
  } else {
    Sil <- NULL
  }

  # calc MSE (for jump statistic)

  # get centers
  if(method == 'kmeans') {
    centers <- km_model@centers
  } else {
    data_or <- as.data.frame(data)
    data_or$cl <- cl
    data_or <- data_or[order(data_or$cl),]
    centers <- matrix(NA, k, ncol(data))
    for(kk in 1:k) centers[kk,] <- colMeans(data_or[data_or$cl==kk, -(ncol(data)+1)])
  }

  if(type=='data') {
    v_mse <- numeric(nrow(data))
    if(k>1) {
      for(i in 1:nrow(data)) {
        diffs <- (data[i,] - centers[km_model@cluster[i],])
        mse <- t(diffs) %*% diffs
        v_mse[i] <- mse
      }
      MSE <- sum(v_mse) / ncol(data)
    } else {
      SE <- apply(data, 1, function(inst){
        diffs <- inst - colMeans(data)
        return(t(diffs) %*% diffs)
      })
      MSE <- sum(SE) / ncol(data)
    }
  } else {
    MSE <- NULL
  }

  outlist <- list('WCD' = WCD, 'Sil'=Sil, 'MSE'=MSE, 'Centers'=centers)
  return(outlist)
}

