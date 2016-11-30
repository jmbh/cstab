
getMeasures <- function(data,
                        k,
                        method = 'kmeans',
                        linkage = 'complete',
                        kmIter = 5,
                        measures = c('wcd','sil','mse'))

{

  # ------ COMPUTE DISTANCE MATRIX

  # calc distance matric
  dmat <- fast_dist(data)


  # ------ APPLY CLUSTERING

  # clustering
  if(k==1) cl <- rep(1,nrow(data))
  if(k > 1) {
    # k-means clustering
    if(method=='kmeans') {
      check_k <- FALSE
      counter <- 0
      while(check_k == FALSE) {

        # terminate if necessary
        counter <- counter + 1
        if(counter>100) stop(paste0('k means solution with ', k, 'centers always converges to solution with empty cluster.'))

        # apply kmeans
        kms <- list() # re-start k means algorithm a couple of times
        for(km_i in 1:kmIter) {
          kms[[km_i]] <- stats::kmeans(data, centers = k, iter.max = 100) #save whole model
          }

        # extract
        wcd      <- unlist(lapply(kms, function(x) x$tot.withinss))
        km_model <- kms[[which.min(wcd)]] # pick k-means clustering with smallest WCD
        cl       <- km_model$cluster
        if(length(unique(cl)) == k) check_k <- TRUE
        }
      }

    # hierarhical clustering
    if(method=='hierarchical') {
      hcobj <- fastcluster::hclust(stats::as.dist(dmat), method = linkage)
      cl <- stats::cutree(hcobj, k)
      }
    }



  # ------ WITHIN CLUSTER DISSIMILARITY

  # calc within cluster dissimilarity
  WCD = NULL
  if('wcd' %in% measures){
    norm_diss <- c()
    for(i in 1:k) norm_diss[i] <- sum(dmat[cl==i, cl==i]) / (2 * sum(cl == i))
    WCD <- sum(norm_diss)
    }


  # ------ SILHOUETTE

  Sil <- NULL
  if('sil' %in% measures) if(k > 1) Sil <- mean(cluster::silhouette(cl, stats::as.dist(dmat))[,3]) else Sil <- 0



  # ------ CLUSTER CENTERS

  # get centers
  centers = NULL
  if('mse' %in% measures | 'centers' %in% measures){
    if(k == 1)  centers <- colMeans(data)
    if(k > 1) {
      if(method == 'kmeans') {
        centers <- km_model$centers
        } else {
        data_or <- as.data.frame(data)
        data_or$cl <- cl
        data_or <- data_or[order(data_or$cl),]
        centers <- matrix(NA, k, ncol(data))
        for(kk in 1:k) centers[kk,] <- colMeans(data_or[data_or$cl==kk, -(ncol(data)+1)])
        }
      }
    }


  # ------ MSE FOR JUMP

  MSE = NULL
  if('mse' %in% measures) {
    if(k == 1) SE <- apply(data, 1, function(inst){diffs <- inst - centers; return(t(diffs) %*% diffs)})
    if(k > 1) {
      SE <- numeric()
      for(i in 1:nrow(data)) {
        diffs <- data[i,] - centers[cl[i],]
        se    <- t(diffs) %*% diffs
        SE[i] <- se
        }
      }
    MSE <- sum(SE) / nrow(data)
    }


  # ------ RETURN

  outlist <- list('WCD' = WCD, 'Sil'=Sil, 'MSE'=MSE, 'Centers'=centers)
  return(outlist)
  }

