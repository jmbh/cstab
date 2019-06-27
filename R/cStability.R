#' Selection of number of clusters via clustering instability
#'
#' Selection of number of clusters via \emph{model-based} or \emph{model-free},
#'   \emph{normalized} or \emph{unnormalized} clustering instability.
#'
#' @param data a n x p data matrix of type numeric.
#' @param kseq a vector with considered numbers clusters k > 1
#' @param nB an integer specifying the number of bootstrap comparisons.
#' @param norm logical specifying whether the instability path should be
#'    normalized. If TRUE, the instability path is normalized, accounting for a
#'    trivial decrease in instability due to a increasing k (see Haslbeck & Wulff,
#'    2016).
#' @param predict boolean specifying whether the model-based or the model-free
#'    variant should be used (see Haslbeck & Wulff, 2016).
#' @param method character string specifying the clustering algorithm. 'kmeans' for
#'    the k-means algorithm, 'hierarchical' for hierarchical clustering.
#' @param linkage character specifying the linkage criterion, in case
#'    \code{type='hierarchical'}. The available options are "single", "complete",
#'    "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median". See
#'    \link[stats]{hclust}.
#' @param kmIter integer specifying the the number of restarts of the k-means algorithm
#'    in order to avoid local minima.
#' @param pbar logical
#'
#' @return a list that contains the optimal k selected by the unnormalized and
#'    normalized instability method. It also includes a vector containing the averaged
#'    instability path (over bootstrap samples) and a matrix containing the instability
#'    path of each bootstrap sample for both the normalized and the unnormalized method.
#'
#' @references
#'    Ben-Hur, A., Elisseeff, A., & Guyon, I. (2001). A stability based method for
#'    discovering structure in clustered data. \emph{Pacific symposium on biocomputing,
#'    7}, 6-17.
#'
#'    Tibshirani, R., & Walther, G. (2005). Cluster validation by prediction strength.
#'    \emph{Journal of Computational and Graphical Statistics, 14}(3), 511-528.
#'
#' @author
#'    Dirk U. Wulff <dirk.wulff@gmail.com>
#'    Jonas M. B. Haslbeck <jonas.haslbeck@gmail.com>
#'
#' @examples
#' \dontrun{
#'   # Generate Data from Gaussian Mixture
#'   s <- .1
#'   n <- 50
#'   data <- rbind(cbind(rnorm(n, 0, s), rnorm(n, 0, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 0, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 0, s)))
#'   plot(data)
#'
#'   # Selection of Number of Clusters using Instability-based Measures
#'   stab_obj <- cStability(data, kseq=2:10)
#'   print(stab_obj)
#'   }
#'
#' @export

cStability_orig <- function(data, # n x p data matrix
                            kseq = 2:20, # sequence of ks tested
                            nB   = 10, # number of bootstrap comparisons
                            norm = TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                            predict = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                            method    = 'kmeans', # or 'hierarchical'
                            linkage = 'complete', # or average, or ...
                            kmIter  = 5,
                            pbar = TRUE) # number of reruns of k-means algorithm
{


  # ---------- Input Checks ----------

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')

  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')

  # On B
  if(round(nB)!=nB) stop('The number of bootstrap comparison has to be a positive integer value.')

  # On type
  if(method %in% c('kmeans', 'hierarchical'))


  # ---------- Create Containers ----------

  n_obj <- nrow(data)
  m_instab <- matrix(NA, nB, length(kseq)) # storage: no normalization
  m_instab_norm <- matrix(NA, nB, length(kseq)) # storage: normalization


  # ---------- Draw bootstrap samples ----------

  ind <- list()
  share <- list()
  for(b in 1:nB) {
    tmp_1 <- sample(1:n_obj, n_obj, replace=T)
    tmp_2 <- sample(1:n_obj, n_obj, replace=T)
    tmp_1 <- tmp_1[order(tmp_1)]
    tmp_2 <- tmp_2[order(tmp_2)]
    ind[[b]] <- list(tmp_1,tmp_2)
    intersct <- intersect(tmp_1,tmp_2)
    share[[b]] <- list(tmp_1 %in% intersct & !duplicated(tmp_1), tmp_2 %in% intersct & !duplicated(tmp_2)) # leave duplicates in ?
    }


  # ---------- Calculate distance matrix ----------

  distm = fast_dist(data) # maybe here?

  # ----- Loop over B comparisons -----

  if(pbar)  pb <- utils::txtProgressBar(min=0, max=nB, style = 2)

  for(b in 1:nB) {
    for(k in kseq) {

        # kmeans
        if(method == 'kmeans') {
          kms_1 <- list()
          kms_2 <- list()
          for(km_i in 1:kmIter) {
            kms_1[[km_i]] = stats::kmeans(data[ind[[b]][[1]],], centers = k)
            kms_2[[km_i]] = stats::kmeans(data[ind[[b]][[2]],], centers = k) # this has to be 2
            }
          km_1 = kms_1[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
          km_2 = kms_2[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
          if(predict == TRUE) {
            cl_1 = kmeans_predict(km_1, data = data)
            cl_2 = kmeans_predict(km_2, data = data)
            } else {
            cl_1 = km_1$cluster[share[[b]][[1]]]
            cl_2 = km_2$cluster[share[[b]][[2]]]
            }
          }

        if(method == 'hierarchical') {
          #t = proc.time()[3]
          hc_1 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[1]],ind[[b]][[1]]]), method = linkage)
          hc_2 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[2]],ind[[b]][[2]]]), method = linkage)
          cl_1 = stats::cutree(hc_1, k)[share[[b]][[1]]]
          cl_2 = stats::cutree(hc_2, k)[share[[b]][[2]]]
          #print(proc.time()[3] - t)
          }

      # check for equality of clusterings
      eq_1 = equal(cl_1)
      eq_2 = equal(cl_2)
      InStab <- mean(eq_1 != eq_2)

      # Normalize = FALSE
      m_instab[b, which(kseq==k)] <- InStab

      # Normalize = TRUE
      norm_val <- instabLookup(table(cl_1), table(cl_2))
      m_instab_norm[b, which(kseq==k)] <- InStab / norm_val

    } # end for k

    if(pbar) utils::setTxtProgressBar(pb, b)

  } # end for B

  # taking the mean
  m_instab_M      = colMeans(m_instab)
  m_instab_norm_M = colMeans(m_instab_norm)

  kopt_instab  = which(m_instab_M == min(m_instab_M))+(min(kseq)-1)
  kopt_instabN = which(m_instab_norm_M == min(m_instab_norm_M))+(min(kseq)-1)

  # replicate function call
  f_call <- list('kseq'=kseq,
                 'nB'=nB,
                 'norm'=norm,
                 'prediction'=predict,
                 'method'=method,
                 'linkage'=linkage,
                 'kmIter'=kmIter)

  outlist <- list("k_instab"=kopt_instab,
                  "k_instab_norm"=kopt_instabN,
                  "instab_path"=m_instab_M,
                  "instab_path_norm"=m_instab_norm_M,
                  "instab_path_matrix"=m_instab,
                  "instab_path_nrom_matrix"=m_instab_norm,
                  'call'=f_call)

  class(outlist) <- c('list', 'cstab', 'cStability')

  return(outlist)

} # EoF

#' Selection of number of clusters via clustering instability
#'
#' Selection of number of clusters via \emph{model-based} or \emph{model-free},
#'   \emph{normalized} or \emph{unnormalized} clustering instability.
#'
#' @param data a n x p data matrix of type numeric.
#' @param kseq a vector with considered numbers clusters k > 1
#' @param nB an integer specifying the number of bootstrap comparisons.
#' @param norm logical specifying whether the instability path should be
#'    normalized. If TRUE, the instability path is normalized, accounting for a
#'    trivial decrease in instability due to a increasing k (see Haslbeck & Wulff,
#'    2016).
#' @param predict boolean specifying whether the model-based or the model-free
#'    variant should be used (see Haslbeck & Wulff, 2016).
#' @param method character string specifying the clustering algorithm. 'kmeans' for
#'    the k-means algorithm, 'hierarchical' for hierarchical clustering.
#' @param linkage character specifying the linkage criterion, in case
#'    \code{type='hierarchical'}. The available options are "single", "complete",
#'    "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median". See
#'    \link[stats]{hclust}.
#' @param kmIter integer specifying the the number of restarts of the k-means algorithm
#'    in order to avoid local minima.
#' @param pbar logical
#'
#' @return a list that contains the optimal k selected by the unnormalized and
#'    normalized instability method. It also includes a vector containing the averaged
#'    instability path (over bootstrap samples) and a matrix containing the instability
#'    path of each bootstrap sample for both the normalized and the unnormalized method.
#'
#' @references
#'    Ben-Hur, A., Elisseeff, A., & Guyon, I. (2001). A stability based method for
#'    discovering structure in clustered data. \emph{Pacific symposium on biocomputing,
#'    7}, 6-17.
#'
#'    Tibshirani, R., & Walther, G. (2005). Cluster validation by prediction strength.
#'    \emph{Journal of Computational and Graphical Statistics, 14}(3), 511-528.
#'
#' @author
#'    Dirk U. Wulff <dirk.wulff@gmail.com>
#'    Jonas M. B. Haslbeck <jonas.haslbeck@gmail.com>
#'
#' @examples
#' \dontrun{
#'   # Generate Data from Gaussian Mixture
#'   s <- .1
#'   n <- 50
#'   data <- rbind(cbind(rnorm(n, 0, s), rnorm(n, 0, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 0, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 0, s)))
#'   plot(data)
#'
#'   # Selection of Number of Clusters using Instability-based Measures
#'   stab_obj <- cStability(data, kseq=2:10)
#'   print(stab_obj)
#'   }
#'
#' @export

cStability <- function(data, # n x p data matrix
                       kseq = 2:20, # sequence of ks tested
                       nB   = 10, # number of bootstrap comparisons
                       norm = TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                       predict = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                       method    = 'kmeans', # or 'hierarchical'
                       linkage = 'complete', # or average, or ...
                       kmIter  = 5,
                       pbar = TRUE) # number of reruns of k-means algorithm
{


  # ---------- Input Checks ----------

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')

  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')

  # On B
  if(round(nB)!=nB) stop('The number of bootstrap comparison has to be a positive integer value.')

  # On type
  if(method %in% c('kmeans', 'hierarchical'))


    # ---------- Create Containers ----------

  n_obj <- nrow(data)
  m_instab <- matrix(NA, nB, length(kseq)) # storage: no normalization
  m_instab_norm <- matrix(NA, nB, length(kseq)) # storage: normalization


  # ---------- Draw bootstrap samples ----------

  ind <- list()
  share <- list()
  for(b in 1:nB) {
    tmp_1 <- sample(1:n_obj, n_obj, replace=T)
    tmp_2 <- sample(1:n_obj, n_obj, replace=T)
    tmp_1 <- tmp_1[order(tmp_1)]
    tmp_2 <- tmp_2[order(tmp_2)]
    ind[[b]] <- list(tmp_1,tmp_2)
    intersct <- intersect(tmp_1,tmp_2)
    share[[b]] <- list(tmp_1 %in% intersct & !duplicated(tmp_1), tmp_2 %in% intersct & !duplicated(tmp_2)) # leave duplicates in ?
  }


  # ---------- Calculate distance matrix ----------

  distm = fast_dist(data) # maybe here?

  # ----- Loop over B comparisons -----

  if(pbar)  pb <- utils::txtProgressBar(min=0, max=nB, style = 2)

  for(b in 1:nB) {
    for(k in kseq) {

      # kmeans
      if(method == 'kmeans') {
        kms_1 <- list()
        kms_2 <- list()
        for(km_i in 1:kmIter) {
          kms_1[[km_i]] = stats::kmeans(data[ind[[b]][[1]],], centers = k)
          kms_2[[km_i]] = stats::kmeans(data[ind[[b]][[2]],], centers = k) # this has to be 2
        }
        km_1 = kms_1[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
        km_2 = kms_2[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
        if(predict == TRUE) {
          cl_1 = kmeans_predict(km_1, data = data)
          cl_2 = kmeans_predict(km_2, data = data)
        } else {
          cl_1 = km_1$cluster[share[[b]][[1]]]
          cl_2 = km_2$cluster[share[[b]][[2]]]
        }
      }

      if(method == 'hierarchical') {
        #t = proc.time()[3]
        hc_1 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[1]],ind[[b]][[1]]]), method = linkage)
        hc_2 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[2]],ind[[b]][[2]]]), method = linkage)
        cl_1 = stats::cutree(hc_1, k)[share[[b]][[1]]]
        cl_2 = stats::cutree(hc_2, k)[share[[b]][[2]]]
        #print(proc.time()[3] - t)
      }

      # check for equality of clusterings
      eq_1 = equal(cl_1)
      eq_2 = equal(cl_2)
      InStab <- mean(eq_1 != eq_2)

      # Normalize = FALSE
      m_instab[b, which(kseq==k)] <- InStab

      # Normalize = TRUE
      M_1 = table(cl_1)
      M_2 = table(cl_2)
      if(predict == FALSE) {
        M_1 = M_1 * (nrow(data) / sum(M_1))
        M_2 = M_2 * (nrow(data) / sum(M_2))
        }
      cs <- instab_simple_var(M_1, M_2)
      m_instab_norm[b, which(kseq==k)] <- .5 * (InStab - cs[1]) / cs[2]

    } # end for k

    if(pbar) utils::setTxtProgressBar(pb, b)

  } # end for B

  # taking the mean
  m_instab_M      = colMeans(m_instab)
  m_instab_norm_M = colMeans(m_instab_norm)

  kopt_instab  = which(m_instab_M == min(m_instab_M))+(min(kseq)-1)
  kopt_instabN = which(m_instab_norm_M == min(m_instab_norm_M))+(min(kseq)-1)

  # replicate function call
  f_call <- list('kseq'=kseq,
                 'nB'=nB,
                 'norm'=norm,
                 'prediction'=predict,
                 'method'=method,
                 'linkage'=linkage,
                 'kmIter'=kmIter)

  outlist <- list("k_instab"=kopt_instab,
                  "k_instab_norm"=kopt_instabN,
                  "instab_path"=m_instab_M,
                  "instab_path_norm"=m_instab_norm_M,
                  "instab_path_matrix"=m_instab,
                  "instab_path_nrom_matrix"=m_instab_norm,
                  'call'=f_call)

  class(outlist) <- c('list', 'cstab', 'cStability')

  return(outlist)

} # EoF

#' Selection of number of clusters via clustering instability
#'
#' Selection of number of clusters via \emph{model-based} or \emph{model-free},
#'   \emph{normalized} or \emph{unnormalized} clustering instability.
#'
#' @param data a n x p data matrix of type numeric.
#' @param kseq a vector with considered numbers clusters k > 1
#' @param nB an integer specifying the number of bootstrap comparisons.
#' @param norm logical specifying whether the instability path should be
#'    normalized. If TRUE, the instability path is normalized, accounting for a
#'    trivial decrease in instability due to a increasing k (see Haslbeck & Wulff,
#'    2016).
#' @param predict boolean specifying whether the model-based or the model-free
#'    variant should be used (see Haslbeck & Wulff, 2016).
#' @param method character string specifying the clustering algorithm. 'kmeans' for
#'    the k-means algorithm, 'hierarchical' for hierarchical clustering.
#' @param linkage character specifying the linkage criterion, in case
#'    \code{type='hierarchical'}. The available options are "single", "complete",
#'    "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median". See
#'    \link[stats]{hclust}.
#' @param kmIter integer specifying the the number of restarts of the k-means algorithm
#'    in order to avoid local minima.
#' @param pbar logical
#'
#' @return a list that contains the optimal k selected by the unnormalized and
#'    normalized instability method. It also includes a vector containing the averaged
#'    instability path (over bootstrap samples) and a matrix containing the instability
#'    path of each bootstrap sample for both the normalized and the unnormalized method.
#'
#' @references
#'    Ben-Hur, A., Elisseeff, A., & Guyon, I. (2001). A stability based method for
#'    discovering structure in clustered data. \emph{Pacific symposium on biocomputing,
#'    7}, 6-17.
#'
#'    Tibshirani, R., & Walther, G. (2005). Cluster validation by prediction strength.
#'    \emph{Journal of Computational and Graphical Statistics, 14}(3), 511-528.
#'
#' @author
#'    Dirk U. Wulff <dirk.wulff@gmail.com>
#'    Jonas M. B. Haslbeck <jonas.haslbeck@gmail.com>
#'
#' @examples
#' \dontrun{
#'   # Generate Data from Gaussian Mixture
#'   s <- .1
#'   n <- 50
#'   data <- rbind(cbind(rnorm(n, 0, s), rnorm(n, 0, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 0, s), rnorm(n, 1, s)),
#'                 cbind(rnorm(n, 1, s), rnorm(n, 0, s)))
#'   plot(data)
#'
#'   # Selection of Number of Clusters using Instability-based Measures
#'   stab_obj <- cStability(data, kseq=2:10)
#'   print(stab_obj)
#'   }
#'
#' @export

cStability_mEst <- function(data, # n x p data matrix
                               kseq = 2:20, # sequence of ks tested
                               nB   = 10, # number of bootstrap comparisons
                               norm = TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                               predict = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                               method    = 'kmeans', # or 'hierarchical'
                               linkage = 'complete', # or average, or ...
                               kmIter  = 5,
                               pbar = TRUE) # number of reruns of k-means algorithm
{


  # ---------- Input Checks ----------

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')

  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')

  # On B
  if(round(nB)!=nB) stop('The number of bootstrap comparison has to be a positive integer value.')

  # On type
  if(method %in% c('kmeans', 'hierarchical'))


    # ---------- Create Containers ----------

  n_obj <- nrow(data)
  m_instab <- matrix(NA, nB, length(kseq)) # storage: no normalization
  m_exp_instab <- numeric(length(kseq)) # storage: no normalization

  # ---------- Draw bootstrap samples ----------

  ind <- list()
  share <- list()
  for(b in 1:nB) {
    tmp_1 <- sample(1:n_obj, n_obj, replace=T)
    tmp_2 <- sample(1:n_obj, n_obj, replace=T)
    tmp_1 <- tmp_1[order(tmp_1)]
    tmp_2 <- tmp_2[order(tmp_2)]
    ind[[b]] <- list(tmp_1,tmp_2)
    intersct <- intersect(tmp_1,tmp_2)
    share[[b]] <- list(tmp_1 %in% intersct & !duplicated(tmp_1), tmp_2 %in% intersct & !duplicated(tmp_2)) # leave duplicates in ?
  }


  # ---------- Calculate distance matrix ----------

  distm = fast_dist(data) # maybe here?

  # ----- Loop over B comparisons -----

  if(pbar)  pb <- utils::txtProgressBar(min=0, max=max(kseq), style = 2)

  for(k in kseq) {

    # container for ms
    m_ms <- matrix(NA, nB * 2, k) # storage: normalization
    new_ns <- rep(0, nB)

    for(b in 1:nB) {

      # kmeans
      if(method == 'kmeans') {
        kms_1 <- list()
        kms_2 <- list()
        for(km_i in 1:kmIter) {
          kms_1[[km_i]] = stats::kmeans(data[ind[[b]][[1]],], centers = k)
          kms_2[[km_i]] = stats::kmeans(data[ind[[b]][[2]],], centers = k) # this has to be 2
        }
        km_1 = kms_1[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
        km_2 = kms_2[[which.min(sapply(kms_1,function(x) mean(x$tot.withinss)))]]
        if(predict == TRUE) {
          cl_1 = kmeans_predict(km_1, data = data)
          cl_2 = kmeans_predict(km_2, data = data)
        } else {
          cl_1 = km_1$cluster[share[[b]][[1]]]
          cl_2 = km_2$cluster[share[[b]][[2]]]
        }
      }

      if(method == 'hierarchical') {
        #t = proc.time()[3]
        hc_1 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[1]],ind[[b]][[1]]]), method = linkage)
        hc_2 = fastcluster::hclust(stats::as.dist(distm[ind[[b]][[2]],ind[[b]][[2]]]), method = linkage)
        cl_1 = stats::cutree(hc_1, k)[share[[b]][[1]]]
        cl_2 = stats::cutree(hc_2, k)[share[[b]][[2]]]
        #print(proc.time()[3] - t)
      }

      # orig n
      new_n = length(cl_1)

      # check for equality of clusterings
      eq_1 = equal(cl_1)
      eq_2 = equal(cl_2)

      # calculate instability
      InStab <- mean(eq_1 != eq_2)

      # Normalize = FALSE
      m_instab[b, which(kseq==k)] <- InStab

      # Normalize = TRUE
      m_1 = table(cl_1)
      m_2 = table(cl_2)
      m_1 = m_1[as.character(1:k)]
      m_2 = m_2[as.character(1:k)]
      m_1[is.na(m_1)] = 0
      m_2[is.na(m_2)] = 0
      m_ms[2*b - 1,] = m_1
      m_ms[2*b,] = m_2

      # store new ns
      new_ns[b] = new_n

    } # end for B

    if(pbar) utils::setTxtProgressBar(pb, k)

    # Gather Ms
    m_ms_M = colMeans(m_ms)
    if(predict == FALSE) m_ms_M = m_ms_M * (nrow(data) / mean(new_ns))
    #print(paste0(sum(m_ms_M),"  -  ",instab_simple(m_ms_M, m_ms_M),"  -  ",paste0(m_ms_M,collapse=' ')))

    # compute instab
    m_exp_instab[which(kseq==k)] = instab_simple(m_ms_M, m_ms_M)

  } # end for B

  # taking the mean
  m_instab_M      = colMeans(m_instab)
  m_instab_norm_M = .5*(m_instab_M - m_exp_instab) / (sqrt(m_exp_instab/2)**2)

  kopt_instab  = which(m_instab_M == min(m_instab_M))+(min(kseq)-1)
  kopt_instabN = which(m_instab_norm_M == min(m_instab_norm_M))+(min(kseq)-1)

  # replicate function call
  f_call <- list('kseq'=kseq,
                 'nB'=nB,
                 'norm'=norm,
                 'prediction'=predict,
                 'method'=method,
                 'linkage'=linkage,
                 'kmIter'=kmIter)

  outlist <- list("k_instab"=kopt_instab,
                  "k_instab_norm"=kopt_instabN,
                  "instab_path"=m_instab_M,
                  "instab_path_norm"=m_instab_norm_M,
                  "instab_path_matrix"=m_instab,
                  "instab_path_norm_matrix"=NA,
                  "c_1"=m_exp_instab,
                  "c_2"=sqrt(m_exp_instab/2)**2,
                  'call'=f_call)

  class(outlist) <- c('list', 'cstab', 'cStability')

  return(outlist)

} # EoF



