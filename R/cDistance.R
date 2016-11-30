#' Selection of number of clusters via distance-based measures
#'
#' Selection of number of clusters via \emph{gap statistic}, \emph{jump statistic}, and
#'    \emph{slope statistic}
#'
#' @param data a n x p data matrix of type numeric.
#' @param kseq a vector with considered numbers clusters k > 1
#' @param method character string indicating the clustering algorithm. 'kmeans' for the
#'    k-means algorithm, 'hierarchical' for hierarchical clustering.
#' @param linkage character specifying the linkage criterion, in case
#'    \code{type='hierarchical'}. The available options are "single", "complete",
#'    "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median". See
#'    \link[stats]{hclust}.
#' @param kmIter integer specifying the the number of restarts of the k-means algorithm
#'    in order to avoid local minima.
#' @param gapIter integer specifying the number of simulated datasets to compute the
#'    \emph{gap statistic} (see Tibshirani et al., 2001).
#'
#' @return a list with the optimal numbers of cluster determined by the \emph{gap statistic}
#'    (Tibshirani et al., 2001), the \emph{jump Statistic} (Sugar & James, 2011) and the
#'    \emph{slope statistic} (Fujita et al., 2014). Along the function returns the \emph{gap},
#'    \emph{jump} and \code{slope} for each k in \code{kseq}.
#'
#' @references
#'   Tibshirani, R., Walther, G., & Hastie, T. (2001). Estimating the number of clusters in a
#'   data set via the gap statistic. \emph{Journal of the Royal Statistical Society: Series B
#'   (Statistical Methodology), 63}(2), 411-423.
#'
#'   Sugar, C. A., & James, G. M. (2011). Finding the number of clusters in a dataset. \emph{Journal
#'   of the American Statistical Association, 98}(463), 750-763,
#'
#'   Fujita, A., Takahashi, D. Y., & Patriota, A. G. (2014). A non-parametric method to estimate
#'   the number of clusters. \emph{Computational Statistics & Data Analysis, 73}, 27-39.
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>
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
#'  # Selection of Number of Clusters using Distance-based Measures
#'  cDistance(data, kseq=2:10)
#'  }
#'
#' @export

cDistance <- function(data, # n x p data matrix
                      kseq, #sequence of ks to be checked
                      method = 'kmeans', # or: hiearchical
                      linkage = 'complete',
                      kmIter = 10, # restarts of k means algorithm
                      gapIter = 10) # number of simulated datasets in gap statistic
{

  # ----- INPUT TESTS

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')

  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')

  # On type
  if(method %in% c('kmeans', 'hierarchical'))


  # ----- HELPERS

  n    = nrow(data)
  dims = ncol(data)
  if(!1 %in% kseq) kseq = c(1,kseq)


  # ----- EVALUATE REAL DATA

  WCD <- Sil <- MSE <- numeric()
  for(k in kseq) {
    obj = getMeasures(data, k, method=method, linkage=linkage, kmIter = kmIter, measures = c('wcd','sil','mse'))
    WCD[k] = obj$WCD
    Sil[k] = obj$Sil
    MSE[k] = obj$MSE
    }


  # ----- EVALUATE SYNTHETIC DATA (Gap-statistic)

  WCD_runs = matrix(NA,nrow=gapIter, ncol=length(kseq))
  for(i in 1:gapIter) {
    data_syn = UniformData(data)
    WCDs = numeric()
    for(j in 1:length(kseq)) {
      k = kseq[j]
      obj = getMeasures(data_syn, k, method=method, linkage=linkage, kmIter = kmIter, measures = c('wcd'))
      WCDs[j] = obj$WCD
      }
    WCD_runs[i,] = WCDs
    }
  WCD_syn = colMeans(WCD_runs)

  # ----- COMPUTE MEASURES

  # Gap Statistic
  WCD_dat_log = log(WCD)
  WCD_syn_log = log(WCD_syn)
  WCD_dat_log = WCD_dat_log - WCD_dat_log[1]
  WCD_syn_log = WCD_syn_log - WCD_syn_log[1]
  gap      = WCD_syn_log - WCD_dat_log
  kopt_gap = kseq[gap == max(gap)]


  # Slope Statistic
  p = 1
  slope = -(Sil[-1] - Sil[-length(Sil)]) * Sil[-1]^p
  kopt_slope = kseq[slope == max(slope)]

  ## Jump Statistic
  MSE_tr = MSE^(- dims/2)
  jump   = (MSE_tr - c(0, MSE_tr[-length(MSE_tr)])) #[-1]
  kopt_jump = kseq[jump == max(jump)]

  outlist <- list('k_Gap'=kopt_gap,
                  'k_Slope'=kopt_slope,
                  'k_Jump'=kopt_jump,
                  'WCD'=WCD,
                  'WCD_syn'=WCD_syn,
                  'Gaps'= gap,
                  'Silhouettes'=Sil,
                  'Slopes'=slope,
                  'Jumps'=jump)

  return(outlist)

} # EoF


