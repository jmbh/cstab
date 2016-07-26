\name{cStability}
\alias{cStability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of number of clusters via clustering instability
}
\description{
Selection of number of clusters via model-based or model-free (normalized) clustering instability
}
\usage{
cStability(data, kseq, B = 10, norm = TRUE, prediction = TRUE,
           type = "kmeans", linkage = "complete", kmIter = 10,
           pbar = TRUE)
}
\arguments{
  \item{data}{
A n x p data matrix
}
  \item{kseq}{
A vector with considered numbers clusters k > 1
}
  \item{B}{
Integer value indicating the number of bootstrap comparisons
}
  \item{norm}{
If TRUE, the instability path is normalized, accounting for a trivial decrease in instability due to a increasing k (see Haslbeck & Wulff, 2016)
}
  \item{prediction}{
If TRUE, the model-based method is used (see e.g. Tibshirani & Walter, 2005), if FALSE the model-free method is used (see e.g. Ben-Hur et al., 2001)
}
  \item{type}{
Character string indicating the clustering algorithm. 'kmeans' for the k-means algorithm, 'hierarchical' for hierarchical clustering.
}
  \item{linkage}{
The linkage criterion, in case type='hierarchical'. For available options are "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median".
}
  \item{kmIter}{
The number of restarts of the k-means algorithm in order to avoid local minima.
}
  \item{pbar}{
If TRUE, a progress bar is shown.
}
}
\value{
The function returns a list that contains the optimal k selected by the unnormalized and normalized instability method. It also includes a vector containing the averaged instability path (over bootstrap samples) and a matrix containing the instability path of each bootstrap sample for both the unnormalized and the normalized method.
}
\references{
Ben-Hur, A., Elisseeff, A., & Guyon, I. (2001). A stability based method for discovering structure in clustered data. In Pacific symposium on biocomputing (Vol. 7, pp. 6-17).

Tibshirani, R., & Walther, G. (2005). Cluster validation by prediction strength. Journal of Computational and Graphical Statistics, 14(3), 511-528.
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>, Dirk Wulff
}

\examples{

\dontrun{

  # Generate Data from Gaussian Mixture
  s <- .1
  n <- 50
  data <- rbind(cbind(rnorm(n, 0, s), rnorm(n, 0, s)),
                cbind(rnorm(n, 1, s), rnorm(n, 1, s)),
                cbind(rnorm(n, 0, s), rnorm(n, 1, s)),
                cbind(rnorm(n, 1, s), rnorm(n, 0, s)))
  plot(data)

  # Selection of Number of Clusters using Instability-based Measures
  stab_obj <- cStability(data, kseq=2:10)
  print(stab_obj)

}






}
