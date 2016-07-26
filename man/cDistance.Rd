\name{cDistance}
\alias{cDistance}
\title{
Selection of number of clusters via distance-based measures
}
\description{
Selection of number of clusters via Gap Statistic, Jump Statistic and Slope Statistic
}
\usage{
cDistance(data, kseq, method = "kmeans", linkage = "complete",
          kmIter = 10, RunsGap = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A n x p data matrix
}
  \item{kseq}{
A vector with considered numbers clusters k > 1
}
  \item{method}{
Character string indicating the clustering algorithm. 'kmeans' for the k-means algorithm, 'hierarchical' for hierarchical clustering.
}
  \item{linkage}{
The linkage criterion, in case type='hierarchical'. For available options are "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" or "median".
}
  \item{kmIter}{
The number of restarts of the k-means algorithm in order to avoid local minima.
}
  \item{RunsGap}{
Number of simulated datasets to compute the Gap Statistic (see Tibshirani et al., 2001)
}
}
\value{
Returns a list with the optimal numbers of cluster determined by the Gap Statistic (Tibshirani et al., 2001), the Jump Statistic (Sugar & James, 2011) and the Slope Statistic (Fujita et al., 2014). It also returns the Gap, Jump and Slope for each k in \code{kseq}.
}
\references{
Fujita, A., Takahashi, D. Y., & Patriota, A. G. (2014). A non-parametric method to estimate the number of clusters. Computational Statistics & Data Analysis, 73, 27-39.

Sugar, C. A., & James, G. M. (2011). Finding the number of clusters in a dataset. Journal of the American Statistical Association.

Tibshirani, R., Walther, G., & Hastie, T. (2001). Estimating the number of clusters in a data set via the gap statistic. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 63(2), 411-423.

}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
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

# Selection of Number of Clusters using Distance-based Measures
cDistance(data, kseq=2:10)


}

}

