#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> kmeans_predict(GenericVector kmeans, NumericMatrix data) {
  NumericMatrix centers = kmeans[1];
  std::vector<double> dist(data.nrow());
  std::vector<int> clust(data.nrow());
  int ncl = data.ncol(), nrw = data.nrow();
  for(int i = 0; i < centers.nrow(); ++i){
    NumericVector center = centers(i,_);
    for(int j = 0; j < nrw; ++j){
      double d = 0;
      NumericVector row = data(j,_);
      for(int k = 0; k < ncl; ++k){
        double dk = row[k] - center[k];
        d += dk * dk;
        }
      if(i == 0){
        dist[j] = d;
        clust[j] = 1;
        } else {
        if(d < dist[j]){
          dist[j] = d;
          clust[j] = i + 1;
          }
        }
      }
    }
  return clust;
  }
