#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix fast_dist(NumericMatrix x,
                       int power = 2) {
  int ni = x.ncol(), nt = x.nrow();
  double xd,d = 0;
  NumericMatrix dist(nt,nt);
  for(int r = 0; r < nt; r++) {     // loop rows
    for(int c = r; c < nt; c++) {   // loop cols (only half due to symmetry)
      d = 0;                        //dummy value of final matrix entry
      if(power == 1){
        for(int i = 0; i < ni; i++) { // loop 1:n
          xd = x(r,i) - x(c,i);
          d += xd;
        }
      } else if (power == 2){
        for(int i = 0; i < ni; i++) { // loop 1:n
          xd = x(r,i) - x(c,i);
          d += xd*xd;
        }
      } else {
        for(int i = 0; i < ni; i++) { // loop 1:n
          xd = x(r,i) - x(c,i);
          for(int j = 1; j < power; j++){
            xd *= xd;
          }
          d += xd;
        }
      }
      if(power == 2){
        d = sqrt(d);
      } else {
        d = exp( log(d) / double(power));
      }
      dist(r,c) = d;   // fill distance matrix
      dist(c,r) = d;   // fill distance matrix
      }
  }
  return dist;
  }


// [[Rcpp::export]]
std::vector<bool> equal(std::vector<int> x) {
  int n_long = x.size(), n_short = x.size()-1;
  std::vector<bool> res((n_long*n_short)/2);
  int ind = 0;
  for(int i = 0; i < n_short; ++i){
    for(int j = (i+1); j < n_long; ++j){
      res[ind] = x[i] == x[j];
      //std::cout << x[i] << '_' << x[j] << '\n';
      ind++;
      }
    }
  return res;
  }

