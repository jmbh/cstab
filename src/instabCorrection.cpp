#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rootChoose(int n, int k, double root){
  double nomin = 1, denom = 1;
  
  // comp nomin
  for(int i = n; i > k; i--){
    nomin *= std::pow(i,1/root);
    }
  for(int i = n-k; i > 0; i--){
    denom *= std::pow(i,1/root);
    }
  return(nomin/denom);
  }

// [[Rcpp::export]]
double rootChooseLookup(int n, int k, std::vector<double> lookup){
  double res, root; 
  int nLookup = lookup.size();
  if(n > (nLookup+1)){
    root = lookup[0];
    double nomin = 1, denom = 1;
    // comp nomin
    for(int i = n; i > k; i--){
      nomin *= std::pow(i,1/root);
      }
    for(int i = n-k; i > 0; i--){
      denom *= std::pow(i,1/root);
      }
    res = nomin/denom;
    } else {
    if(n != k & k > 0){
      res = lookup[n] / (lookup[k]*lookup[n-k]);
      } else {
      res = 1;
    }
  }
  return(res);
  }


double rootFact(int n, double root){
  double v = 1;
  for(int i = n; i > 0; i--){
    v *= pow(i,1/root);
    }
  return(v);
  }


// [[Rcpp::export]]
std::vector<double> lookup(int n, double root){
  std::vector<double> lookupTab;
  lookupTab.push_back(root);
  for(int i = 0; i < (n+1); i++){
    lookupTab.push_back(rootFact(i+1,root));
    }
  return(lookupTab);
  }


// [[Rcpp::export]]
double rootCombLookup(std::vector<double> ns,  std::vector<double> lookup){
  double res = 1;
  int i, n = 0;
  for(i = 0; i < ns.size(); i++) n += ns[i];
  for(i = 0; i < ns.size(); i++){ 
    res *= rootChooseLookup(n,ns[i],lookup);
    //std::cout << res << '\n';
    n -= ns[i];
    }
  return(res);
  }

// [[Rcpp::export]]
double stabExp(std::vector<double> ns,  std::vector<double> lookup){
  double res = 0;
  int i, k, n = 0;
  std::vector<double> nstmp;
  double totComb = rootCombLookup(ns,lookup);
  for(i = 0; i < ns.size(); i++) n += ns[i];  
  for(i = 0; i < ns.size(); i++){
    k = round(ns[i]);
    if(k > 1){
      for(int j = 0; j < ns.size(); j++) {
        if(j != i){
          nstmp.push_back(ns[j]);
          }
        }
      double rootCo     = rootChooseLookup(n-2,k-2,lookup) * rootCombLookup(nstmp,lookup);
      //std::cout << rootCo << '\n';
      double normRootCo = rootCo / totComb;
      double normCo     = std::pow(normRootCo,lookup[0]);
      res += normCo;
      nstmp.clear();
    }
  }
  return(res);
}



