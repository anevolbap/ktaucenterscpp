#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
List distance_to_centers(NumericMatrix data, NumericMatrix centers) {
  /* data is a matrix with observations(rows) and variables, 
   centers is a matrix with cluster centers coordinates (rows) 
   */
    
  int k = centers.rows(); // number of centers
  int n = data.rows(); // number of observations
  int p = data.cols(); // number of variables
  IntegerVector membership(n); // membership vector
  NumericVector min_distance(n); // distance vector to closest center
  NumericMatrix distance_matrix(n, k); // return value  
  
  for(int n_iter = 0; n_iter < n; n_iter++){
    double min_dist_aux = R_PosInf;
    for(int k_iter = 0; k_iter < k; k_iter++){
      double dist=0;
      for(int p_iter = 0; p_iter < p; p_iter++){
        dist = dist + pow(data(n_iter, p_iter) - centers(k_iter, p_iter), 2);
      }
      dist = sqrt(dist);
      distance_matrix(n_iter, k_iter) = dist;
      if (dist < min_dist_aux) {
        min_distance(n_iter) = dist;
        membership(n_iter) = k_iter + 1;
        min_dist_aux = dist;
      }
    }
  }
  return(List::create(_["distance_matrix"]=distance_matrix, 
                      _["membership"]=membership,
                      _["min_distance"]=min_distance));
}