#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_simplexCENTER_cpp(CharacterVector simplex, NumericMatrix vertices) {
  std::string simplex_str = as<std::string>(simplex[0]);
  std::vector<int> vert_ids;
  std::istringstream ss(simplex_str);
  std::string token;
  
  while (ss >> token) {
    vert_ids.push_back(std::stoi(token));
  }
  
  // Convert to Armadillo matrix for fast operations
  arma::mat vert_mat(vertices.begin(), vertices.nrow(), vertices.ncol(), false);
  
  arma::vec center = arma::zeros<arma::vec>(3);
  int count = 0;
  
  for (int id : vert_ids) {
    int idx = id - 1;
    if (idx >= 0 && idx < vertices.nrow()) {
      center += vert_mat.row(idx).t();  // Fast row access and addition
      count++;
    }
  }
  
  if (count > 0) {
    center /= count;  // Vectorized division
    return NumericVector::create(center(0), center(1), center(2));
  } else {
    return NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  }
}