# include <RcppArmadillo.h>
# include "dp-norm.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List gtf_cpp_solver (int M, vec y, vec x, int n, vec theta, vec u, vec u0, vec w, 
            double lambda, double rho, mat D, mat Max_inv, double tol = 1e-3){
  
  int iter;
  double r_norm, s_norm;
  // start of iteration: 
  for(iter = 0; iter < M; iter++) { 

    // update primal variable: 
    vec xy (n);
    for(int i=0; i<n; i++){
      xy[i] = x[i] * y[i];
    }
    theta = Max_inv * (2 * rho * trans(D) * (u + w) + xy);

    // update alternating variable:   
    vec y_dp = D * theta - w;
    double lam_dp = lambda / (2*rho);
    u = prox_dp_norm(y_dp, lam_dp); 
    
    // update dual variable: 
    w += u - D * theta;

    // stopping criteria check: 
    vec r = D * theta - u;
    r_norm = sqrt(sum(square(r)));
    // dual residuals: 
    vec s = u0 - u;
    s_norm = sqrt(sum(square(s)));

    if( (r_norm < tol) && (s_norm < tol) ){
      iter++;
      break;
    }
    // auxiliary variables update: 
    u0 = u;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("theta") = wrap(theta), Rcpp::Named("u") = wrap(u), Rcpp::Named("w") = wrap(w), 
    Rcpp::Named("prim_res") = r_norm, Rcpp::Named("dual_res") = s_norm, Rcpp::Named("iter_num") = iter
  );
}