# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List ptf_cpp_solver (int M, vec y, vec x, int n, vec theta, vec u, vec u0, vec w, double lambda, double rho, mat D, 
  double tol = 1e-3){
  
  int iter;
  double r_norm;
  double s_norm;
  mat MatD (n,n);
  mat W (n,n); // Hessian
  vec V (n); // gradients
  mat Max_inv (n,n);

  MatD = 2 * rho * trans(D) * D;

  // start of iteration: 
  for(iter = 0; iter < M; iter++) {

    // update primal variable:  
    // -- matrix inverse: 
    for(int i=0; i<n; i++){
      W(i,i) = pow(x[i],2) * exp(x[i] * theta[i]) / n;
    }
    Max_inv = MatD - W;
    Max_inv = inv(Max_inv);
    // -- gradients:
    for(int i=0; i<n; i++){
        V[i] = y[i] * x[i] / n - x[i] * exp(x[i] * theta[i]) / n;
    }
    // -- update:
    theta = Max_inv * (2 * rho * trans(D) * (u + w) - W * theta - V);
    
    // update alternating variable:   
    vec D_u = D * theta - w;
    vec vecNull (n-2, fill::value(0.0));
    vec v2max = abs(D_u) - lambda / (2 * rho);
    u = sign(D_u) % arma::max(vecNull, v2max); 
    
    // update dual variable: 
    w += u - D * theta;

    // stopping criteria check: 
    // primal residuals: 
    vec r = D * theta - u;
    r_norm = sqrt(sum(square(r)));
    // dual residuals: 
    vec s = u0 - u;
    s_norm = sqrt(sum(square(s)));

    if( (r_norm < tol) && (s_norm < tol) ){
      iter++;
      break;
    }
    // auxiliary variables (for checking stopping criteria) update: 
    u0 = u;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("theta") = wrap(theta), Rcpp::Named("u") = wrap(u), Rcpp::Named("w") = wrap(w), 
    Rcpp::Named("prim_res") = r_norm, Rcpp::Named("dual_res") = s_norm, Rcpp::Named("iter_num") = iter
  );
}