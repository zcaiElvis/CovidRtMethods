Rcpp::sourceCpp(here::here('model/tf/c++/gtf.cpp'))
Rcpp::sourceCpp(here::here('model/tf/c++/ptf.cpp'))
Rcpp::sourceCpp(here::here('model/tf/c++/dp-norm.cpp'))
source(here::here('model/tf/r/basics.R'))

# ADMM equivalent solvers for 4 opt problems
gtf_solver <- function(max_iter, y, x, k, u_init=NULL, w_init=NULL, para = 0.01, tol = 1e-3){
  
  n = length(y)
  # initial parameters: 
  lambda = para
  rho = para/2
  theta = numeric(n)
  u = u0 = w = numeric(n-k)
  if(!is.null(u_init)) {
    if(length(u_init) == length(u)) {
      u = u0 = u_init
    } else {
      stop("Error: Imcompatible length of initial u (u_init).")
    }
  }
  if(!is.null(w_init)) {
    if(length(w_init) == length(w)) {
      w = w_init
    } else {
      stop("Error: Imcompatible length of initial w (w_init).")
    }
  }
  
  if (k==0){
    theta = prox_dp_norm(y=y, lam=lambda)
    conver = TRUE
  } else if (k>0){
    D = D_generator(n, k-1) #order k-1
    #matrix inverse (to be used in the theta step): 
    Max_inv = solve(diag(x^2) + 2 * rho * t(D) %*% D)
    mod = gtf_cpp_solver(M=max_iter, y=y, x=x, n=n, theta=theta, u=u, u0=u0, w=w, 
                         lambda=lambda, rho=rho, D=D, Max_inv=Max_inv, tol=tol)
    theta = mod$theta
    conver = (mod$iter_num < max_iter) # check convergence before reaching the max iterate number
  }
  
  return(list(
    theta = theta, convergence = conver
    #u = mod$u, w = mod$w, prim_res = mod$prim_res, dual_res = mod$dual_res, iter_num = mod$iter_num 
  ))
}

ptf_solver <- function(max_iter, y, x, k, u_init=NULL, w_init=NULL, para = 0.01, tol = 1e-3){
  
  n = length(y)
  # initial parameters: 
  lambda = para
  rho = para/2
  D = D_generator(n, k) 
  theta = numeric(n)
  u = u0 = w = numeric(n-k-1)
  if(is.null(u_init) == FALSE) {
    if(length(u_init) == length(u)) {
      u = u0 = u_init
    } else {
      stop("Error: Imcompatible length of initial u (u_init).")
    }
  }
  if(!is.null(w_init)) {
    if(length(w_init) == length(w)) {
      w = w_init
    } else {
      stop("Error: Imcompatible length of initial w (w_init).")
    }
  }
  
  mod = ptf_cpp_solver(M=max_iter, y=y, x=x, n=n, theta=theta, u=u, u0=u0, w=w, 
                         lambda=lambda, rho=rho, D=D, tol=tol)
  conver = (mod$iter_num < max_iter) # check convergence before reaching the max iterate number
  
  return(list(
    theta = mod$theta, convergence = conver
  ))
}
