


### Build the n-1 by n difference matrix
build_D <- function(dat_length){
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  return(D)
}


### Function: calculate ridge regression from closed-form solution
### Input: vector of iwt and incidence
### Return: Rt estimate

get_r <- function(W, Y, lambda=10){
  dat_length = length(Y)
  D = build_D(dat_length)
  W = diag(W)
  return(solve((t(W)%*%W-lambda*t(D)%*%D))%*%t(W)%*%Y)
} 


### Function: LOOCV for ridge regression
### Input: vector of iwt and incidence
### Return: score and lambda from best (lowest score) to worst 

CV <- function(W, Y, lambdas = exp(seq(0.1,10,0.2))){
  ### Length of data
  dat_length = length(Y)
  
  ### Get Difference matrix
  D = build_D(dat_length)
  
  ### Make W into matrix
  W = diag(W)
  
  ### Get grid of lambda
  cv_scores = c()
  best = 0
  
  I = diag(rep(1, dat_length))
  
  for (lambda in lambdas){
    L = solve((t(W)%*%W-lambda*t(D)%*%D))%*%t(W)
    H = I-t(L)%*%W
    H_tilde = diag(diag(H))
    HHY = H%*%solve(H_tilde)%*%Y
    E_cv = 1/(dat_length)*sum(HHY^2)
    cv_scores = c(cv_scores, E_cv)
  }
  output = list()
  ordered = order(cv_scores)
  lambdas = lambdas[ordered]
  

  return(data.frame(scores = cv_scores[ordered], lambdas = lambdas[ordered]))
}