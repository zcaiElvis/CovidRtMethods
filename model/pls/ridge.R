

### Function: Build the n-1 by n difference matrix
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
  return(solve((t(W)%*%W+lambda*t(D)%*%D))%*%t(W)%*%Y)
}



### Get Loss function value
get_loss <- function(W, Y, r, lambda){
  dat_length = length(Y)
  D = build_D(dat_length)
  W = diag(W)
  
  loss = sum((Y-W%*%r)^2/dat_length) + lambda*sum((D%*%r)^2)
  return(loss)
}




### Function: get hat matrix
### H = W'(W'W+lam*D'D)W
get_hat <- function(W, D, lambda){
  H = W %*% solve(t(W)%*%W + lambda*t(D)%*%D)%*%t(W)
  return(H)
}


get_score <- function(W, Y, D, I, lambda){
  H = get_hat(W, D, lambda)
  score = mean((((I-H)%*%Y)/(1-diag(H)))^2)
  return(score)
}



### Function: LOOCV for ridge regression
### Input: vector of iwt and incidence
### Return: score and lambda from best (lowest score) to worst 

CV <- function(W, Y, lambdas = exp(seq(0.1,10,0.2))){
  
  ### Formatting data
  dat_length = length(Y)
  W = diag(W)
  D = build_D(dat_length)
  I = diag(dat_length)

  ### Get grid of lambda
  cv_scores = c()
  
  ## Loop through all possible lambdas
  for (lambda in lambdas){
    E_cv = get_score(W, Y, D, I, lambda)
    cv_scores = c(cv_scores, E_cv)
  }
  
  ## Sort the lambda based on size of cv_scores
  output = list()
  ordered = order(cv_scores)
  lambdas = lambdas[ordered]
  
  ### Return data frame of sorted cv_scores and lambdas
  return(data.frame(scores = cv_scores[ordered], lambdas = lambdas[ordered]))
}