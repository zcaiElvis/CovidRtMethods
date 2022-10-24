
### Function: Calculating R from closed-from ridge regression solution
### Return: Rt and predicted rss
pls_ridge <- function(iwt, Y, lambda){
  dat_length = length(Y)
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  
  W <- diag(iwt)
  

  R <- solve((t(W)%*%W-lambda*t(D)%*%D))%*%t(W)%*%Y
  rss <- sum(Y-W%*%R)^2 + sum(lambda*(D%*%R)^2) 
  
  output = list()
  output$r = R
  output$rss = rss
  
  return(output)
}

### Function: LOOCV for ridge regression, return cv score and lambda by cv ascending order
### Return 
ridge_cv <- function(iwt, Y, lambdas = exp(seq(0.1,10,0.2))){
  ### Length of data
  dat_length = length(Y)
  
  Y <- as.matrix(Y, ncol=dat_length)
  
  ### Construct difference matrix
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  
  
  ### Get grid of lambda
  cv_scores = c()
  
  I = diag(rep(1, dat_length))
  
  W = diag(iwt)
  
  for (lambda in lambdas){
    L = solve((t(W)%*%W-lambda*t(D)%*%D))%*%t(W)
    H = I-t(L)%*%W
    H_tilde = diag(diag(H))
    HHY = H%*%solve(H_tilde)%*%Y
    E_cv = 1/(dat_length)*sum(HHY^2)
    cv_scores = c(cv_scores, E_cv)
  }
  
  ordered = order(cv_scores)
  lambdas = lambdas[ordered]
  
  output= list()
  output$score = ordered
  output$lambdas = lambdas
  output$optim_lambda = lambdas[1]
  
  return(output)
}


# plot(ridge_cv(1:10, 2*c(1:10))$score)
# plot(pls_ridge(1:10, 2*c(1:10), lambda = 220)$r)