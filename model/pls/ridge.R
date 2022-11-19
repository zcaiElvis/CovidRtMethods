library(dplyr)
library(zoo)

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


### Function: return score of a lambda
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
  
  return(data.frame(lambdas = lambdas, scores = cv_scores))
}









cv_loss <- function(w, Y, lambdas = exp(seq(0.1,10,0.2))){
  losses = c()
  for(lambda in lambdas){
    r = get_r(w, Y, lambda)
    loss = mean((Y - w*r)^2 + sum((diff(r))^2))
    losses = c(losses, loss)
  }
  
  best_lambda = lambdas[which.min(losses)]
  output = list()
  output$best_lambda = best_lambda
  output$losses = losses
  return(output)
}


cv_genlasso <- function(w, Y, k = 3, lambdas = exp(seq(0.1, 10, 0.2))){
  dat_length = length(Y)
  kindex = c(1:dat_length)%%k + 1
  data = data.frame(idx = 1:dat_length, kindex = kindex, w = w, Y = Y)
  scores = matrix(0, nrow = k, ncol = length(lambdas))
  
  for(i in 1:k){
    
    for(j in 1:length(lambdas)){
      
      train = data[data$kindex != i,]
      test = filter(data, kindex == i)
      
      pred_r <- get_r(W= train$w, Y = train$Y, lambda = lambdas[j])
      
      r <- rep(0, dat_length)
      r[train$idx] = pred_r
      r[test$idx] = NA
      r <- na.fill(r, "extend")
      # scores[i,j] = sum(Y-w*r)^2 + lambda*sum(diff(r)^2)
      filled_r <- r[test$idx]
      
      scores[i,j] = sum(test$Y - test$w*filled_r)^2 + lambdas[j]*sum(diff(filled_r)^2)
    }
  }
  scores = apply(scores, MARGIN = 2, mean)
  
  output = data.frame(lambdas = lambdas, scores = scores)
  return(output)
}

# d = read.csv("data/processed/d2.csv")
# genl <- cv_genlasso(d$iwt, d$y, k = 5, lambda = exp(1:30))
