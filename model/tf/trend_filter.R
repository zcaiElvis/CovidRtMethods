library(testthat)
library(progress)
### Function: calculate r given w and y
### Parameters:
###   w, y: weighted case count and case count
###   k: degree of the difference matrix
###   rho: penalty of smoothness
###   num_iter: number of iterations
###   atol: absolute tolerance for convergence

trend_filter <- function(w, y, rho = 0.1, num_iter = 20){
  k=0
  dat_length = length(y)
  
  r = rep(1, dat_length) # initialize r
  rnext = r # next r

  D = build_D(dat_length)
  DD = t(D)%*% D
  
  mu <- (svd(t(D)%*%D)$d[1])^2 # largest eigen value of D^TD
  
  ### Initialize z and u
  z = matrix(D%*%r)
  u = as.matrix(rep(1, dat_length-k-1))
  
  counter = 0
  while(T){
    rnext = update_r(r, w, y, z, u, rho, mu, D, DD, dat_length)
    znext = update_z(rnext, w, z, u, D)
    unext = update_u(rnext, y, znext, u, D)
    
    counter = counter + 1
    if(counter > num_iter){break}
    
    r = rnext
    z = znext
    u = unext
  }
  
  return(rnext)
}

### Function: update r
update_r <- function(r_now, w, y, z, u, rho, mu, D, DD, dat_length){
  a = rep(mu, dat_length)
  b = rho*(DD%*%r_now - t(D)%*%z + t(D)%*%u) - mu*r_now +  w/dat_length
  c = -y/dat_length
  r = (-b + sqrt(b^2 - 4*a*c))/(2*a)
  return(r)
}



### Function: update z
update_z <- function(r, w, z_now, u, D){
  z_new = sign(z_now)*(abs(z_now) - (D%*%r-u))  ## TODO: Dlog(r)-z for Poisson
}

### Function: update u
update_u <- function(r, y, z, u, D){
  return(u+ D%*%r-z) ## TODO: Dlog(r)-z for Poisson
}


### Function: helper, build D matrix
build_D <- function(dat_length){
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  return(D)
}

### Cross validation for trendfilter
# Can take any function that compute r
# Can take any form of loss
cv_genlasso <- function(w, Y, k = 3, rhos = exp(seq(-10, -1, 0.5))){
  dat_length = length(Y)
  kindex = c(1:dat_length)%%k + 1
  data = data.frame(idx = 1:dat_length, kindex = kindex, w = w, Y = Y)
  scores = matrix(0, nrow = k, ncol = length(rhos))
  pb <- progress_bar$new(total = length(rhos))
  
  for(i in 1:k){
    
    for(j in 1:length(rhos)){
      
      train = data[data$kindex != i,]
      test = filter(data, kindex == i)
      
      pred_r <- trend_filter(train$w, train$Y, rho = rhos[j], num_iter=500)
      
      r <- rep(0, dat_length)
      r[train$idx] = pred_r
      r[test$idx] = NA
      r <- na.fill(r, "extend")
      # scores[i,j] = sum(Y-w*r)^2 + lambda*sum(diff(r)^2)
      filled_r <- r[test$idx]
      
      # scores[i,j] = sum(test$Y - test$w*filled_r)^2 + rhos[j]*sum(abs(diff(filled_r)))
      scores[i,j] = sum(test$w*filled_r - test$Y*log(test$w*filled_r))/dat_length + rhos[j]*sum(abs(diff(filled_r)))
    }
    pb$tick()
  }
  scores = apply(scores, MARGIN = 2, mean)
  
  output = data.frame(rhos = rhos, scores = scores)
  return(output)
}


# caresult <- trend_filter(ca$iwt, ca$y, rho = 0.01, num_iter = 1000)
# plot(caresult)

# d = read.csv("../data/processed/d2.csv")

# r = trend_filter(d$iwt, d$y, rho = 0.03, num_iter = 1000)

# cv_tf <- cv_genlasso(d$iwt, d$y, k = 3)
# 
# ggplot(data.frame(idx = d$idx, true = d$r, pred = r))+
#   geom_line(aes(x=idx, y = true), color = "red")+
#   geom_line(aes(x=idx, y = pred), color = "blue")
