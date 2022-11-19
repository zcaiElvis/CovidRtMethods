library(testthat)

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
  z = D%*%r
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
  expect(dat_length == length(r_now), "length of data not equal to length of R")
  
  a = rep(mu, dat_length)
  b = rho*(DD%*%r_now - t(D)%*%z + t(D)%*%u)/r_now - mu*r_now +  w/dat_length
  c = -y/dat_length
  
  r = (-b + sqrt(b^2 - 4*a*c))/(2*a)
  return(r)
}



### Function: update z
update_z <- function(r, w, z_now, u, D){
  z_new = sign(z_now)*(abs(z_now) - (D%*%r-u))
}

### Function: update u
update_u <- function(r, y, z, u, D){
  return(u+ D%*%r-z)
}


### Function: helper, build D matrix
build_D <- function(dat_length){
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  return(D)
}



d = read.csv("d2.csv")

r = trend_filter(d$iwt, d$y, rho = 0.03, num_iter = 1000)

ggplot(data.frame(idx = d$idx, true = d$r, pred = r))+
  geom_line(aes(x=idx, y = true), color = "red")+
  geom_line(aes(x=idx, y = pred), color = "blue")
