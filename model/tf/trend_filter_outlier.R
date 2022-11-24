library(testthat)

### Function: calculate r given w and y
### Parameters:
###   w, y: weighted case count and case count
###   k: degree of the difference matrix
###   rho: penalty of smoothness
###   num_iter: number of iterations
###   atol: absolute tolerance for convergence

trend_filter_ol <- function(w, y, rho = 0.1, gamm=0.1, num_iter = 1, atol=0.01){
  k=0
  dat_length = length(y)
  
  r = rep(1, dat_length)  # initialize r
  rnext = r # next r
  
  D = build_D(dat_length)
  DD = t(D)%*% D
  
  mu <- (svd(t(D)%*%D)$d[1])^2 # largest eigen value of D^TD
  
  ### Initialize z and u
  z = D%*%r
  u = as.matrix(rep(1, dat_length-k-1))
  
  ### Initialize o
  o = rep(0, dat_length)
  
  counter = 1
  while(T){
    rnext = update_r_ol(r, w, y, z, u, o, rho, mu, D, DD, dat_length)
    onext = update_o_ol(rnext, w, y, z, u, D, gamm, dat_length)
    znext = update_z_ol(rnext, w, z, u, D)
    unext = update_u_ol(rnext, y, znext, u, D)
    

    counter = counter + 1
    if(counter > num_iter){break}
    # if(abs(sum(rnext-r)) < atol){break}
    
    r = rnext
    o = onext
    z = znext
    u = unext
  }
  
  output = list()
  output$r = rnext
  output$o = onext
  output$z = znext
  output$u = unext
  output$num_iter = counter
  
  return(output)
}


### Function: update r
# update_r <- function(r_now, w, y, z, u, rho, mu, D, DD, dat_length){
#   expect(dat_length == length(r_now), "length of data not equal to length of R")
# 
#   a = rep(mu, dat_length)
#   b = rho*(DD%*%r_now - t(D)%*%z + t(D)%*%u)/r_now - mu*r_now +  w/dat_length
#   c = -y/dat_length
# 
#   r = (-b + sqrt(b^2 - 4*a*c))/(2*a)
#   return(r)
# }

### Function: update r
update_r_ol <- function(r_now, w, y, z, u, o, rho, mu, D, DD, dat_length){
  # expect(dat_length == length(r_now), "length of data not equal to length of R")

  l = rho*(DD%*%r_now - t(D)%*%z + t(D)%*%u) - mu*r_now +  w/dat_length
  a = mu*w
  b = (-mu*o + l*w)
  c = -l*o - y*w/dat_length

  r = (-b + sqrt(b^2 - 4*a*c))/(2*a)
  return(r)
}



### Function: update z
update_z_ol <- function(r, w, z_now, u, D){
  z_new = sign(z_now)*(abs(z_now) - (D%*%log(r)-u))  ## TODO: Dlog(r)-z for Poisson
  return(z_new)
}

### Function: update u
update_u_ol <- function(r, y, z, u, D){
  return(u+ D%*%log(r)-z) ## TODO: Dlog(r)-z for Poisson
}


### Function: update o
update_o_ol <- function(r, w, y, z, u, D, gamm, dat_length){
  a = -gamm
  b = (gamm*w*r+1/dat_length)
  c = -w*r/dat_length + y/dat_length
  o = (-b + sqrt(b^2 - 4*a*c))/(2*a)
  # print((b^2)[1:50])
  # print((4*a*c)[1:50])
  return(o)
}


### Function: helper, build D matrix
build_D <- function(dat_length){
  D = diag(dat_length)
  D[row(D) == col(D)-1] = -1
  D = D[1:(nrow(D)-1),]
  return(D)
}

cv_genlasso_ol <- function(w, Y, k = 3, rhos = exp(seq(-10, -1, 2))){
  dat_length = length(Y)
  kindex = c(1:dat_length)%%k + 1
  data = data.frame(idx = 1:dat_length, kindex = kindex, w = w, Y = Y)
  # scores = matrix(0, nrow = k, ncol = length(rhos))
  scores = array(rep(0, length(rhos)*length(rhos)*k), dim=c(length(rhos), length(rhos), k))
  pb <- progress_bar$new(total = length(rhos))
  
  for(j in 1:length(rhos)){
    
    for(l in 1:length(rhos)){
      
      for(i in 1:k){
        train = data[data$kindex != i,]
        test = filter(data, kindex == i)
        
        pred_r <- trend_filter_ol(train$w, train$Y, rho = rhos[j], gamm =rhos[l], num_iter=500)
        pred_r = pred_r$r
        
        r <- rep(0, dat_length)
        r[train$idx] = pred_r
        r[test$idx] = NA
        r <- na.fill(r, "extend")
        # scores[i,j] = sum(Y-w*r)^2 + lambda*sum(diff(r)^2)
        filled_r <- r[test$idx]
        
        # scores[i,j] = sum(test$Y - test$w*filled_r)^2 + rhos[j]*sum(abs(diff(filled_r)))
        scores[j, l, i] = sum(test$w*filled_r - test$Y*log(test$w*filled_r))/dat_length + rhos[j]*sum(abs(diff(filled_r)))
      }
      
    }
    pb$tick()
  }
  # print(scores)
  # scores = apply(scores, MARGIN = 2, mean)
  
  # output = data.frame(rhos = rhos, scores = scores)
  # return(output)
  print(scores)
}


# caresult= trend_filter_ol(ca$iwt, ca$y, rho = 1, gamm=1, num_iter = 200)


d = read.csv("../data/processed/2022-11-20 21:43:18/d1.csv")

result = trend_filter_ol(d$iwt, d$y, rho = 0.1, gamm=0.000001, num_iter = 2000)


# result$num_iter

# cv_ol <- cv_genlasso_ol(d$iwt, d$y, k=2)
# 
# ggplot(data.frame(idx = d$idx, true = d$r, pred = result$r))+
#   geom_line(aes(x=idx, y = true), color = "red")+
#   geom_line(aes(x=idx, y = pred), color = "blue")
# 
# plot(d$y, type = "l")
# plot(result$o, type = "l")
# plot(d$iwt*result$r - d$y)
