library(tidyverse)
source("model/pls/ridge.R")




### Generate dataset
set.seed(11121)
dat_length = 100
sd = 100
slope = 10

x <- 1:dat_length
y <- slope*x + rnorm(dat_length, 0, sd = sd)
plot(y, type = "l")



### get cv score
cv_result <- CV(x, y, lambda = exp(1:20))

best_lambda_cv = cv_result$lambdas[which.min(cv_result$scores)]

best_r_cv <- get_r(x, y, best_lambda_cv)


### Test
# lambda = cv_result$lambdas[1]
# r <- get_r(x, y, lambda)
# 
# ### Calculate cv score of one lambda
# h = hat_matrix(diag(x), build_D(dat_length),lambda = lambda)
# 
# 
# scores <- function(y, h){return(mean(((y - h%*%y)/(1-diag(h)))^2))}
# 
# 
# hat_matrix <- function(W, D, lambda){
#   return(W %*% solve(t(W)%*%W + lambda*t(D)%*%D) %*% t(W))
# }
# 
# scores(y, h)

losses = c()
for(lambda in cv_result$lambdas){
  r = get_r(x, y, lambda)
  loss = mean((y - r*x)^2) + sum((diff(r))^2)
  losses = c(losses, loss)
}

best_r_loss <- get_r(x,y,lambda = losses[which.min(losses)])


ggplot(data.frame(idx=1:dat_length, cv = best_r_cv, loss = best_r_loss, true = rep(slope, dat_length)), aes(x=idx))+
  geom_line(aes(x=idx, y = cv), color="red")+
  geom_line(aes(x=idx, y = loss), color = "blue")+
  geom_line(aes(x=idx, y = true))+
  ylim(c(0, slope+sqrt(sd)+2))





test_r <- get_r(d$iwt, d$y, 10000000)
plot(test_r, type = "l")