library("pomp")



### Transition equation

tran_lognormal <- discrete_time(
  function(x, sdlog, ...){
    xnext <- rlnorm(1, meanlog = log(x), sdlog = sdlog)
    c(x = xnext)
  },
  delta.t = 1
)



### Observation equation

obs_pois <- function(t, x, T, Y, y, g_shape, g_scale, ..., log){
  if(t > 1){
    lambda = x * sum(Y[(t-1):1] * disc_gamma(1:(t-1), shape = g_shape, scale = g_scale))
    dpois(round(y), lambda = lambda)
  }else{
    dpois(round(y), lambda = 1)
  }
}

### Synthetic dataset

d1 <- read.csv("data/processed/d1.csv")
plot(d1$y)




pf <- pfilter(
  Np=100,
  times = "idx",
  t0 = 1,
  data = d1
)