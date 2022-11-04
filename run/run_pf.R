library("pomp")
source("constant/constant.R")
source("function/disc_gamma.R")
source("function/make_plot.R")


### Transition equation

tran_lognormal <- discrete_time(
  function(x, sdlog, ...){
    xnext <- rlnorm(1, meanlog = log(x), sdlog = sdlog)
    c(x = xnext)
  },
  delta.t = 1
)



### Observation equation

meas_pois_interval <- function(t, x, T, Y, y, g_shape, g_scale, ..., log){
  if(t > 1){
    lambda = x * sum(Y[(t-1):1] * disc_gamma(1:(t-1), shape = g_shape, scale = g_scale))
    # sum(dpois(round(0.98*y):round(1.02*y), lambda = lambda))
    y*log(lambda)-lambda
  }else{
    dpois(round(y), lambda = 1)
  }
}

### Synthetic data import

d4 <- read.csv("data/processed/d2.csv")
plot(d4$y)
plot(d4$r)

### Initial values

init_vals = c(sdlog = 0.05, x0 = 1.5,
              g_shape = sid_covid_shape, g_scale = sid_covid_scale)



### Run particle filter default

pf <- pfilter(
  Np=1000,
  times = "idx",
  t0 = 1,
  data = d4,
  rinit = function(x0, ...){c(x=x0)},
  params = init_vals,
  rprocess = tran_lognormal,
  dmeasure = meas_pois_interval,
  T=t,
  Y = d4$y,
  statenames = "x",
  filter.mean = TRUE,
  
)

plot(pf)


### Make comparison plot



