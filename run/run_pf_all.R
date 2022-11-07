library("pomp")
# source("constant/constant.R")
# source("function/disc_gamma.R")
# source("function/make_plot.R")


### Transition equation

tran <- discrete_time(
  function(x, s, sdlog, ...){
    snext <- rlnorm(1, meanlog = log(s), sdlog = sdlog)
    xnext <- rlnorm(1, meanlog = log(x), sdlog = snext)
    c(x = xnext, s = snext)
  },
  delta.t = 1
)



### Observation equation

meas_pois_interval <- function(t, x, s, T, Y, y, g_shape, g_scale, ..., log){
  if(t > 1){
    lambda = x[1] * sum(Y[(t-1):1] * disc_gamma(1:(t-1), shape = g_shape, scale = g_scale))
    y*log(lambda)-lambda
  }else{
    dpois(round(y), lambda = 1)
  }
}


### Run particle filter

run_pf <- function(d, sdlog){
  init_vals = c(sdlog = sdlog,
                g_shape = sid_covid_shape, g_scale = sid_covid_scale)
  
  pf <- pfilter(
    Np=1000,
    times = "idx",
    t0 = 1,
    data = d,
    params = init_vals,
    rinit = function(x0, ...){c(x=1.5, s = 0.1)},
    rprocess = tran,
    dmeasure = meas_pois_interval,
    T=t,
    Y = d$y,
    statenames = c("x"),
    filter.mean = TRUE,
    filter.traj = TRUE,
    # pred.mean = TRUE,
    # pred.var = TRUE,
  )
  
  # result = as.data.frame(pf)
  
  return(result)
}


source("function/gen_syn_data.R")
source("constant/constant.R")
souce("function/disc_gamma.R")
d <- read.csv("data/processed/d2.csv")

result <- run_pf(d, 0.1)

meaniwt <- result$filter.mean.x*disc_gamma(c())





### Synthetic data import

# d4<- read.csv("data/processed/d4.csv")
# plot(d4$y)
# plot(d4$r)

### Initial values

# init_vals = c(sdlog = 0.1,
#               g_shape = sid_covid_shape, g_scale = sid_covid_scale)

### Run particle filter default

# pf <- pfilter(
#   Np=1000,
#   times = "idx",
#   t0 = 1,
#   data = d4,
#   params = init_vals,
#   rinit = function(x0, ...){c(x=1, s = 0.05)},
#   rprocess = tran,
#   dmeasure = meas_pois_interval,
#   T=t,
#   Y = d4$y,
#   statenames = c("x", "s"),
#   filter.mean = TRUE,
#   filter.traj = TRUE,
#   # pred.mean = TRUE,
#   # pred.var = TRUE,
# )

# plot(pf)
# 
# result = as.data.frame(pf)
# 
# 
# ### Make comparison plot
# diag_pf = diag_plots(d4$r, result$r, d4$iwt, d4$y, cap=0)
# 
# diag_pf$rt
# 
# diag_pf$oneday
# 
# ggplot(data=data.frame(idx = d4$idx, true= d4$r, pred = result$filter.mean.x), aes(x=idx))+
#   geom_line(aes(x = idx, y = true), color="blue")+
#   geom_line(aes(x=idx, y = pred), color = "red")+
#   theme_bw()
# 
# 
# plot(result$filter.mean.s, type = "l")


