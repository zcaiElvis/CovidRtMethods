# source("function/disc_gamma.R")
# source("constant/constant.R")
library("zoo")


gen_with_rt <- function(rt, init_i = 1, g_shape = 2, g_scale = 2, cyc_amp = 1){
  rt_len <- length(rt)
  i <- rep(0, rt_len)
  i[1] <- init_i
  for(j in 2:rt_len){
    lamb <- rt[j]*sum(i[(j-1):1]*disc_gamma(1:(j-1), shape = g_shape, scale = g_scale))
    i[j] <- rpois(1, lamb)
    if(cyc_amp != 0){i[j] <- ceiling(max(0, (i[j]/cyc_amp)*sin(7*j)+i[j]))}
  }
  
  return(i)
}


smooth_rt <- function(rt, times){
  for(i in 1:times){
    rt <- rollmean(rt,10, fill=NA, align = "left")
    rt <- na.locf(rt)
  }
  return(rt)
}


if (sys.nframe() == 0) {
  plot(smooth_rt(c(rep(1,20), rep(1.5, 20), rep(3, 10))))
  smooth_rt(c(rep(1,20), rep(1.5, 20), rep(3, 10)))
  plot(gen_with_rt(rep(1, 1000), init_i=10, cyc_amp=4), type="l")
}


### Add mobility because it is confounding. R should be the contact rate in a mixing condition




