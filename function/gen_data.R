source("function/disc_gamma.R")


### Function: Generating Incidence given Rt.
### g_shape, g_scale are gamma shape and scale params

gen_with_rt <- function(rt, g_shape = 2, g_scale = 2){
  rt_len <- length(rt)
  i <- rep(0, rt_len)
  i[1] <- 1
  for(j in 2:rt_len){
    lamb <- rt[j]*sum(i[(j-1):1]*disc_gamma(1:(j-1), shape = g_shape, scale = g_scale))
    i[j] <- rpois(1, lamb)
  }
  return(i)
}

### Function: Reparametrize Gamma mean/sd to shape/scale
gamma_reparam <- function(g_mean, g_sd){
  g_shape = (g_mean/g_sd)^2
  g_scale = 1/(g_mean/g_sd^2)
  return(c(g_shape, g_scale))
}



shape_scale <- gamma_reparam(15.3, 9.3)
r1 <- c(rep(2, 100), rep(0.5, 200))

plot(gen_with_rt(r1, shape_scale[1], shape_scale[2]))