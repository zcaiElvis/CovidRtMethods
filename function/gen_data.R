source("function/disc_gamma.R")
source("constant/constant.R")


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


### Corresponds to simulation 2A in Epifilter
shape_scale <- gamma_reparam(sid_ebola_mean, sid_ebola_sd)
r1 <- c(rep(2, 100), rep(0.5, 200))
y1 <- gen_with_rt(r1, shape_scale[1], shape_scale[2])
write.csv(data.frame(r=r1, y=y1, idx=length(y1)), row.names = FALSE, file = "data/processed/a.csv")
plot(y1)


### Corresponds to simulation 2B in Epifilter
r2 <- c(rep(4, 40), rep(0.6, 40), rep(2, 70), rep(0.2, 150))
y2 <- gen_with_rt(r2, shape_scale[1], shape_scale[2])
write.csv(data.frame(r=r2, y=y2, idx=length(y2)), row.names = FALSE, file = "data/processed/b.csv")
plot(y2)


