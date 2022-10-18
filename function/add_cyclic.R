





### Function: generate y based on Rt, with added cyclic effect
gen_with_rt_cycle <- function(rt, g_shape = 5, g_scale = 5){
  rt_len <- length(rt)
  i <- rep(0, rt_len)
  i[1] <- 1
  for(j in 2:rt_len){
    lamb <- rt[j]*sum(i[(j-1):1]*disc_gamma(1:(j-1), shape = g_shape, scale = g_scale))
    i[j] <- rpois(1, lamb)
    i[j] <- ceiling(max(0, i[j]/5*sin(7*j)+i[j]))
  }

  return(i)
}



if (sys.nframe() == 0) {
  # shape_scale <- gamma_reparam(sid_ebola_mean, sid_ebola_sd)
  # r3 <- c(rep(2, 100), rep(0.5, 200))
  # y3 <- gen_with_rt_cycle(r1, shape_scale[1], shape_scale[2])
  # write.csv(data.frame(r=r3, y=y3, idx=length(y3)), row.names = FALSE, file = "data/processed/c.csv")
  # plot(y3, type="l")
  # 
  # 
  shape_scale <- gamma_reparam(sid_ebola_mean, sid_ebola_sd)
  r4 <- c(rep(2, 100), rep(0.5, 200))
  y4 <- gen_with_rt_cycle(r2, shape_scale[1], shape_scale[2])
  write.csv(data.frame(r=r4, y=y4, idx=length(y4)), row.names = FALSE, file = "data/processed/d.csv")
  plot(y4, type="l")
  plot(r4)
  
  # shape_scale <- gamma_reparam(sid_ebola_mean, sid_ebola_sd)
  # r5 <- c(rep(1.2, 100), rep(2, 100), rep(0.8, 100), rep(1, 100), rep(1.5, 100), rep(0.6, 100))
  # y5 <- gen_with_rt_cycle(r5, shape_scale[1], shape_scale[2])
  # write.csv(data.frame(r=r5, y=y5, idx=length(y5)), row.names = FALSE, file = "data/processed/e.csv")
  # plot(y5, type="l")
  # 
}


