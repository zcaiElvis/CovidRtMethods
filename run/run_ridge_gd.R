source("model/pls/penalties_smooth.R")



ridge_obj <- function(data, par, loss_func, iwt = iwt, smooth_func, penalties, pen_func = identity, ...){
  dat_length = nrow(data)
  loss = loss_func(z=data$y, iwt = iwt, r = par)
  r_pen <- penalties$r* smooth_func(par)
  obj_value = sum(loss+r_pen)
  return(obj_value)
}


run_ridge_gd <- function(data, lambda){
  ### Initialize
  dat_length = length(data$y)
  init_r = rep(1, dat_length)
  nlm(f=ridge_obj, p = init_r, iterlim =5000, print.level = 0, data=data, penalties = list("r"=lambda), 
      loss_func = normal_loss, smooth_func = r_smooth_penalty, iwt = data$iwt)
}

