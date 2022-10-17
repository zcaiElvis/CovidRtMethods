source("model/pls/penalties_smooth.R")
source("function/disc_gamma.R")
source("function/get_iwt.R")
source("constant/constant.R")
source("function/make_plot.R")




run_pl <- function(data, p1, p2, loss_func, penalty = list("rl"=50)){
  dat_length = nrow(data)
  maxit = 100
  data$iwt <- get_iwt(data$y, disc_gamma(x=1:dat_length, shape = p1, scale = p2))
  init_r <- rep(1, dat_length)
  
  result <- nlm(f=smooth_loss, p = init_r, iterlim =2000, print.level = 0, data=data, penalties = penalty, loss_func = rmse_loss)
  return(result)
  
}


if (sys.nframe() == 0) {
  sim_a <- read.csv("data/processed/a.csv")
  iwt_a <- get_iwt(sim_a$y, disc_gamma(x=1:nrow(sim_a), shape = sid_ebola_shape, scale = sid_ebola_scale))
  pred_a <- run_pl(sim_a, sid_ebola_shape, sid_ebola_scale, penalty=list("rl" = 20))
  diag_a <- diag_plots(true_r = sim_a$r, pred_r = pred_a$estimate, iwt = iwt_a, i = sim_a$y)
  
  diag_a$rt
  ggsave("a_pois.png", plot = diag_a$rt, path="plot/pls")
  diag_a$oneday
  ggsave("a_1day_pois.png", plot=diag_a$oneday, path = "plot/pls")
  
  
  
  sim_b <- read.csv("data/processed/b.csv")
  iwt_b <- get_iwt(sim_b$y, disc_gamma(x=1:nrow(sim_b), shape = sid_ebola_shape, scale = sid_ebola_scale))
  pred_b <- run_pl(sim_b, sid_ebola_shape, sid_ebola_scale, penalty = list("rl" = 20))
  diag_b <- diag_plots(true_r = sim_b$r, pred_r = pred_b$estimate, iwt = iwt_b, i = sim_b$y)
  
  diag_b$rt
  ggsave("b_pois.png", plot = diag_b$rt, path="plot/pls")
  diag_b$oneday
  ggsave("b_1day_pois.png", plot=diag_b$oneday, path = "plot/pls")
  
  
}





