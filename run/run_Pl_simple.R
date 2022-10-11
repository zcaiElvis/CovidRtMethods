source("model/pls/penalties_smooth.R")
source("function/disc_gamma.R")
source("function/get_iwt.R")
source("constant/constant.R")
source("function/make_plot.R")




run_pl <- function(data, p1, p2, loss_func){
  dat_length = nrow(data)
  maxit = 100
  data$iwt <- get_iwt(data$y, disc_gamma(x=1:dat_length, shape = p1, scale = p2))
  init_r <- rep(1, dat_length)
  penalty <- list("rl" = 20)
  
  result <- nlm(f=smooth_loss, p = init_r, iterlim =2000, print.level = 1, data=data, penalties = penalty, loss_func = rmse_loss)
  return(result$estimate)
  
}



sim_a <- read.csv("data/processed/a.csv")
sim_b <- read.csv("data/processed/b.csv")


true_r1 <- sim_a$r
pred_r1 <- run_pl(sim_a, sid_ebola_shape, sid_ebola_scale, loss_func = pois_loss)
compare1 <- compare_rt(true_r1, pred_r1)
ggsave("a_pois.png",plot = compare1, path = "plot/pls")


true_r2 <- sim_b$r
pred_r2 <- run_pl(sim_b, sid_ebola_shape, sid_ebola_scale, loss_func = pois_loss)
compare2 <- compare_rt(true_r2, pred_r2)
ggsave("b_pois.png",plot = compare2, path = "plot/pls")



