library("EpiNow2")
source("constant/constant.R")
source("function/make_plot.R")


run_epinow <- function(data, p1, p2){
  data_size <- nrow(data)
  rdate <- data.frame(date =read.csv("data/processed/rdate.csv"))
  rdate <- rdate$x[1:data_size]
  
  generation_time = list(mean = p1, mean_sd = 0.1, sd = p2, sd_sd = 0.1, max=30)
  result_a <- epinow(reported_cases = data.frame(date = as.Date(rdate), confirm=data$y), generation_time = generation_time,
                     rt = rt_opts(prior=list(mean=1.5, sd=1)),
                     stan = stan_opts(cores=2, samples = 1000))

  return(result_a$estimates$summarised$mean)
}


sim_a <- read.csv("data/processed/a.csv")
sim_b <- read.csv("data/processed/b.csv")


true_r1 <- sim_a$r
pred_r1 <- run_epinow(sim_a, sid_ebola_mean, sid_ebola_sd)
compare1 <- compare_rt(true_r1, pred_r1)
ggsave("a.png",plot = compare1, path = "plot/epinow")

true_r2 <- sim_b$r
pred_r2 <- run_epinow(sim_b, sid_ebola_mean, sid_ebola_sd)
compare2 <- compare_rt(true_r2, pred_r2)
ggsave("b.png",plot=compare2, path="plot/epinow")




