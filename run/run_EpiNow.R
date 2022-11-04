library("EpiNow2")
source("constant/constant.R")
source("function/make_plot.R")


run_epinow <- function(data, p1, p2){
  data_size <- nrow(data)
  rdate <- data.frame(date =read.csv("data/processed/rdate.csv"))
  rdate <- rdate$x[1:data_size]
  
  generation_time = list(mean = p1, mean_sd = 0.01, sd = p2, sd_sd = 0.01, max=30)
  
  result_a <- epinow(reported_cases = data.frame(date = as.Date(rdate), confirm=data$y), 
                     generation_time = generation_time,
                     rt = rt_opts(prior=list(mean=1.5, sd=3)),
                     stan = stan_opts(cores=4, samples = 1000))

  return(result_a)
}


d <- read.csv("data/processed/d2.csv")

epinow <- run_epinow(d, sid_covid_mean, sid_covid_sd)

ggsave("epinow_d2.png", plot =plot(summary(epinow, type = "parameters", params = "R")$median),
       path = "plot/epinow")

write.csv(summary(epinow, type = "parameters", params = "R"), file = "data/results/epinow2_d2.csv")



