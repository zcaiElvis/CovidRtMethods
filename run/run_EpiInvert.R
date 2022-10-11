library("EpiInvert")



### Dont run EpiInvert for now, crashes


setwd("~/Desktop/School/research/covid_rt_comparison/")


sim_a <- read.csv("data/processed/a.csv")
sim_b <- read.csv("data/processed/b.csv")

epiinvert_result <- EpiInvert(sim_a$y,last_incidence_date = "2022-04-01", "1000-01-01",
                              select_params(list(
                                max_time_interval = 300)))
