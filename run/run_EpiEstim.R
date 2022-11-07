library("tidyverse")
library("dplyr")
library("EpiEstim")

### Function: Run epiestim with default parameters
### p1, p2 are mean and sd of si. Here si is gamma, reparametrized to mean and sd
run_epiestim <- function(data, p1, p2){
  epiestim_result <- estimate_R(data$y, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = p1, 
                                  std_si = p2)))
  
  median_r <- epiestim_result$R$`Median(R)`
  
  median_r <- c(rep(mean(median_r[1:10]),7), median_r)
  
  return(median_r)
}


if (sys.nframe() == 0) {
  sim_a <- read.csv("data/processed/a.csv")
  sim_b <- read.csv("data/processed/b.csv")
  
  
  true_r1 <- sim_a$r
  pred_r1 <- run_epiestim(sim_a, p1=sid_ebola_mean, p2=sid_ebola_sd)
  compare1 <- compare_rt(true_r1, pred_r1)
  ggsave("a.png",plot = compare1, path = "plot/epiestim")
  
  true_r2 <- sim_b$r
  pred_r2 <- run_epiestim(sim_b, p1=sid_ebola_mean, p2=sid_ebola_sd)
  compare2 <- compare_rt(true_r2, pred_r2)
  ggsave("b.png",plot = compare2, path = "plot/epiestim")
  
}



