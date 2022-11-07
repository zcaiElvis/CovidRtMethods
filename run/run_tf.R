library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(testthat)



source("model/tf/r/admm_solvers.R")

d <- read.csv("data/processed/d2.csv")


r <- gtf_solver(1000, y = d$y, x = d$iwt, k=0, para = 0.01)


plot(c(r$theta/d$iwt)[30:500])



