library("tidyverse")
library("dplyr")
source("constant/constant.R")
source("function/disc_gamma.R")
source("function/make_plot.R")
source("model/pls/ridge.R")


### Importing Rt estimation pacakges
library(EpiEstim)
source("run/run_EpiEstim.R")

source("model/pls/ridge.R")

source("run/run_pf_all.R")

### Importing data, no cyclic effect

d = read.csv("data/processed/d2.csv")
plot(d$y, type = "l")
plot(d$r, type = "l")


### Run methods

# EpiEstim
epiestim <- run_epiestim(d, p1 = sid_covid_mean, p2 = sid_covid_sd)
plot(epiestim$`Median(R)`)

# PLS
loocv <- CV(d$iwt, d$y)
ridge <- get_r(d$iwt, d$y, lambda = loocv$lambdas[1])
plot(ridge)
get_loss(d$iwt, d$y, ridge, loocv$lambdas[1])

# NLP GD
ridge_gd <- run_ridge_gd(d, loocv$lambdas[1])
plot(ridge_gd$estimate)
get_loss(d$iwt, d$y, ridge_gd$estimate, loocv$lambdas[1])

# Particle Filter
pf <- run_pf(d)
plot(pf$filter.mean.x)

# EpiNow
epin <- read.csv("data/results/epinow2_d2.csv")
plot(epin$median)
