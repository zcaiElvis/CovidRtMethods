library(glmgen)
source("constant/constant.R")
source("function/get_iwt.R")



s1 <- read.csv("data/processed/a.csv")
s1_iwt <- get_iwt(s1$y, disc_gamma(1:nrow(s1), sid_ebola_shape, sid_ebola_scale))
y1 <- s1$y


out <- trendfilter(y = y1, x = s1_iwt, k = 1, family="poisson")
result <- predict(out, lambda=3.267028e+00, type="response")

plot(out$beta[,30])
# plot(result, type="l")
