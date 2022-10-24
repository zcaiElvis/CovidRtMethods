source("function/get_iwt.R")
source("function/disc_gamma.R")
source("constant/constant.R")
library("glmgen")
### Import s1
s1 <- read.csv("data/processed/a.csv")
s1_iwt <- get_iwt(s1$y, disc_gamma(1:nrow(s1), sid_ebola_shape, sid_ebola_scale))



tf <- glmgen::trendfilter(x = s1_iwt, y = s1$y, k = 0, family="gaussian")

lambdas = summary(tf)[,2]
rss = summary(tf)[,3]
optim_lambda = lambdas[order(rss)[1]]

betas <- data.frame(tf$beta)

plot(tf$beta[,15])