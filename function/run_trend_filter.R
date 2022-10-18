library(glmgen)
source("constant/constant.R")
source("function/get_iwt.R")



s1 <- read.csv("data/processed/a.csv")
s1_iwt <- get_iwt(s1$y, disc_gamma(1:nrow(s1), sid_ebola_shape, sid_ebola_scale))
y1 <- s1$y


out <- trendfilter(y = y1, x = s1_iwt, k = 1, family="poisson")
result <- predict(out, lambda=3.267028e+00, type="response")

plot(out$beta[,49])
# plot(result, type="l")




r_penalized <- as.data.frame(out$beta)

# colnames(r_penalized) <- 1:50

r_penalized$idx <- 1:300


ggplot(data= r_penalized, aes(x=idx))+
  geom_line(aes(x=idx, y = `3.015`), color = "darkred")+
  geom_line(aes(x=idx, y = `50.558`), color = "darkgreen")+
  geom_line(aes(x=idx, y = `1715.53`), color = "blue") +
  geom_line(aes(x=idx, y = `7024.95`), color = "purple") +
  theme_bw()