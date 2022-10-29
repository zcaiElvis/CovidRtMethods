library("pomp")
library("truncnorm")


tran_normal_unconst <- discrete_time(
  function(x, sigma1, ...){
    xnext <- rnorm(1, x, sd = sigma1)
    c(x = xnext)
  },
  delta.t = 1
)


tran_normal_const <- discrete_time(
  function(x, sigma1, ...){
    xnext <- rtruncnorm(1, a = 0, b = 2, mean = x, sd = sigma1)
    c(x = xnext)
  },
  delta.t = 1
)

tran_normal_epifilter <- discrete_time(
  function(x, sigma1, ...){
    xnext <- rnorm(1, x, sd = 0.1*sqrt(abs(x)))
    c(x=xnext)
  },
  delta.t = 1
)

### 1. Variance of log-normal changes with mu
### 2. R_t follows a lognormal distribution

tran_lognormal <- discrete_time(
  function(x, sdlog, ...){
    xnext <- rlnorm(1, meanlog = log(x), sdlog = sdlog)
    c(x = xnext)
  },
  delta.t = 1
)


ln_var <- function(logmu, logsd){
  return((exp(logsd^2)-1)*exp(2*logmu+logsd^2))
}

ln_sd <- function(logmu, logsd){
  return(sqrt(ln_var(logmu, logsd)))
}


