library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(testthat)

# load ADMM solvers ----
source(here::here("model/tf/r/admm_solvers.R"))

# -------------------------------------------------------------------
# generate signals (no noise) ----
make_signal_v <- function(n) { 
  n <- as.integer(round(n))
  n2 <- n %/% 2
  r <- seq(1, n2 + 1, length.out = n2 + 1)
  l <- rev(r)
  r <- r[-1]
  if (n %% 2 == 0L) l <- l[-(n2 + 1)]
  return(c(l, r) / n) # change *scalar
}
n <- 10
y <- make_signal_v(n) #+ rnorm(n,0,.1)
plot(y)
x = rep(1,10)

# modeling ----
max_iter = 1e3
k=1 #k>=1

D = D_generator(n,k)
lam_max = max(solve(D %*% t(D)) %*% D %*% y); lam_max
para=0.1

mod_gaus = gtf_solver(max_iter=max_iter, y=y, x=x, k=k, para=para)
mod_pois = ptf_solver(1, y=y, x=x, k=k, para=para)

# visualization ----
res <- tibble(
  signal = 1:n,
  gaus_est = mod_gaus$theta[,1],
  pois_est = mod_pois$theta[,1]
)

cbPalette <- c("#E69F00", "#D55E00", "#009E73", "#0072B2")
fig <- res %>%
  pivot_longer(!signal, names_to = "Algorithm", values_to = "value") %>%
  group_by(Algorithm) %>%
  ggplot(aes(x = signal)) +
  geom_point(aes(y=rep(y,each=2)), size=1) +
  geom_line(aes(y = value, color=Algorithm, linetype=Algorithm), size=.7) +
  labs(x = "Signal", y = "Trend") +
  scale_colour_manual(values=cbPalette) +
  scale_x_continuous(breaks = c(1,4,7,10)) +
  theme_bw()
ggsave(here::here("fig.png"), plot=fig)
