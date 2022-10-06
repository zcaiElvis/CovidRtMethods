library(extraDistr)

disc_gamma <- function(x, shape = 2, scale = 2){
  pgm <- pgamma(x, shape = shape, scale = scale)
  pgm <- c(0, pgm)
  pgm <- diff(pgm)
  pgm <- pgm/sum(pgm)
  return(pgm)
}

