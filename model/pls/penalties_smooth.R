library(Metrics)


# Function: KL loss of Poisson and true cases
# Parameters:
  # z - incidence case
  # iwt - poisson mean
  # r - effective reproduction number
  # dat_length - length of data
normal_loss <- function(z, iwt, r){
  return(sum((z-r*iwt)^2))
}


pois_loss <- function(z, iwt, r){
  p = r*iwt
  # loss = dpois(z, p)
  loss = dpois(z, p, log=TRUE)
  return(-sum(loss))
}


# Function: Smoothness loss of R
# Parameters:
  # r - effective reproduction number
  # dat_length: length of data

r_smooth_penalty <- function(r){
  return(sum(diff(r)^2))
}

r_smooth_penalty_varies <- function(r, i, penalty){
  return(penalty*(i/10)*sum(diff(r)^2))
}


pos_penalty <- function(r){
  return(as.numeric(any(r<0))*1000)
}


# Function: Combining all losses of smooth loss functions
# Parameters:
  # data - dataframe
    # column y - incidence
    # column iwt - sum(Iw)
  # par - vector of r
  # penalties - vector of lambda penalty

smooth_loss <- function(data, par, penalties, loss_func){
  dat_length <- length(data)
  
  kl <- loss_func(z = data$y, iwt = data$iwt, r = par)
  rl <- penalties$rl* r_smooth_penalty(par)
  all_losses = kl+rl
  
  return(sum(all_losses))
}


# Function: Combining all losses of smooth loss functions
# Parameters:
# data - dataframe
# column y - incidence
# column iwt - sum(Iw)
# par - vector of r
# penalties - vector of lambda penalty

smooth_loss <- function(data, par, penalties, loss_func, pen_change = FALSE){
  dat_length <- length(data)
  
  kl <- loss_func(z = data$y, iwt = data$iwt, r = par)
  if(pen_change){
    rl <- r_smooth_penalty_varies(par, data$y, penalties$rl)
  }
  else{
    rl <- penalties$rl* r_smooth_penalty(par)
  }
  all_losses = kl+rl
  
  return(sum(all_losses))
}



