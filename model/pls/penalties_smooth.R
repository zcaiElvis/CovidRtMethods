library(Metrics)
# Function: KL loss of Poisson and true cases
# Parameters:
  # z - incidence case
  # iwt - poisson mean
  # r - effective reproduction number
  # dat_length - length of data
rmse_loss <- function(z, iwt, r){
  p = r*iwt
  loss = rmse(z, p)
  return(loss)
}


pois_loss <- function(z, iwt, r){
  p = r*iwt
  loss = dpois(z, p)
  return(sum(loss))
}


# Function: Smoothness loss of R
# Parameters:
  # r - effective reproduction number
  # dat_length: length of data

# r_smooth_loss <- function(r, dat_length){
#   penalty = rep(0, dat_length-2)
#   for(i in 2:(dat_length-1)){
#     penalty[i] = abs((0.5*r[i-1]-r[i]+0.5*r[i+1])) # this was a L1 norm in the paper
#   }
#   return(sum(penalty))
# }

r_smooth_loss <- function(r, dat_length){
  penalty = rep(0, dat_length-1)
  for(j in 2:dat_length){
    penalty[j] = (r[j]-r[j-1])^2
  }
  return(sum(penalty))
}


pos_smooth_loss <- function(r, dat_length){
  pos_penalty = 0
  for(k in 1:dat_length){
    if(r[k] < 0){
      pos_penalty = pos_penalty + 1000
    }
  }
  return(pos_penalty)
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
  rl <- penalties$rl* r_smooth_loss(par, dat_length = dat_length)
  all_losses = kl+rl
  
  return(sum(all_losses))
}


