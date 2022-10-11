
# Function: Calculate KL penalty
# Parameters:
#   z - reported incidence case
#   iwt - mean of poisson
#   dat_length- length of data
#   o - outlier, whether outlier is allowed to change in KL loss

kl_loss <- function(z, iwt, r, dat_length, o=c(0)){
  p = r*iwt + o
  l = rep(0, dat_length)
  for(i in 1:dat_length){
    if(z[i]>0 & p[i]>0){
      l[i] = z[i]*log(z[i]/p[i])+p[i]-z[i]
      
    }else if(z[i]==0 & p[i] >= 0){
      l[i] = p[i]
    }else{
      l[i] = 50000000 # this should be infinity
    }
  }
  
  return(sum(l))
}

# Function: Calculate smoothness penalty
# Parameters:
  # r - effective reproduction number
  # dat_length: length of data
r_smooth_loss <- function(r, dat_length){
  penalty = rep(0, dat_length-2)
  for(i in 2:(dat_length-1)){
    penalty[i] = abs(0.5*r[i-1]-r[i]+0.5*r[i+1])
  }
  return(sum(penalty))
}


# Function: Calculate outlier size penalty
# Parameters:
  # o - vector of outliers
o_size_loss <- function(o){
  return(sum(abs(o)))
}

# Function: Calculate negativity penalty
# Parameters:
  # r - effective reproduction number
pos_loss <- function(r, dat_length){
  pl <- rep(0, dat_length)
  for(j in 1:dat_length){
    if(r[j] < 0){
      pl[j] = 50000000 # this should be infinity
    }
  }
  return(sum(pl))
}



# Function: Combining all losses
# Parameters:
  # data - dataframe
    # column y - incidence
    # column iwt - sum(Iw)
  # par - matrix
    # column 1 - r
    # column 2 - o
  # penalty - list
    # p_r - penalty of r smoothness
    # p_o - penalty of o size
  # loss - vector of (0,1), controller to firing up kl, rl, ol, pl

all_loss <- function(data, par, dat_length, penalty, loss, model_o){
  par <- matrix(par, dat_length, 2)
  if(model_o){
    kl <- kl_loss(z= data$y, iwt = data$iwt, r = par[,1], dat_length =dat_length, o=par[,2])
  }else{
    kl <- kl_loss(z= data$y, iwt = data$iwt, r = par[,1], dat_length =dat_length)
  }
  rl <- penalty$p_r*r_smooth_loss(par[,1], dat_length)
  ol <- penalty$p_o*o_size_loss(par[,2])
  pl <- pos_loss(par[,1], dat_length)
  
  loss <- loss*c(kl, rl, ol, pl)
  
  return(sum(loss))
}
