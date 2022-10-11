get_iwt <- function(I, w){
  t <- length(I)
  iwt <- rep(0, t)
  iwt[1] <- I[1]
  
  for(a in 2:t){
    iwt[a] <- sum(I[(a-1):1]*w[1:(a-1)])
  }
  return(iwt)
}
