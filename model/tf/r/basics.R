# divided difference matrices of dim (n-2)*n:
D_generator <- function(n, k){
  
  D = matrix(0, nrow=n-1, ncol=n)
  for (i in 1:(n-1)){
    D[i, i:(i+1)] = matrix(c(-1,1), 1, 2)
  }
  
  if (k==0){
    return(D)
  } else if (k>0){
    for(i in 1:k){
      Di = matrix(0, nrow=n-i-1, ncol=n-i)
      for (d in 1:(n-i-1)){
        Di[d, d:(d+1)] = matrix(c(-1,1), 1, 2)
      }
      D = Di %*% D
    }
  }
  
  return(D)
}
