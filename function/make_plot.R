library("tidyverse")
library("dplyr")


### Function: Plotting two Rt

compare_rt <- function(r1, r2){
  rsize = length(r1)
  if(length(r1) != length(r2)){
    rsize = min(length(r1),length(r2))
    r1 <- tail(r1, rsize)
    r2 <- tail(r2, rsize)
  }

  plot<- ggplot(data.frame(idx=1:rsize, r1=r1, r2=r2), aes(x=idx))+
    geom_line(aes(x=idx, y = r1), color = "blue")+
    geom_line(aes(x=idx, y = r2), color = "red")+
    theme_bw()+
    ylab("R")+
    xlab("Time")
  
  return(plot)
}



dfaom <- compare_rt(1:10, 4:13)