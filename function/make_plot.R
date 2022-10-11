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


### Function: Plot one-day ahead prediction

onedayhead <- function(r, iwt, i){
  dat_length <- length(r)
  pred <- r*iwt
  plot<- ggplot(data.frame(idx=1:dat_length, pred = pred, i = i), aes(x=idx))+
    geom_line(aes(x=idx, y = i), color = "blue")+
    geom_line(aes(x=idx, y = pred), color = "red")+
    theme_bw()+
    ylab("Case")+
    xlab("Time")
  return(plot)
}



diag_plots <- function(true_r, pred_r, iwt, i){
  output = list()
  compare_r <- compare_rt(true_r, pred_r)
  oneday_r <- onedayhead(pred_r, iwt, i)
  
  output$rt = compare_r
  output$oneday = oneday_r
  return(output)
}


# 
# dfaom <- compare_rt(1:10, 4:13)