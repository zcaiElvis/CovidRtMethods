library("tidyverse")
library("dplyr")

# Function: Read from owid dataset
# Parameters:
  # country
  # data_loc
  # data_interval: interval of date to extract from, default to all data
get_owid_data <- function(country = "Canada", data_loc = "data/raw/owid_Sep5.csv",
                          data_interval = c()){
  data <- read.csv(data_loc, header=TRUE, sep=",")
  covid <- data %>%
    filter(location == country)%>%
    select(c(date, new_cases)) %>%
    rename(y = new_cases)
  
  if(length(data_interval) != 0){
    covid <- data_from_interval(covid, data_interval)
  }
  return(covid)
}

# Function: Helper function, select columns within the date interval
# Parameters:
  # covid: dataset
  # date_interval: interval of date to extract from
data_from_interval <- function(covid, date_interval){
  dat <- covid
  dat$date <- as.Date(as.character(dat$date))
  dat <- subset(dat, date >= date_interval[1] & date <= date_interval[2])
  return(dat)
}

# ca <- get_owid_data(data_interval = c("2021-09-10", "2021-09-12"))