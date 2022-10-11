source("function/process_data.R")

### Save date
save_date <- function(){
  data <- get_owid_data()
  date <- data$date
  write.csv(date, "data/processed/rdate.csv", row.names = FALSE)
}
# save_date()
