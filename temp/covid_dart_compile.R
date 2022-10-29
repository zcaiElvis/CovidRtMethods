source("function/process_data.R")
library(zoo)


ca <- get_owid_data()
ca$y <- na.locf(ca$y)

write.csv(ca, file="data/processed/ca.csv", row.names=FALSE)
