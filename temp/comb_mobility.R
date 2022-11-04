library("tidyverse")
library("dplyr")


m2020 <- read.csv("data/raw/2020_CA_Region_Mobility_Report.csv")
m2021 <- read.csv("data/raw/2021_CA_Region_Mobility_Report.csv")
m2022 <- read.csv("data/raw/2022_CA_Region_Mobility_Report.csv")


m2020 %<>%
  filter(sub_region_1 == "") %>%
  filter(sub_region_2 == "")

m2021 %<>%
  filter(sub_region_1 == "") %>%
  filter(sub_region_2 == "")

m2022 %<>%
  filter(sub_region_1 == "") %>%
  filter(sub_region_2 == "")

m <- rbind(m2020, m2021, m2022)

m <- select(m, date, retail_and_recreation_percent_change_from_baseline, 
         grocery_and_pharmacy_percent_change_from_baseline,
         parks_percent_change_from_baseline,
         transit_stations_percent_change_from_baseline,
         workplaces_percent_change_from_baseline,
         residential_percent_change_from_baseline)

