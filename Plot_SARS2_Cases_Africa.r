
#load packages
library(tidyverse)
library(readxl)
library(xlsx)



##### LOAD DATA ---------------------------------------------------------------
## Set working directory
setwd("~/choose/directory")


## Load data on reported cases
Reported_Cases <- read_excel("owid-covid-data.xlsx") %>%
# select only African data
  filter(continent == "Africa")

Cases_selected <- Reported_Cases %>%
# select needed collumns
  dplyr::select(3, 4, 7) %>%
  pivot_wider(names_from = location, values_from = new_cases_smoothed) 

# Replace NAs by 0
Cases_selected[is.na(Cases_selected)] <- 0

# Summarise cases
Cases_count <- Cases_selected %>% 
  rowwise() %>%
  mutate(Cases_total = sum(Algeria:Sudan, na.rm = TRUE)) %>%
  dplyr::select(1,58) %>%
  mutate(Date = as.Date(date))




##### PLOT DATA ---------------------------------------------------------------
ggplot(data = Cases_count, aes(x = Date, y = Cases_total))+
  geom_bar(stat = 'identity') + 
  theme_classic() +
  scale_x_date(breaks = seq.Date(from = as.Date("yyyy-mm-dd"), 
                                 to = as.Date("yyyy-mm-dd"), by = dd)) 



ggsave("Cases_Africa.pdf", width = 7, height = 4)