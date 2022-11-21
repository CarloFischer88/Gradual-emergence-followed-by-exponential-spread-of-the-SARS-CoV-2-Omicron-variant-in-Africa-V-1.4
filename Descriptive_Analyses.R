#load packages
library(tidyverse)
library(readxl)
library(xlsx)
library(scico)
library(viridis)
library(rcompanion)
library(lubridate)
library(ciTools)
library(stats)
library(data.table)

Sys.setenv(TZ = "Europe/London")



#### LOAD DATA ----------------------------------------------------------------
## Set working directory
setwd("~/working/directory")


## Load dataset
data <- read.csv("PCR_data.csv",
                 header = TRUE,        # Whether to read the header or not
                 sep = ",",            # Separator of the values
                 dec = ".")
## Change language settings
Encoding(data$Location) <- "UTF-8"




#### TIDY DATA ----------------------------------------------------------------
data_tidy <- data %>%
  # Change ages to completed life years
  mutate(Age = floor(as.numeric(Age))) %>% 
  # Remove SARS-CoV-2 negative samples (= NA in Variant collumn)
  filter(!is.na(Variant)) %>%
  # Change date and create binned dates
  mutate(Collection_date2 = as.POSIXct(Collection_date, format = "%Y-%m-%d")) 
# Filter data
data_tidy <- data_tidy %>%
  # Remove samples collected before study timeframe
  filter(
    Collection_date2 >= as.POSIXct("2021-06-01", format = "%Y-%m-%d"))
## Add sampling month (year-month) for each sample
data_tidy$Month <- format(data_tidy$Collection_date2, format="%Y-%m")

#### Identify first BA.1 marker positive samples ------------------------------
First_Omicron <- as.data.frame(data_tidy %>% 
                                 filter(Omicron == 1) %>%
                                 group_by(Country) %>%
                                 summarise(min(Collection_date2))) %>%
  rename(Collection_date = 2) 





#### CALCULATE MEDIAN AGE -----------------------------------------------------
Median_Age <- as.data.frame.matrix(groupwiseMedian(
  Age ~ Country,
  data = data_tidy %>% filter(!is.na(Age)),
  conf = 0.95,
  bca = FALSE,
  R = 1000,
  exact = TRUE,
  digits = 2)) %>%
  rename(Lower = 5, Upper = 6) %>%
  mutate('Median (95% CI)' = 
           paste(Median, ' (', Lower, ' - ', Upper, ')', sep = "")) %>%
  dplyr::select(1, 2, 7)

# Export table 
write.csv(Median_Age,
          "~/directory\\Median_Age.csv")



#### CALCULATE SEX FRACTIONS --------------------------------------------------
Sex <- as.data.frame(table(data_tidy$Country,
                           data_tidy$Sex,
                           useNA = "always")) %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  rename(Country = 1,
         Female = 2,
         Male = 3,
         Unknown = 4) %>%
  mutate(Total = Female + Male + Unknown) %>%
  filter(!is.na(Country)) %>%
  mutate(Female = round((Female/Total)*100, digits = 1),
         Male = round((Male/Total)*100, digits = 1),
         Unknown = round((Unknown/Total)*100, digits = 1)) %>%
  dplyr::select(-Total)

# Export table 
write.csv(Sex,
          "~/directory\\Sex.csv")




#### PLOT AVAILABLE SAMPLES PER MONTH BY COUNTRY ------------------------------
## Count samples
data_count <- data_tidy %>%
  group_by(Country, Month) %>%
  summarise(Count = n()) %>%
  filter(Month >= "2021-06")

## Plot sample number by month
ggplot(data = data_count, 
               aes(x = Month, 
                   y = factor(Country,
                              # Reorder y axis
                              levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = log10(Count))) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#E9EFF5", 
                       mid = "#F5CFBA", 
                       high = "#610D45", 
                       midpoint = 1.5,
                       labels = c("1", "10", "100", "1,000")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "Samples")


# Save plot
ggsave("Samples_Heatmap.pdf", height = 4, width = 7)





#### CALCULATE MONTHLY OMICRON AND DELTA FRACTION -----------------------------
data_fraction <- data_tidy %>%
  mutate(Omicron_dummy = ifelse(Variant == "Omicron", 1, 0),
         Delta_dummy = ifelse(Variant == "Delta", 1, 0)) %>%
  group_by(Country, Month) %>%
  summarise(Omicron_fraction = mean(Omicron_dummy),
            Delta_fraction = mean(Delta_dummy))



### Plot Omicron fraction
ggplot(data = data_fraction, 
       aes(x = Month, 
           y = factor(Country,
       # Reorder y axis
           levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = Omicron_fraction)) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#DAE8E8", 
                       mid = "#9788B3", 
                       high = "#5E0523", 
                       midpoint = .6,
                       #labels = c("1", "10", "100", "1,000")
  ) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "Omicron fraction")

# Save plot
ggsave("Omicron_Heatmap.pdf", height = 4, width = 7)




### Plot Delta fraction
ggplot(data = data_fraction, 
       aes(x = Month, 
           y = factor(Country,
                      # Reorder y axis
                      levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = Delta_fraction)) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#DAE8E8", 
                       mid = "#9788B3", 
                       high = "#5E0523", 
                       midpoint = .6,
                       #labels = c("1", "10", "100", "1,000")
  ) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "Delta fraction")

# Save plot
ggsave("Delta_Heatmap.pdf", height = 4, width = 7)










