### ### ### ### Predict R0 of BA.1 among African countries ### ### ### ###

# NOTE: BA.1 R0/RT, BA.1 cases, Delta cases, and total cases are calcualte/shown
# relative to the first BA.1 signal for each country

#load packages
library(tidyverse)
library(readxl)
library(xlsx)
library(EpiEstim)
library(incidence2)
library(data.table)
library(matrixStats)



##### LOAD DATA ---------------------------------------------------------------
## Set working directory
setwd("working/directory")

## Load results on PCR-based BA.1 testing which already include data on first
## BA.1 cases acceptable for each country
data <- read.csv("PCR_data.csv",
                 header = TRUE,
                 sep = ",",
                 dec = ".")

## Define countries represented in the PCR data
## (Gabon is excluded as no BA.1 cases were observed in the sample)
Countries_Africa <- c("Algeria", "Angola", "Benin", "Burkina Faso", "Cameroon", "Ethiopia",
                      "Gambia", "Ghana", "Kenya", "Madagascar", "Mali", "Morocco", "Mozambique",
                      "Namibia", "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe")
## Load Our world in data dataset on reported SARS-CoV-2 cases
Reported_Cases <- read.csv("owid-covid-data.csv") %>%
  # Rename Country variable
  rename(Country = location) %>%
  # Remove data on countries not represented in the PCR data
  filter(Country %in% Countries_Africa)



##### ESTIMATE OMICRON FRACTION OF REPORTED SARS-CoV-2 CASES ------------------
### DEFINE FUNCTION
MODEL_OMICRON <- function(dataset, CountryX, FirstCase) {
  # Select Omicron data
  data_loaded <- dataset %>%
    filter(Country == CountryX) %>%
    # Calculate difference between first BA.1 signal and collection date
    mutate(Diff_1st_Case = 
             as.numeric(difftime((as.POSIXct(Collection_date, "%Y-%m-%d")),
                                 as.POSIXct(FirstCase, "%Y-%m-%d"),
                                 units= "days")))
  # Select data on reported cases
  data_cases <- Reported_Cases %>%
    filter(Country == CountryX) %>%
    # Select collumns with needed data
    select(1,2,6,8) %>% 
    # Transform dates into posixct format
    mutate(Collection_date = as.POSIXct(date, "%Y-%m-%d")) %>%
    mutate(Start_date = (as.POSIXct(FirstCase, "%Y-%m-%d"))) %>%
    # Calculate difference between first BA.1 signal and reporting date 
    mutate(Diff_1st_Case = as.numeric(difftime(Collection_date, Start_date, units= "days"))) %>%
    # Remove samples which were calculated more than 50 days before the first BA.1 signal
    filter(Diff_1st_Case >= -50) %>%
    # Select needed collumns
    select(Diff_1st_Case, new_cases_smoothed_per_million)
  # Calculate model
  glmO = glm(Omicron ~ Diff_1st_Case, data = data_loaded, 
             family = binomial(link = "logit"))
  # Create dummy dataset 
  newdata = data.frame(Diff_1st_Case = 
                         as.numeric(c(min(data_loaded$Diff_1st_Case, na.rm = TRUE):
                                        200)))
  # Predict BA.1 fraction over time
  predicted <- as.data.frame(predict(glmO, newdata, 
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]      
  predicted <- predicted %>% 
    rename(Diff_1st_Case = 1) %>%
    mutate(Diff_1st_Case = as.numeric(Diff_1st_Case)) %>%
    # Correct day difference as predicted begins with 1 instead of e.g. -50
    mutate(Diff_1st_Case= Diff_1st_Case + min(data_loaded$Diff_1st_Case) -1
    )
  # Merge data on cases and BA.1 fraction
  MERGED <- merge(predicted, data_cases, by = "Diff_1st_Case") %>%
    # Calculate BA.1 cases by multiplying BA.1 fraction with reported cases
    mutate(Omicron_Cases = round((Value*new_cases_smoothed_per_million
    ), 
    digits = 3))
  # Return dataframe
  return(MERGED)
}




##### APPLY FUNCTION TO ESTIMATE BA.1 CASES -----------------------------------
# Algeria
CALC_ALG <- MODEL_OMICRON(dataset = data,
                        CountryX = "Algeria",
                        FirstCase = "yyyy-mm-dd")
# Plot values for quality control
ggplot(data = CALC_ALG, aes(x = Diff_1st_Case, y = Value))+
  geom_point()

# Angola
CALC_ANG <- MODEL_OMICRON(dataset = data,
                        CountryX = "Angola",
                        FirstCase = "yyyy-mm-dd")
# Plot values for quality control
ggplot(data = CALC_ANG, aes(x = Diff_1st_Case, y = Value))+
  geom_point()

# ... (repreat for every included country)


#### Merge predicted values for BA.1 ------------------------------------------
## Join predicted country data in one list to join them in a data frame
Predicted_list <- list(CALC_ALG,
                       CALC_ANG,
                       CALC_BEN,
                       CALC_BUF,
                       CALC_CAM,
                       CALC_ETP,
                       CALC_GAM,
                       CALC_GHA,
                       CALC_KEN,
                       CALC_MAL,
                       CALC_MAR,
                       CALC_NMB,
                       CALC_NIG,
                       CALC_SEN,
                       CALC_SAF,
                       CALC_TOG,
                       CALC_UGA,
                       CALC_MOZ,
                       CALC_COG,
                       CALC_MAD,
                       CALC_ZIM)

## Merge all data frames in the list for predicted BA.1 cases
Predicted_merged <- Predicted_list %>% reduce(full_join, by='Diff_1st_Case') %>%
  # Select needed data
  select(-contains("smoothed"),
         -contains("Value"))


## Calculate median of predicted BA.1 cases over time 
Predicted_merged$Median <- round(rowMedians(as.matrix(Predicted_merged[,2:ncol(Predicted_merged)]), na.rm = TRUE), digits = 3)
## Calculate mean of predicted BA.1 cases over time
Predicted_merged$Mean <- round(rowMeans(Predicted_merged[,2:(ncol(Predicted_merged)-1)], na.rm = TRUE), digits = 3)



## Calculate mean total SARS-CoV-2 cases
Total_merged <- Predicted_list %>% reduce(full_join, by='Diff_1st_Case') %>%
  select(-contains("Value"),
         -contains("Omicron")) %>%
  rename(Time = 1) 
# Calculate mean
Total_merged$Mean <- round(rowMeans(Total_merged[,2:ncol(Total_merged)], na.rm = TRUE), digits = 3)





##### ESTIMATE DELTA FRACTION OF REPORTED SARS-CoV-2 CASES --------------------
### DEFINE FUNCTION
MODEL_DELTA <- function(dataset, CountryX, FirstCase) {
  # Select Omicorn data
  data_loaded <- dataset %>%
    # Select country-specific data
    filter(Country == CountryX) %>%
    # Calculate difference between first BA.1 signal and collection date
    mutate(Diff_1st_Case = 
             as.numeric(difftime((as.POSIXct(Collection_date, "%Y-%m-%d")),
                                 as.POSIXct(FirstCase, "%Y-%m-%d"),
                                 units= "days")))
  # Select data on reported SARS-CoV-2 cases
  data_cases <- Reported_Cases %>%
    # Select country-specific data
    filter(Country == CountryX) %>%
    # Select needed collumns
    select(1,2,6,8) %>% 
    # Transform dates into posixct format
    mutate(Collection_date = as.POSIXct(date, "%Y-%m-%d")) %>%
    mutate(Start_date = (as.POSIXct(FirstCase, "%Y-%m-%d"))) %>%
    # Calculate difference between first BA.1 signal and reporting date
    mutate(Diff_1st_Case = as.numeric(difftime(Collection_date, Start_date, units= "days"))) %>%
    filter(Diff_1st_Case >= -50) %>%
    # Select needed data
    select(Diff_1st_Case, new_cases_smoothed_per_million)
  # Calculate model
  glmD = glm(Delta_recode ~ Diff_1st_Case, data = data_loaded, 
             family = binomial(link = "logit"))
  # Create dummy dataset
  newdata = data.frame(Diff_1st_Case = 
                         as.numeric(c(min(data_loaded$Diff_1st_Case, na.rm = TRUE):
                                        200)))
  # Predict Delta fraction over time
  predicted <- as.data.frame(predict(glmD, newdata, 
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]      
  predicted <- predicted %>% 
    rename(Diff_1st_Case = 1) %>%
    # Make day difference numeric
    mutate(Diff_1st_Case = as.numeric(Diff_1st_Case)) %>%
    # Correct day difference as predicted begins with 1 instead of e.g. -50
    mutate(Diff_1st_Case= Diff_1st_Case + min(data_loaded$Diff_1st_Case) -1
    )
  # Merge data on cases and Omicron fraction
  MERGED <- merge(predicted, data_cases, by = "Diff_1st_Case") %>%
    mutate(Delta_Cases = round((Value*new_cases_smoothed_per_million
    ), 
    digits = 3))
  
  # Return dataframe
  return(MERGED)
}



##### APPLY FUNCTION TO ESTIMATE DELTA CASES ----------------------------------
# Algeria
DELTA_ALG <- MODEL_DELTA(dataset = data,
                         CountryX = "Algeria",
                         FirstCase = "yyyy-mm-dd")
# Angola
DELTA_ANG <- MODEL_DELTA(dataset = data,
                         CountryX = "Angola",
                         FirstCase = "yyyy-mm-dd")

# ... (repreat for every included country)


#### Merge predicted values for Delta -----------------------------------------
## Join predicted country data in one list to join them in a data frame
DELTA_list <- list(DELTA_ALG,
                   DELTA_ANG,
                   DELTA_BEN,
                   DELTA_BUF,
                   DELTA_CAM,
                   DELTA_ETP,
                   DELTA_GAM,
                   DELTA_GHA,
                   DELTA_KEN,
                   DELTA_MAL,
                   DELTA_MAR,
                   DELTA_NMB,
                   DELTA_NIG,
                   DELTA_SEN,
                   DELTA_SAF,
                   DELTA_TOG,
                   DELTA_UGA,
                   DELTA_ZIM,
                   DELTA_MAD,
                   DELTA_COG,
                   DELTA_MOZ)

# Merge country-specific Delta cases into one dataframe
DELTA_merged <- DELTA_list %>% reduce(full_join, by='Diff_1st_Case') %>%
  select(-contains("smoothed"),
         -contains("Value")) %>%
  rename(Time = Diff_1st_Case)


## Calculate median of predicted Delta cases
DELTA_merged$Median <- round(rowMedians(as.matrix(DELTA_merged[,2:ncol(DELTA_merged)]), na.rm = TRUE), digits = 3)
## Calculate mean of predicted Delta cases
DELTA_merged$Mean <- round(rowMeans(DELTA_merged[,2:(ncol(DELTA_merged)-1)], na.rm = TRUE), digits = 3)




#### ESTIMATE R0 FOR BA.1 ---------------------------------------------------
## Define the sliding window size for estimation of R0
t_start <- seq(2, nrow(Predicted_merged)-9)   
t_end <- t_start + 9

## Estimate R0
res_parametric_si <- estimate_R(Predicted_merged$Mean, 
                                # Select prediction method
                                method="parametric_si",
                                config = make_config(list(
                                  # define serial interval length
                                  mean_si = 3.3, 
                                  # define serial interval deviation
                                  std_si = 2.4, 
                                  t_start = t_start, 
                                  t_end = t_end
                                )))

## Extract needed results from R0 estimation results
R_results <- as.data.frame(res_parametric_si$R) %>%
  rename(Mean = 3, Q005 = 6, Q0025= 5, Q025 = 7, Median =8, Q075 = 9, 
         Q095 = 10, Q0975 = 11) %>%
  # Correct time as R0 prediction starts with day 0 instead of e.g., -50
  mutate(Time = t_end - 50)

## Extract Omicron incidence data
I_results <- as.data.frame(res_parametric_si$I) %>%
  rename(Incidence = 1)
I_results$row_names <- row.names(I_results) 
I_results <- as.data.frame(I_results) %>%
  rename(Time = 2) %>%
  mutate(Time = as.numeric(Time)-50)

## Create test plot of R data
plot(res_parametric_si, legend = FALSE)


## Merge BA.1 and Delta incidence data
I_Delta_joined <- inner_join(I_results, DELTA_merged %>% select(1,ncol(.)),
                             by = "Time") %>%
  rename(Omicron = Incidence, Delta = Mean) %>%
  pivot_longer(cols = -2, names_to = "Variant")
I_Delta_joined$Variant <- 
  factor(I_Delta_joined$Variant, 
         levels = c("Delta", "Omicron")) 



#### PLOT R0 ESTIMATIONS AND CASES --------------------------------------------
ggplot(data = R_results, aes(x = Time, y = Median))+
  # Plot daily Delta and BA.1 cases
  geom_bar(data = I_Delta_joined, 
           aes(x = Time, y = value/5, fill = Variant), stat = "identity",
           alpha = 0.3)+
  scale_fill_manual(values=c("#512DAD", "#3F9BBA"))+
  # add trend line for mean daily SARS-CoV-2 cases
  geom_smooth(data = Total_merged, aes(x = Time, y = Mean/5), 
              color = "orange",
              method = loess, fullrange = TRUE, se = FALSE)+
  # add horizontal line to show R0/RT = 1
  geom_hline(yintercept=1, linetype="dotted", color = "black", size=0.5) +
  # add RT confidence intervals
  geom_ribbon(aes(ymin=Q005, ymax=Q095), linetype=0, alpha=0.2)+
  geom_ribbon(aes(ymin=Q025, ymax=Q075), linetype=0, alpha=0.5)+
  # add point at RT maximum
  geom_point(aes(x = R_results$Time[R_results$Median == max(R_results$Median)], 
                 y = max(Median)), size = 5, colour = "#F5345E")+
  geom_line(size = 1)+
  # set axis scales
  scale_x_continuous(limits = c(-2, 50.5), expand = c(0,0))+
  scale_y_continuous(name="Median R0", sec.axis = sec_axis(~ 5*., name="BA.1 incidence/1 million"),
                     limits= c (0,12), expand =c(0,0)) +
  # add theme
  theme_classic() +
  ggtitle("R0 of Omicron in Africa") +
  xlab("Days after first Omicron detection")+
  labs(tag = paste("maximum R0:", round(max(R_results$Median), digits = 1)), 
       sep = " ") +
  theme(plot.tag.position = c(0.1, 0.6),
        plot.tag = element_text(size = 12, hjust = 0, colour = "#B81135"),
        panel.spacing.y = unit(1.5, "lines")) 



## Save plot
ggsave("R0_16-09-2022.pdf", width = 7, height = 4)