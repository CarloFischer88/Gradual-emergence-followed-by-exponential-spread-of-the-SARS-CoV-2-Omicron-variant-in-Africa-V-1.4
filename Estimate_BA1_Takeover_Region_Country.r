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


## Load population data
pop <- read_excel("Population_data.xlsx") 
## Define countries
# (Gabon is excluded as no BA.1 cases were observed in the sample)
Countries_Africa <- c("Algeria", "Angola", "Benin", "Burkina Faso", "Cameroon", "Ethiopia",
                      "Gambia", "Ghana", "Kenya", "Madagascar", "Mali", "Morocco", "Mozambique",
                      "Namibia", "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe")
# Mutate population data
pop_selection <- pop %>% rename(Country = 3) %>% 
  # Select countries
  filter(Country %in% Countries_Africa ) %>%
  # Select msot recent year (2021)
  filter(Year == max(Year)) %>%
  # Select needed ollumns
  select(Country, Year, 113) %>%
  # Show population in millions
  mutate(Total = Total*1000)




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
  
  
## Identify first BA.1 marker positive samples ------------------------------
First_Omicron <- as.data.frame(data_tidy %>% 
                                   filter(Omicron == 1) %>%
                                   group_by(Country) %>%
                                   summarise(min(Collection_date2))) %>%
    rename(Collection_date = 2) 


## Calcualte difference between first BA.1 marker positive sample and collection date
data_tidy <-   data_tidy %>%
  mutate(Day_Difference = case_when(
    Country == "Algeria" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Angola" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Benin" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Burkina Faso" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Cameroon" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Congo" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Ethiopia" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Gabon" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Gambia" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Ghana" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Kenya" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Mali" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Morocco" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Mozambique" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Namibia" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Niger" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Senegal" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "South Africa" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Togo" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Uganda" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Zimbabwe" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days"),
    Country == "Madagascar" ~ difftime(Collection_date2, "yyyy-mm-dd", units= "days")
  ),
  Day_Difference = round(Day_Difference, digits = 1),
  Days_post_June = as.numeric(round(difftime(Collection_date2, "2021-06-01", units= "days"), digits = 0)))






#### IDENTIFY DATE FOR OMICRON DOMINANCE --------------------------------------
### Define function to identify date of dominance
MODEL_DOMINANCE <- function(dataset, CountryX) {
  # Fit glm
  glmX = glm(Omicron ~ Days_post_June, data = dataset %>% filter(Country == CountryX), 
             family = binomial(link = "logit"))
  # Create empty dataset
  newdata = data.frame(Days_post_June = 0:320) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
  predicted <- as.data.frame(predict(glmX, newdata, #se.fit=TRUE,
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]   
  predicted <- predicted %>% 
    rename(Day_Difference = 1) %>%
    filter(Value > 0.5) %>%
    filter(Value == min(Value)) %>%
    mutate(Country = CountryX)
  # Return dataframe
  return(predicted)
}

dom.Alg <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Algeria")  
dom.Ang <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Angola") 
dom.Ben <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Benin")  
dom.BuF <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Burkina Faso") 
dom.Cam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Cameroon") 
dom.Eth <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ethiopia") 
dom.Gab <- MODEL_DOMINANCE(dataset = data_grouped, CountryX = "Gabon")  
dom.Gam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gambia") 
dom.Gha <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ghana") 
dom.Ken <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Kenya") 
dom.Mal <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mali") 
dom.Mor <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Morocco") 
dom.Moz <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mozambique") 
dom.Nam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Namibia") 
dom.Nig <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Niger") 
dom.Cog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Congo") 
dom.Sen <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Senegal") 
dom.SAf <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "South Africa") 
dom.Tog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Togo") 
dom.Uga <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Uganda") 
dom.Mad <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Madagascar") 
dom.Zim <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Zimbabwe") 


## Summarize takeover modelling results
Dominance_Modelling <- rbind(dom.Alg, dom.Ang, dom.Ben, dom.BuF, dom.Cam, dom.Eth,
                             dom.Gam, dom.Gha, dom.Ken, dom.Mal, dom.Mor, dom.Moz,
                             dom.Nam, dom.Nig, dom.Cog, dom.Sen, dom.SAf, dom.Tog,
                             dom.Uga, dom.Mad, dom.Zim) %>%
  # Reorder and select data
  select(3, 1)


# Merge modelling data and population
Dominance_Modelling_pop <- merge(Dominance_Modelling, pop_selection, by = "Country") %>%
  # Make days numeric
  mutate(Day_Difference = as.numeric(Day_Difference)) %>%
  # Calculate first expected cases in days post June
  mutate(First_case = round(Day_Difference - (log10(0.5/(1/Total))/log10(2))*3), digits = 0) %>%
  # Calclulate date for first expected BA.1 cases
  mutate(First_Case_date = as.POSIXct("yyyy-mm-dd", format = "%Y-%m-%d") + lubridate::days(First_case),
         # Add 30 days to calculated dates of first expected BA.1 cases
         First_Case_date_30 = First_Case_date - lubridate::days(30))





#### MERGE FIRST EXPECTED CASES AND PCR RESULTS DATA --------------------------
# Select needed collumns on first expected cases
Dominance_Modelling_pop_sel <- Dominance_Modelling_pop %>% select(1, 7)
# Merge datasets
data_cleaned <- merge(data_tidy, Dominance_Modelling_pop_sel, by = "Country")
# Remove BA.1 marker positive samples collected before first expected case
data_cleaned2 <- data_cleaned %>%
  filter(Omicron == 1 & Collection_date2 >= First_Case_date | Omicron == 0)




#### GROUP PCR DATA BY AFRICAN REGION -----------------------------------------
## Define Regions
WA = c("Benin", "Burkina Faso", "Gambia", "Ghana", "Mali", "Niger", "Senegal", "Togo")
# Gabon was ignored as the short sampling timeframe did not allow to observe
# BA.1 over time
NA <- c("Algeria", "Morocco")
CA <- c("Angola", "Cameroon", "Congo")
EA <- c("Ethiopia", "Kenya", "Mozambique", "Uganda") 
# Madagascar ignord as it is not main land
# Zimbabwe was ignored as no BA.1 cases were observed in the sample
SA <- c("Namibia", "South Africa")

## Group data
data_grouped <- data_cleaned2 %>%
  mutate(Region = ifelse(Country %in% WA, "Western Africa", 
                  ifelse(Country %in% NAf, "Northern Africa",
                  ifelse(Country %in% CA, "Central Africa",
                  ifelse(Country %in% EA, "Eastern Africa",
                  ifelse(Country %in% SA, "Southern Africa", NA)))))) %>%
  filter(!is.na(Region)) %>%
  # Calculate difference between collection date and June
  mutate(DayDiff = difftime(Collection_date2, "2021-06-01", units= "days")) %>%
  mutate(DayDiff = round(DayDiff, digits = 0),
         DayDiff = as.numeric(DayDiff))




#### MODELL TAKEOVER ----------------------------------------------------------
ggplot() +
  geom_smooth(data = data_grouped,
              aes(x = Collection_date2, y = Omicron, colour = Region, fill = Region),
              method.args=list(family="binomial"),
              method = glm,
              fullrange = TRUE, se = TRUE) +
  # add theme
  theme_classic() +
  # define x axis limits
  expand_limits(x = as.POSIXct(c("2021-06-01", "2022-08-01")))+
  # add line for report of Omicron emerence
  geom_vline(xintercept = as.POSIXct("2021-11-27"))

## Save plot
ggsave("Takeover_Curves.pdf", height = 4, width = 7)





#### IDENTIFY DATE FOR OMICRON DOMINANCE PER AFRICAN REGION -------------------
## Define function to identify date of dominance
MODEL_DOMINANCE_REGION <- function(dataset, RegionY) {
  # Fit glm
  glmY = glm(Omicron ~ Days_post_June, data = dataset %>% filter(Region == RegionY), 
             family = binomial(link = "logit"))
  # Create empty dataset
  newdata = data.frame(Days_post_June = 0:320) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
  predicted <- add_ci(newdata, 
                      glmY, 
                      names = c("lwr", "upr"), alpha = 0.1) 
  predicted <- setDT(predicted, keep.rownames = TRUE)[]   
  pred_pred <- predicted %>% filter(pred >= 0.5) %>% filter(pred == min(pred)) %>%
    mutate(Date_pred = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    select(6)
  pred_lwr <- predicted %>% filter(lwr >= 0.5) %>% filter(lwr == min(lwr))%>%
    mutate(Date_lwr = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    select(6)
  pred_upr <- predicted %>% filter(upr >= 0.5) %>% filter(upr == min(upr))%>%
    mutate(Date_upr = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    select(6)
  predicted_dates <- bind_cols(Region = RegionY, pred_pred, pred_lwr, pred_upr)
  # Return dataframe
  return(predicted_dates)
}

       

## Apply function to identify dominance per Region
Dom_WA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                        Region = "Western Africa")
Dom_CA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                        Region = "Central Africa")
Dom_EA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                        Region = "Eastern Africa")
Dom_SA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                        Region = "Southern Africa")
Dom_NA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                        Region = "Northern Africa")
       
## Merge Modelling results
Dom_Merged <- bind_rows(Dom_CA, Dom_EA, Dom_NA, Dom_SA, Dom_WA) %>%
         # Calculate difference in days
         mutate(lwr_days = difftime(Date_lwr, Date_pred, units= "days"),
                upr_days = difftime(Date_upr, Date_pred, units= "days"))
       # Export table 
       write.csv(Dom_Merged,
                 "~/directory\\Dom_Merged.csv")
       
 
       
       
#### ANALYSE THE TIME BETWEEN FIRST CASES AND DOMINANCE -----------------------
## Calculate first cases within acceptable range
First_Cases_filtered <- merge(Dominance_Modelling, data_tidy, by = "Country") %>%
         rename(Days_Dominant = 2) %>%
         filter(Omicron == 1 & Days_post_June <= Days_Dominant) %>% 
         group_by(Country) %>%
         summarise(min(Days_post_June)) %>%
         rename(First_Case = 2) 
       
## Calculate days until dominance and remove samples with to short timeframe
Data_T_Dominance <- merge(First_Cases_filtered, Dominance_Modelling, by = "Country") %>%
         mutate(Days_until_Dominance = as.numeric(Day_Difference) - First_Case) %>%
         # remove samples with too short time between dominance and first case
         mutate(Days_10 = round(Days_until_Dominance - (log10(0.5/(0.1))/log10(2))*3, digits = 0)) %>%
         filter(Days_10 >= 0)
       
## Calculate median days until dominance
       Median_Days_Dom <- as.data.frame.matrix(groupwiseMedian(
         Days_until_Dominance ~ 1,
         data = Data_T_Dominance,
         conf = 0.95,
         bca = FALSE,
         R = 1000,
         exact = TRUE,
         digits = 2)) %>%
         rename(Lower = 5, Upper = 6)
       
       
## Plot median time until dominance
ggplot(data = Data_T_Dominance, aes(x = 1, y = as.numeric(Days_until_Dominance)))+
         geom_jitter(width = 0.1, height = 0, alpha = 0.5, size = 5, colour = "orange")+
         geom_point(size = 2, aes(x = 1, y = median(Days_until_Dominance))) +
         geom_errorbar(aes(ymin = Median_Days_Dom$Lower, ymax = Median_Days_Dom$Upper))+
         scale_x_continuous(limits = c(0.5, 1.5))+
         expand_limits(y=c(0, 100)) +
         theme_classic()

## Save plot
ggsave("Median_Difference_Dominance.pdf", height = 4, width = 7)
       
       
       