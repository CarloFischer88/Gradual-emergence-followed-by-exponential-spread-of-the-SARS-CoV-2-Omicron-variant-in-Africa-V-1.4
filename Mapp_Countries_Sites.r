# The goal of this script is to visualize all countries and sampling sites
# that were represented in the study


#load packages
library(tidyverse)
library(readxl)
library(xlsx)
library(raster)
library(sf)
library(ggspatial)



#### LOAD DATA ----------------------------------------------------------------
## Set working directory
setwd("~/choose/directory")
## Load dataset
data <- read.csv("dataset_PCR.csv",
                 header = TRUE,        # Whether to read the header or not
                 sep = ",",            # Separator of the values
                 dec = ".")
# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"




#### TIDY DATA ----------------------------------------------------------------
data_tidy <- data %>%
  # Change ages to completed life years
  mutate(Age = floor(as.numeric(Age))) %>%
  # Remove SARS-CoV-2 negative samples
  filter(!is.na(Variant)) %>%
  # Change date and create binned dates
  mutate(Collection_date2 = as.POSIXct(Collection_date, "%Y-%m-%d"))
# Calcualte month Year date
data_tidy$Month <- format(data_tidy$Collection_date2, format="%Y-%m") 
data_tidy <- data_tidy %>%
  # Remove tooo old samples
  filter(Month >= "2021-06")




#### SELECT DATA --------------------------------------------------------------
data_selected <- data_tidy %>%
  dplyr::select(10, 12, 13) %>%
  #remove duplicates (keep one entry per sampling site)
  distinct()

# Check for countries and entries
Sites <- as.data.frame(table(data_selected$Country))


#### LOAD SHAPEFILES ----------------------------------------------------------
setwd("~/PhD/R/Shapefiles_Geotifs")
World <- st_read("gadm36_0.shp")

## Select African and bordering countries
Africa <- World %>% 
  filter(NAME_0 %in% c("Nigeria", "Ethiopia","Egypt", 
                       "Democratic Republic of the Congo", "South Africa", 
                       "Tanzania", "Kenya", "Algeria", "Sudan", "Morocco",
                       "Uganda", "Mozambique", "Ghana", "Angola", 
                       "Côte d'Ivoire", "Madagascar", "Cameroon", "Niger",
                       "Burkina Faso", "Mali", "Malawi", "Zambia", "Somalia",
                       "Senegal", "Chad", "Zimbabwe", "Rwanda", "Tunisia",
                       "Guinea", "Benin", "Burundi", "South Sudan", "Togo", 
                       "Eritrea", "Sierra Leone", "Libya", "Republic of Congo",
                       "Central African Republic", "Liberia", "Mauritania", 
                       "Namibia", "Botswana", "Lesotho", "Gambia", "Gabon",
                       "Guinea-Bissau", "Mauritius", "Equatorial Guinea", 
                       "Swaziland", "Djibouti", "Comoros", "Western Sahara",
                       "Saudi Arabia", "Oman", "Israel", "Jordan", "Kuwait", 
                       "Lebanon", "Qatar", "United Arab Emirates", "Yemen",
                       "Iran", "Syria", "Iraq", "Turkey", "Greece", "Italy",
                       "Spain", "Portugal", "Albania", "Palestina", "Azerbaijan"))

## Dummycode which countries are represented in the study
Africa <- Africa %>% 
  mutate(Study = ifelse(NAME_0 %in% c("Algeria",
                                      "Angola",
                                      "Benin", 
                                      "Burkina Faso",
                                      "Cameroon",
                                      "Congo",
                                      "Ethiopia",
                                      "Gabon",
                                      "Gambia",
                                      "Ghana", 
                                      "Kenya",
                                      "Mali",
                                      "Morocco",
                                      "Madagascar",
                                      "Mozambique",
                                      "Namibia", 
                                      "Niger", 
                                      "Senegal",
                                      "South Africa",
                                      "Togo",
                                      "Uganda",
                                      "Zimbabwe"),
                        1, 0))





#### Plot represented countries and sampling sites ----------------------------
ggplot() +
  # plot countries and colour represented countries
  geom_sf(data = Africa, aes(fill = factor(Study)), color = "#CAD2D4", 
          size = .005, show.legend = FALSE)+
  scale_fill_manual(values=c("#F0F0EE", "#EDA26D"))+
  # add theme
  theme_classic()+
  # colour background
  theme(panel.background = element_rect(fill = "#CAD2D4",
                                        colour = "#CAD2D4",
                                        size = 0.0, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  # Plot sampling sites
  geom_point(data = data_selected, aes(x = lon, y = lat),
             fill = "#0A4766", color = "#0A4766", alpha = 0.5, size = 3,
             shape = 16) +
  # add north south arrow
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "bl", width_hint = 0.2) +
  # Choose area to be plotted
  coord_sf(xlim = c(-20, 53), ylim = c(-37, 38), expand = FALSE)+
  # Add tag
  labs(tag = "South 
Atlantic
Ocean") +
  theme(plot.tag.position = c(0.27, 0.4),
        plot.tag = element_text(size = 10, hjust = 0, colour = "#707980", face="italic"),
        panel.spacing.y = unit(1, "lines"))

## Save plot
ggsave("Map_Sites", height = 5, width = 8, dpi = 300)

