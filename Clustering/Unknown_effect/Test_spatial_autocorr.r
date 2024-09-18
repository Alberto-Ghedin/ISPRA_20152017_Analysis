library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(tibble)
library(jsonlite)
library(tidyverse)
library(vegan)
library(stats)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")


# Load the raster
rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/chl.nc", sep = "/"))
rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/rast_test.grd", sep = "/"))

is_water <- rast != 0
is_water[is.na(is_water)] <- 1

env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) %>% dplyr::select(-c(Secchi_depth)) %>% dplyr::filter(!Region %in% c("Lazio", "Basilicata", "Calabria")) %>% na.omit()

covariates <- c(
  "T",
  "Salinity",
  "O_sat",
  "pH",
  "Chla",
  "NO3",
  "NO2",
  "NH4",
  "TN",
  "PO4",
  "TP",
  "SiO4"
)


env_data_std <- env_data %>% 
    dplyr::select(all_of(covariates)) %>% 
    dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dist(method = "euclidean")


clust <- hclust(env_data_std, method = "ward.D2")



env_data$cluster <- cutree(clust, k = 15)
env_data %>% dplyr::select(Region, cluster) %>% table() %>% as.data.frame() %>% 
ggplot() + geom_tile(aes(x = Region, y = cluster, fill = Freq)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Clusters per region", x = "Region", y = "Cluster")
