library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(tibble)
library(jsonlite)
library(tidyverse)
library(vegan)
library(stats)
library(gdistance)
library(raster)
library(spdep)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")


# Load the raster
rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/chl_phyc.nc", sep = "/"))
#rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/rast_test.grd", sep = "/"))

df_sites <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Stations_info.csv", sep = "/")) 
sites_coord <- df_sites %>% dplyr::select(Longitude, Latitude) %>% as.matrix()


is_water <- rast
is_water[rast != 0] <- 1
is_water[is.na(is_water)] <- 0



fun_perm <- function(x){

  if(x[1] + x[2] == 2){return(1)}

  else {return(0)}

}

tr <- transition(is_water, fun_perm, directions = 16)
trc <- geoCorrection(tr, type = "c", scl = "FALSE")
cosDist <- gdistance::costDistance(trc, sites_coord) / 1000
dist <- as.matrix(cosDist)

cells <- cellFromXY(is_water, df_sites %>% filter(dist[,1] == Inf)  %>% dplyr::select(all_of(c("Longitude", "Latitude"))))
is_water[cells] <- 1

tr <- transition(is_water, fun_perm, directions = 16)
trc <- geoCorrection(tr, type = "c", scl = "FALSE")
cosDist <- gdistance::costDistance(trc, sites_coord) / 1000
dist <- as.matrix(cosDist)

wij <- 1 / dist
diag(wij) <- 0



env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) %>% dplyr::select(-c(Secchi_depth)) %>% dplyr::filter(!Region %in% c("Lazio", "Basilicata", "Calabria"))
params <- fromJSON(paste(HOME_ , "/ISPRA_20152017_Analysis/params.json", sep = "/"))
seasons <- params[["seasons"]]
seasons <- rep(names(seasons), each = 3)

#create season column 
env_data$Date <- as.Date(env_data$Date, format = "%Y-%m-%d")
env_data <- env_data %>% mutate(Season = factor(seasons[as.integer(format(Date, "%m"))]), Year = factor(format(Date, "%Y"))) %>% mutate(Season_year = paste(Season, Year, sep = "_"))

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


season_means <- env_data %>% 
  dplyr::select(-c(Date, Region, Year, Season_year)) %>% 
  group_by(Season, id) %>% 
  summarise_all(mean, na.rm = TRUE) %>%
  ungroup() %>%
  split(.$Season)

season_col <- rep(names(season_means), each = length(covariates))
covariate_col <- rep(covariates, length(names(season_means))) 
statistic_col <- rep(0, length(covariates) * length(names(season_means)))
p_vale_col <- rep(0, length(covariates) * length(names(season_means)))
i <- 1
for (season in names(season_means)){
  print(paste("Season: ", season))
  for (covariate in covariates) {
    print(paste("Covariate: ", covariate))
    values <-  season_means[[season]][[covariate]] %>% na.omit() %>% as.numeric()
    not_na <- !is.na(season_means[[season]][[covariate]])
    ids <- season_means[[season]]$id[not_na]
    wij_subset <- wij[ids, ids]
    lw <- spdep::mat2listw(wij_subset, row.names = ids, style = "W", zero.policy = TRUE)
    moran_test_result <- spdep::moran.test(values, lw, randomisation = FALSE, alternative = "two.sided")
    statistic_col[i] <- moran_test_result$statistic
    p_vale_col[i] <- moran_test_result$p.value
    i <- i + 1
    }
}

moran_statistics <- data.frame(Season = season_col, Covariate = covariate_col, Statistic = statistic_col, P_value = p_vale_col) 

moran_statistics %>% write.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/seasonal_spatial_autocorr.xlsx", sep = "/"), row.names = FALSE)