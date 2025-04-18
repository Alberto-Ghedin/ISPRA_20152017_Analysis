library(ggplot2)
library(openxlsx)
library(jsonlite)
library(ncdf4)
library(raster)
library(gdistance)
library(dplyr)
library(terra)
library(adespatial)

HOME_ <- "."
rast <- raster(paste(path.expand("~"), "PHD/ISPRA_20152017_Analysis/Clustering/Unknown_effect/rast_test.grd", sep = "/"))

is_water <- rast != 0
is_water[is.na(is_water)] <- 1

plot(is_water)

fun_perm <- function(x){

  if(x[1] + x[2] == 2){return(1)}

  else {return(0)}

}

tr <- transition(is_water, fun_perm, directions = 16)
trc <- geoCorrection(tr, type = "c", scl = "FALSE")

stat_info <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% 
dplyr::select(c(Longitude, Latitude, Basin, Region)) %>% mutate(
    New_basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Region == "PUG" & Basin == "SouthAdr" ~ "SA",
        Region == "PUG" & Basin == "Ion" ~ "SM",
        Region == "BAS" ~ "SM",
        Region == "CAL" & Basin == "Ion" ~ "SM",
        Region == "SIC" ~ "SIC", 
        Region == "CAL" & Basin == "SouthTyr" ~ "ST",
        Region == "CAM" ~ "ST",
        Region == "LAZ" & Basin == "SouthTyr" ~ "ST",
        Region == "LAZ" & Basin == "NorthTyr" ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
    ) %>% distinct()

cosDist_per_basin <- sapply(
    stat_info$New_basin %>% unique(), function(x) {
        sites_coord <- stat_info %>% dplyr::filter(New_basin == x) %>% dplyr::select(c(Longitude, Latitude)) %>% as.matrix()
        gdistance::costDistance(trc, sites_coord) / 1000
    }
)
names(cosDist_per_basin) <- stat_info$New_basin %>% unique()

dist_per_basin  <- sapply(cosDist_per_basin, as.matrix)

dist_per_basin[[1]] 
morans_per_basin <- 
sapply(
    dist_per_basin, function(x) {
       morans <- adespatial::dbmem(x, MEM.autocor = "all")
    }
)

n <- 5  # Replace with the desired size

# Generate a random matrix
random_matrix <- matrix(runif(n * n, min = 0, max = 10), nrow = n, ncol = n)

# Make the matrix symmetric
symmetric_matrix <- (random_matrix + t(random_matrix)) / 2
morans <- adespatial::dbmem(symmetric_matrix, MEM.autocor = "all")
