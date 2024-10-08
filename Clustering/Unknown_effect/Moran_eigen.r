library(ggplot2)
library(openxlsx)
library(jsonlite)
library(ncdf4)
library(raster)
library(gdistance)
library(dplyr)
library(terra)
library(adespatial)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

raster
# Load the raster
rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/chl.nc", sep = "/"))
rast <- raster(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/rast_test.grd", sep = "/"))

is_water <- rast != 0
is_water[is.na(is_water)] <- 1

plot(is_water)


fun_perm <- function(x){

  if(x[1] + x[2] == 2){return(1)}

  else {return(0)}

}

tr <- transition(is_water, fun_perm, directions = 16)
trc <- geoCorrection(tr, type = "c", scl = "FALSE")

df_sites <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Stations_info.csv", sep = "/")) 
sites_coord <- df_sites %>% select(Longitude, Latitude) %>% as.matrix()


plot(read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Stations_info.csv", sep = "/")) %>% select(Longitude, Latitude))
cosDist <- gdistance::costDistance(trc, sites_coord) / 1000

dist <- as.matrix(cosDist)

cosDist
morans <- adespatial::dbmem(cosDist, MEM.autocor = "all")

morans
typeof(cosDist)

names(morans)
attributes(morans)
morans$values
barplot(attr(morans, "values"))

library(sp)
library(ade4)

s.value(sites_coord, morans[, c(1,2,3,4,5,6)])

dim(sites_coord)

plot(morans[, c(1,2,3,4,5)], SpORcoords = sites_coord, symbol = "circle")

sites_coord

length(attr(morans, "values"))
df_mor <- data.frame(id = c("eigen", df_sites$id), rbind(attr(morans, "values"), as.matrix(morans)))

write.csv(df_mor, paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/morans.csv", sep = "/"), row.names = FALSE)
