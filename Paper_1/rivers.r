library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(tibble)
library(openxlsx)


rivers <- openxlsx::read.xlsx("./EFAS_rivers.xlsx")
phyto_abund <- read.csv("./phyto_abund.csv")
phyto_abund %>% head()
stations <- phyto_abund %>% dplyr::distinct(id, Longitude, Latitude)
rivers_mouths <- rivers %>% dplyr::select(rivername, lon_mouth, lat_mouth)

library(terra)
stations_vect <- terra::vect(stations, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
rivers_vect <- vect(rivers, geom = c("lon_mouth", "lat_mouth"), crs = "EPSG:4326")

stations_proj <- project(stations_vect, "EPSG:3857")
rivers_proj <- project(rivers_vect, "EPSG:3857")

buffers <- buffer(stations_proj, width = 20000) 

station_hits <- relate(buffers, rivers_proj, "intersects")

station_hits %>% dim()

result <- tibble(
  station_id = stations$id,
  nearby_rivers = lapply(station_hits, function(idxs) rivers$river_name[idxs])
)

result_long <- data.frame(station_id = c(), river_name = c())

result_long

for (i in seq_along(length(stations$id))) {
  sid <- stations$id[i]
  nearby <- rivers$rivername[station_hits[i, ]]
  
  if (length(nearby) > 0) {
    result_long <- bind_rows(result_long, data.frame(
      station_id = sid,
      river_name = nearby
    ))
  }
}

result_long

