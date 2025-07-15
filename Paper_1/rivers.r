library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(tibble)
library(openxlsx)
library(ggspatial)
library(sf) 

HOME_ <- "./Paper_1"
rivers <- openxlsx::read.xlsx(paste(HOME_, "EFAS_rivers.xlsx", sep = "/"))
phyto_abund <- read.csv(paste(HOME_,"phyto_abund.csv", sep ="/"))
stations <- phyto_abund %>% dplyr::distinct(id, Longitude, Latitude)
rivers_mouths <- rivers %>% dplyr::select(rivername, lon_mouth, lat_mouth)

italy <- st_read(paste(HOME_, "Italy_shp/Italy.shp", sep = "/"))
surroundings <- st_read(paste(HOME_,"Surrounding_shp/Surrounding.shp", sep = "/"))
basins <- st_read(paste(HOME_,"Basins_shp/Basins.shp", sep = "/"))

ordered_transect <- c(
    "SMTS", "SMLG", "VENEZIA", "ROSOLINA", "PORTO_GARIBALDI", "CESENATICO", "RIMINI", "Chienti", "Esino", "GU",
    "VA", "R14001_B2", "FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO", "PORTO_CESAREO", "PUNTA_RONDINELLA", "SINNI", "Villapiana",
    "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche", "Isole_Ciclopi", "Plemmirio", "Isola_Correnti", "San_Marco", "Isole_Egadi", "Capo_Gallo","Vibo_marina", 
    "Cetraro", "Cilento", "Salerno", "Napoli", "Domizio", "m1lt01", "m1lt02",  "m1rm03",  "m1vt04", "Collelungo","Carbonifera",
    "Donoratico", "Fiume_Morto", "Mesco" , "Portofino"  ,  "Voltri"  , "Quiliano" , 
    "Olbia", "Arbatax", "Villasimius", "Cagliari", "Oristano", "Alghero", "Porto_Torres"

)
sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))

stations <- merge(
    stations, 
    sea_depth %>% dplyr::select(id, SeaDepth, Transect)
)
stations$Transect <- factor(stations$Transect, levels = ordered_transect, ordered = TRUE)
stations <- stations %>% group_by(Transect) %>% 
    summarise(
        Longitude = mean(Longitude, na.rm = TRUE),
        Latitude = mean(Latitude, na.rm = TRUE), 
        .groups = "drop"
    )

library(ggrepel)
ggplot() + 
geom_sf(data = italy, fill = "lightgrey", color = "black") + 
geom_point(data = stations, aes(x = Longitude, y = Latitude), color = "blue", size = 2) + 
geom_text_repel(
  data = stations %>% group_by(Transect) %>% 
    summarise(
      Longitude = mean(Longitude, na.rm = TRUE),
      Latitude = mean(Latitude, na.rm = TRUE), 
      .groups = "drop"
    ),
  aes(x = Longitude, y = Latitude, label = Transect),
  size = 3
) +
geom_point(data = rivers_mouths, aes(x = lon_mouth, y = lat_mouth), color = "red", size = 1) +
geom_text_repel(
  data = rivers_mouths,
  aes(x = lon_mouth, y = lat_mouth, label = rivername),
  size = 3
) +
labs(title = "Stations and Rivers Mouths",
     x = "Longitude",
     y = "Latitude") +
theme_minimal()

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


reg_stations <- phyto_abund %>% dplyr::filter(Region == "PUG") %>% dplyr::distinct(id, Longitude, Latitude)
reg_stations
# Convert to sf objects for distance calculation
reg_sf <- st_as_sf(reg_stations, coords = c("Longitude", "Latitude"), crs = 4326)
rivers_sf <- st_as_sf(rivers_mouths, coords = c("lon_mouth", "lat_mouth"), crs = 4326)

# Find nearest river for each CAM station
nearest_river_idx <- st_nearest_feature(reg_sf, rivers_sf)
reg_sf$nearest_river <- rivers_mouths$rivername[nearest_river_idx]
reg_sf$riv_dist <- as.numeric(st_distance(reg_sf, rivers_sf[nearest_river_idx, ], by_element = TRUE) / 1000)

riv_dis <- reg_sf %>% as.data.frame() %>% 
merge(
  rivers %>% dplyr::select(rivername, `MEAN_2011_2023.m3s-1`), 
  by.x = "nearest_river",
  by.y = "rivername"
) %>% 
mutate(
  riv_discharge = `MEAN_2011_2023.m3s-1` / riv_dist ^ 2
) 

sample_abund %>% dplyr::filter(Transect %in% c("Cilento", "Salerno")) %>% pull(id) %>% unique()

phyto_abund %>% dplyr::filter(Region == "CAM") %>% 
dplyr::group_by(Date, id) %>% 
dplyr::summarise(
  sample_abund = sum(Num_cell_l), 
  .groups = "drop"
) %>% 
merge(
  stations, 
  by = "id",
) %>%
merge( 
cam_sf %>% dplyr::select(id, riv_dist) %>% as.data.frame(), 
  by = "id",
) %>% 
lm(
  log10(sample_abund) ~ riv_dist, 
  data = .
) %>% summary()

# Compute distance from each station in reg_stations to Tevere river mouth
tevere_mouth <- rivers_mouths %>% filter(rivername == "Tevere")
tevere_sf <- st_as_sf(tevere_mouth, coords = c("lon_mouth", "lat_mouth"), crs = 4326)

reg_stations$tevere_dist_km <- as.numeric(st_distance(reg_sf, tevere_sf[rep(1, nrow(reg_sf)), ], by_element = TRUE) / 1000)
reg_stations$tevere_dist_km <- as.numeric(st_distance(reg_sf, tevere_sf, by_element = TRUE) / 1000)
reg_stations


library(stringr)
chem_phys %>% filter(str_detect(id, "IT_m1rm03"))

chem_phys %>% merge(
  sample_abund %>% dplyr::select(id, Transect, Region),
  by = "id"
) %>% dplyr::select(Region == "LAZ") %>% pull(id) %>% unique()

ids <- phyto_abund %>% dplyr::filter(Region == "LAZ") %>% pull(id) %>% unique()

ids
