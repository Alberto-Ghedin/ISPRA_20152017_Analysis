library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(ggspatial)

IMG_FORMAT <- "svg"
italy <- st_read("./Italy_shp/Italy.shp")
surroundings <- st_read("./Surrounding_shp/Surrounding.shp")
basins <- st_read("./Basins_shp/Basins.shp")
from_region_to_abbreviation <- list(
    "Friuli-Venezia-Giulia" = "FVG",
    "Veneto" = "VEN", 
    "Emilia-Romagna" = "EMR",
    "Marche" = "MAR",
    "Abruzzo" = "ABR",
    "Molise" = "MOL",
    "Puglia" = "PUG",
    "Basilicata" = "BAS",
    "Calabria" = "CAL",
    "Campania" = "CAM", 
    "Lazio" = "LAZ",
    "Toscana" = "TOS",
    "Liguria" = "LIG",
    "Sicilia" = "SIC",
    "Sardegna" = "SAR"
)
phyto_abund <- read.csv("./phyto_abund.csv") %>% 
mutate(
    Basin = case_when(
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
)

basins_centers <- basins

# Compute centroids after transforming to EPSG:3857 and back to original CRS
basins_centers$geometry <- st_transform(basins$geometry, 3857) %>%
  st_centroid() %>%
  st_transform(st_crs(basins))

region_centers <- italy %>%
  filter(region %in% names(from_region_to_abbreviation))

region_centers <- region_centers %>%
  mutate(geometry = st_transform(geometry, 3857) %>%
                      st_centroid() %>%
                      st_transform(st_crs(italy)))

# Map region names to abbreviations
region_centers <- region_centers %>%
  mutate(region = from_region_to_abbreviation[region])

abbrv_positions <- list(
  "ABR" = c(13.7, 42.23197), 
  "PUG" = c(16.48, 40.9), 
  "BAS" = c(15.9, 40.48), 
  "CAL" = c(16.34784, 39.07882),
  "CAM" = c(14.6, 40.8658),
  "EMR" = c(11.03106, 44.53113),
  "FVG" = c(12.7, 46.16188),
  "LAZ" = c(12.4, 41.98468),
  "LIG" = c(9.15, 43.65),
  "MAR" = c(12.9, 43.35484),
  "MOL" = c(14.2, 41.6),
  "SAR" = c(8.8, 40.09751),
  "SIC" = c(14.14632, 37.59335),
  "TOS" = c(11.12342, 43.45912),
  "VEN" = c(11.75, 45.6)
)

# Replace geometry column with the new points from abbrv_positions
region_centers <- region_centers %>%
  mutate(geometry = st_sfc(lapply(abbrv_positions, function(pos) st_point(pos)), crs = st_crs(italy)))

p <- ggplot() +
  geom_sf(data = surroundings, fill = "grey", color = "black") +
  geom_sf(data = italy, fill = "lightgrey", color = "black") + 
   geom_text(data = region_centers, aes(x = st_coordinates(geometry)[,1], 
                                             y = st_coordinates(geometry)[,2], 
                                             label = region),
                  size = 5, fontface = "bold") + 
   scale_x_continuous(breaks = seq(8, 18.5, 1.5), limits = c(8, 18.5)) +
   scale_y_continuous(breaks = seq(36, 47, 1.5), limits = c(36, 47)) +
   ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
     height = unit(3, "cm"), width = unit(3, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"
    )) + 
   annotation_scale(location = "bl", width_hint = 0.2) +
   geom_point(data = phyto_abund %>% distinct(Longitude, Latitude, Basin), 
             aes(x = Longitude, y = Latitude, fill = Basin),
             color = "black", 
             size = 2, shape = 21, stroke = 0.5) +
 labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    legend.position = "right", 
  legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )


ggsave(
    filename = "station_basin_map",
    device = IMG_FORMAT,
    plot = p,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
)
