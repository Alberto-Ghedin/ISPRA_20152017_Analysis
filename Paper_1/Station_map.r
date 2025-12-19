library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(ggspatial)
library(rjson)

IMG_FORMAT <- "pdf"
HOME_ <- "."
source(file.path(HOME_, "utils.r"))
italy <- st_read(paste(HOME_, "Italy_shp/Italy.shp", sep = "/"))
surroundings <- st_read(paste(HOME_, "Surrounding_shp/Surrounding.shp", sep = "/"))
basins <- st_read(paste(HOME_, "Basins_shp/Basins.shp", sep = "/"))
sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% dplyr::select(id, Transect,SeaDepth)
)
phyto_abund$Region <- from_region_to_abreviation[as.character(phyto_abund$Region)]
phyto_abund$Transect <- factor(phyto_abund$Transect, levels = ordered_transect, ordered = TRUE)
phyto_abund <- phyto_abund %>% mutate(
    Basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Transect %in% c("FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO") ~ "SA",
        Transect %in% c("PORTO_CESAREO", "PUNTA_RONDINELLA") ~ "SM",
        Region == "BAS" ~ "SM",
        Transect %in% c("Villapiana", "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche") ~ "SM",
        Region == "SIC" ~ "SIC", 
        Transect %in% c("Vibo_marina", "Cetraro") ~ "ST",
        Region == "CAM" ~ "ST",
        Transect %in% c("m1lt01", "m1lt02") ~ "ST",
        Transect %in% c("m1rm03", "m1vt04") ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
)
phyto_abund$Basin <- factor(phyto_abund$Basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)


basins_centers <- basins

# Compute centroids after transforming to EPSG:3857 and back to original CRS
basins_centers$geometry <- st_transform(basins$geometry, 3857) %>%
  st_centroid() %>%
  st_transform(st_crs(basins))

region_centers <- italy %>%
  filter(region %in% names(from_region_to_abreviation))

region_centers <- region_centers %>%
  mutate(geometry = st_transform(geometry, 3857) %>%
                      st_centroid() %>%
                      st_transform(st_crs(italy)))

# Map region names to abbreviations
region_centers <- region_centers %>%
  mutate(region = from_region_to_abreviation[region])

abbrv_positions <- list(
  "ABR" = c(13.8, 42.23197), 
  "PUG" = c(16.6, 41), 
  "BAS" = c(16.2, 40.48), 
  "CAL" = c(16.5, 39.3),
  "CAM" = c(14.7, 40.95),
  "EMR" = c(11.03106, 44.53113),
  "FVG" = c(12.95, 46.16188),
  "LAZ" = c(12.55, 41.98468),
  "LIG" = c(9.12, 44.5),
  "MAR" = c(13.2, 43.37),
  "MOL" = c(14.5, 41.6),
  "SAR" = c(9.1, 40.09751),
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
   geom_text(data = region_centers %>% dplyr::filter(!(region %in% c("LIG", "MOL"))), aes(x = st_coordinates(geometry)[,1], 
                                             y = st_coordinates(geometry)[,2], 
                                             label = region),
                  size = 5, fontface = "bold") + 
  geom_text_repel(
    data = region_centers %>% dplyr::filter(region %in% c("LIG", "MOL")),
    aes(x = st_coordinates(geometry)[,1], 
        y = st_coordinates(geometry)[,2], 
        label = region),
    size = 5, fontface = "bold",
    nudge_x = c(-1, 1),
    nudge_y = c(-1, 1),
    segment.color = "black"
  ) + 
   scale_x_continuous(breaks = seq(8, 18.5, 1.5), limits = c(8, 18.5)) +
   scale_y_continuous(breaks = seq(36, 47, 1.5), limits = c(36, 47)) +
   ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
     height = unit(3, "cm"), width = unit(3, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20",
      #text_family = "ArcherPro Book"
    )) + 
   annotation_scale(location = "bl", width_hint = 0.2) +
  geom_point(data = phyto_abund %>% distinct(Longitude, Latitude, Basin), 
         aes(x = Longitude, y = Latitude, fill = Basin),
         color = "black", 
         size = 2, shape = 21, stroke = 0.5) +
  scale_fill_manual(
    values = setNames(
      colorBlindness::paletteMartin[c(seq(1, 9, by = 2), seq(2, 8, by = 2))],
      c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR")
    ),
    name = "Basin"
  ) +
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
p
ggsave(
    filename = paste(HOME_, paste("station_basin_map", IMG_FORMAT, sep = "."), sep = "/"),
    plot = p,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
)


p <- ggplot() +
  geom_sf(data = italy, fill = "lightgrey", color = "black") + 
  geom_text(data = region_centers %>% dplyr::filter(!(region %in% c("LIG", "MOL"))), aes(x = st_coordinates(geometry)[,1], 
                                             y = st_coordinates(geometry)[,2], 
                                             label = region),
                  size = 5, fontface = "bold") + 
  geom_point(data = phyto_abund %>% dplyr::filter(Region %in% c("CAL", "CAM"), Basin == "ST") %>% group_by(Transect) %>% 
  mutate(Longitude = mean(Longitude), Latitude = mean(Latitude)),
         aes(x = Longitude, y = Latitude),
         fill = "black",
         color = "black",
         size = 2) + 
  geom_text_repel(
    data = phyto_abund %>% dplyr::filter(Region %in% c("CAL", "CAM"), Basin == "ST") %>% group_by(Transect) %>% 
  summarise(Longitude = mean(Longitude), Latitude = mean(Latitude)) %>% 
  mutate(Transect = ifelse(Transect == "Vibo_marina", "Vibo Marina", as.character(Transect))),
         aes(x = Longitude, y = Latitude, label = Transect),
    size = 5, fontface = "bold",
    segment.color = "black"
  ) + 
  labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    legend.position = "right", 
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  ) +
  xlim(13.5, 16.5) + 
  ylim(38.5, 41) 
ggsave(
    filename = paste(HOME_, paste("station_map_ST_cal_cam", IMG_FORMAT, sep = "."), sep = "/"),
    plot = p,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
)
  

phyto_abund %>% dplyr::filter(Region %in% c("CAL", "CAM"), Basin == "ST") %>% 
  group_by(Transect) %>% 
  summarise(Longitude = mean(Longitude), Latitude = mean(Latitude)) %>% 
  ungroup() %>%
  mutate(Transect = ifelse(Transect == "Vibo_marina", "Vibo Marina", as.character(Transect)))
