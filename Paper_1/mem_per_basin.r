library(ggplot2)
library(openxlsx)
library(jsonlite)
library(ncdf4)
library(raster)
library(gdistance)
library(dplyr)
library(terra)
library(adespatial)
library(sf)


HOME_ <- "."
italy <- st_read("./Italy_shp/Italy.shp")
rast <- raster(paste(path.expand("~"), "PHD/ISPRA_20152017_Analysis/Clustering/Unknown_effect/chl.nc", sep = "/"))
italy <- st_transform(italy, crs = crs(rast))
is_water <- rasterize(italy, rast, field = 0, background = 1)

fun_perm <- function(x){

  if(x[1] + x[2] == 2){
    return(1) #Connected 
    }

  else {
    return(0) #Not connected
  }
}

tr <- transition(is_water, fun_perm, directions = 16)
trc <- geoCorrection(tr, type = "c", scl = "FALSE")

stat_info <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% 
dplyr::select(c(id, Longitude, Latitude, Basin, Region)) %>% mutate(
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


morans_per_basin <- 
sapply(
    cosDist_per_basin, function(x) {
       morans <- adespatial::dbmem(x, MEM.autocor = "all")
    }
)

pos_eigens <- sapply(
    morans_per_basin, function(x) {
        sum(attr(x, "values") > 0)
    }
)

morans_per_basin[[1]][, seq(1, pos_eigens[[1]])]
selected_eigens <- lapply(seq_along(morans_per_basin), function(i) {
    morans_per_basin[[i]][, seq(1, pos_eigens[[i]])]
})

names(selected_eigens) <- names(morans_per_basin)

selected_eigens
output_list <- sapply(
    stat_info %>% pull(New_basin) %>% unique(),
    function(name) {
    data.frame(id = stat_info %>% dplyr::filter(New_basin == name) %>% pull(id), selected_eigens[[name]])
}
)

wb <- createWorkbook()
for (i in seq_along(output_list)) {
    addWorksheet(wb, names(output_list)[i])
    writeData(wb, names(output_list)[i], output_list[[i]])
}

saveWorkbook(wb, "./MEMs_per_basin.xlsx", overwrite = TRUE)
