library(openxlsx)
library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(patchwork)
library(ggplot2)

HOME_ <- "."
IMAGE_FORMAT <- "pdf"
source(file.path(HOME_, "utils.r"))

sheets <- getSheetNames(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"))
IndVal <- sapply(sheets, function(sheet) {
    data <- openxlsx::read.xlsx(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), sheet = sheet, colNames = TRUE)
    colnames(data)[1] <- "Taxon"
    data <- data %>% dplyr::filter(Taxon != "Other phytoplankton")
    return(data)
}, simplify = FALSE)
names(IndVal) <- sheets

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% select(id, Transect,SeaDepth)
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


abund_only_genera <- phyto_abund %>% 
    dplyr::filter(Det_level %in% c("Genus", "Species")) %>%
    group_by(Date, id, Basin, Genus) %>% 
    summarise(Abund = sum(Num_cell_l), Season = first(Season), Class = first(Class), .groups = "drop")


quantile_abund_genera <- abund_only_genera %>% 
    group_by(Season, Basin, Genus) %>%
    summarise(
        Abund =median(Abund),
    .groups = "drop"
)

pie_charts <- sapply(
    names(IndVal),
    function(basin) {
        selected_genera <- IndVal[[basin]] %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
        ungroup() %>% pull(Taxon)
        dominant_genera <- quantile_abund_genera %>% dplyr::filter(Genus %in% selected_genera & Basin == basin) %>% 
            group_by(Season, Basin) %>% 
            mutate(
                rel_abund = Abund / sum(Abund), 
                .groups = "drop"
            ) %>% dplyr::filter(rel_abund > 0.05) %>% pull(Genus) %>% unique()
        pie_data <- quantile_abund_genera %>%
            dplyr::filter(Genus %in% dominant_genera & Basin == basin) %>%
            mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
            group_by(Season) %>%
            mutate(rel_abund = Abund / sum(Abund))

        pie_charts <- lapply(levels(pie_data$Season), function(season) {
            season_data <- pie_data %>% filter(Season == season)
            label_pos <- season_data %>% 
            mutate(csum = rev(cumsum(rev(rel_abund))), 
            pos = rel_abund/2 + lead(csum, 1),
            pos = if_else(is.na(pos), rel_abund/2, pos))

            ggplot(season_data, aes(x = "", y = rel_abund, fill = Genus)) +
                geom_bar(stat = "identity", width = 1, color = "black") +
                coord_polar("y") +
                labs(title = season, x = NULL, y = NULL) +
                geom_label_repel(
                    data = label_pos,
                    aes(y = pos, label = Genus),
                    size = 4.5, 
                    nudge_x = 1, 
                    fill = "white", 
                    color = "black",
                    show.legend = FALSE
                    ) +

                theme_void() +
                scale_fill_manual(values = unname(colorBlindness::paletteMartin)) +
                theme(legend.position = "none")
        })
        names(pie_charts) <- levels(pie_data$Season)
        return(pie_charts)
    }, 
    simplify = FALSE
)
names(pie_charts) <- names(IndVal)



library(ggspatial)
italy <- st_read(paste(HOME_, "Italy_shp/Italy.shp", sep = "/"))

italy["geometry"]
italy$center <- st_centroid(italy$geometry)

lat_long_df <- st_coordinates(italy$center) %>% 
    as.data.frame() %>% 
    setNames(c("lon", "lat")) %>% head(9) 


library(ggtree)
half_size <- 1.2 
for (season in c("Winter", "Spring", "Summer", "Autumn")) {
    p <- ggplot() + 
        geom_point(data = lat_long_df, aes(x = lon, y = lat), fill = "lightgrey", color = "black")
    for (i in seq(1,9)) {
    pie_grob <- ggplotGrob(pie_charts[[i]][[season]])
    p <- p + annotation_custom(
        grob = pie_grob,
        xmin = lat_long_df$lon[i] - half_size, xmax = lat_long_df$lon[i] + half_size,
        ymin = lat_long_df$lat[i] - half_size, ymax = lat_long_df$lat[i] + half_size
    )
}
}


ggsave(file.path(HOME_, "phyto_pie_charts.pdf"), p, width = 10, height = 8, device = pdf, units = "in", dpi = 300)
