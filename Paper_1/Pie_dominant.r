library(openxlsx)
library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(tidyverse)

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


abund_only_genera <- phyto_abund %>% 
    dplyr::filter(Det_level %in% c("Genus", "Species")) %>%
    group_by(Date, id, Basin, Genus) %>% 
    summarise(Abund = sum(Num_cell_l), Season = first(Season), Class = first(Class), .groups = "drop")


quantile_abund_genera <- abund_only_genera %>% 
    group_by(Season, Basin, Genus) %>%
    summarise(
        Abund =mean(Abund),
    .groups = "drop"
)

relevant_genera <- sapply(
    IndVal, 
    function(df) {
        genera <- df %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
        pull(Taxon)
        return(genera)
    }
)
names(relevant_genera) <- names(IndVal)

pie_charts <- sapply(
    names(IndVal),
    function(basin) {
        genera_season <- IndVal[[basin]] %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
        pivot_longer(
            cols = where(is.numeric),
            names_to = "Season",
            values_to = "IndVal"
        ) %>% 
        group_by(Taxon) %>% 
        summarise(
        Season = Season[which.max(IndVal)],
        ) %>% rename(Genus = "Taxon")


        pie_data <- quantile_abund_genera %>% dplyr::filter(Basin == basin) %>%
                merge(
                    genera_season, 
                    by = c("Genus", "Season")
                )  %>%
                    mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
                    group_by(Season) %>%
                    mutate(rel_abund = Abund / sum(Abund)) %>% 
                    dplyr::filter(rel_abund > 0.05)

        pie_data <- abund_only_genera %>% dplyr::filter(Genus %in% relevant_genera[[basin]] & Basin == basin) %>% 
                        dplyr::select(Date, id, Genus, Season, Abund) %>%
                        complete(Genus, nesting(Date, id, Season), fill = list(Abund = 0)) %>% 
                        group_by(
                            Season, Genus
                        ) %>% 
                        summarise(
                            abund = mean(Abund)
                        ) %>% mutate(
                            rel_abund = abund / sum(abund)
                        ) %>%
                        group_by(Season) %>%
                        arrange(desc(rel_abund)) %>%
                        mutate(cum_abund = cumsum(rel_abund)) %>%
                        dplyr::filter(cum_abund <= 0.99 | row_number() == 1) %>% 
                        mutate(
                            rel_abund = abund / sum(abund)
                        ) %>%
                        ungroup()
        pie_data <- pie_data %>% mutate(
            Genus = case_when(
                Genus == "Thalassionema" ~ "Tham", 
                Genus == "Thalassiosira" ~ "Thar",
                TRUE ~ Genus
            )
        )
        pie_data$Genus <- substr(pie_data$Genus, 1, 4)
        pie_charts <- lapply(unique(pie_data$Season), function(season) {

            season_data <- pie_data %>% dplyr::filter(Season == season) %>% 
            head(10)
            label_pos <- season_data %>% 
            mutate(csum = rev(cumsum(rev(abund))), 
                   pos = abund/2 + lead(csum, 1),
                   pos = if_else(is.na(pos), abund/2, pos))

            ggplot(season_data, aes(x = "", y = abund, fill = fct_inorder(Genus))) +
                geom_bar(stat = "identity", width = 1, color = "black") +
                coord_polar(theta = "y") +
                labs(x = NULL, y = NULL) +
                geom_label_repel(
                    data = label_pos,
                    aes(y = pos, label = Genus),
                    size = 4.5, 
                    nudge_x = 1, 
                    fill = "white", 
                    color = "black",
                    show.legend = FALSE, 
                    force = 100
                    ) + 
                theme_void() +
                #scale_fill_grey(start = 0.3, end = 0.9) +
                scale_fill_manual(values = unname(colorBlindness::paletteMartin[-1])) + 
                theme(legend.position = "none")
        })
        names(pie_charts) <- unique(pie_data$Season)
        return(pie_charts)
    }, 
    simplify = FALSE
)
names(pie_charts) <- names(IndVal)


library(ggspatial)
library(sf)
italy <- st_read(paste(HOME_, "Italy_shp/Italy.shp", sep = "/"))
surroundings <- st_read(paste(HOME_, "Surrounding_shp/Surrounding.shp", sep = "/"))

lat_long_df <- phyto_abund %>% distinct(id, Longitude, Latitude, Basin) %>% 
group_by(Basin) %>%
summarise(
    lon = mean(Longitude, na.rm = TRUE),
    lat = mean(Latitude, na.rm = TRUE)
) %>% column_to_rownames(var = "Basin")

basins <- c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR")
half_size <- 1.2 

for (season in c("Winter", "Spring", "Summer", "Autumn")) {

    p <- ggplot() +
      geom_sf(data = surroundings, fill = "grey", color = "black") +
      geom_sf(data = italy, fill = "lightgrey", color = "black") +
      scale_x_continuous(breaks = seq(8, 18.5, 1.5), limits = c(8, 18.5)) +
      scale_y_continuous(breaks = seq(37, 46, 1), limits = c(37, 46)) +
      xlim(8, 18.5) +
      ylim(37, 46) +
      ggspatial::annotation_north_arrow(
        location = "tr", which_north = "true",
        pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
        height = unit(3, "cm"), width = unit(3, "cm"),
        style = ggspatial::north_arrow_nautical(
          fill = c("grey40", "white"),
          line_col = "grey20"
        )
      ) +
      annotation_scale(location = "bl", width_hint = 0.2) +
      labs(x = "Longitude", y = "Latitude", title = paste("Characteristic genera across Italy in", season, sep = " ")) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )


       for (i in basins) {
        if (!is.null(pie_charts[[i]][[season]])) {
            pie_grob <- ggplotGrob(pie_charts[[i]][[season]])
            p <- p + annotation_custom(
            grob = pie_grob,
            xmin = lat_long_df[i, "lon"] - half_size, xmax = lat_long_df[i, "lon"] + half_size,
            ymin = lat_long_df[i, "lat"] - half_size, ymax = lat_long_df[i, "lat"] + half_size
            )
        }
       }
    ggsave(file.path(HOME_, paste0("phyto_pie_charts_", season, ".pdf")), p, width = 18, height = 13, device = pdf, units = "in", dpi = 300)
}


basin <- "NA"
pie_data <- abund_only_genera %>% dplyr::filter(Genus %in% relevant_genera[[basin]] & Basin == basin) %>% 
                        dplyr::select(Date, id, Genus, Season, Abund) %>%
                        complete(Genus, nesting(Date, id, Season), fill = list(Abund = 0)) %>% 
                        group_by(
                            Season, Genus
                        ) %>% 
                        summarise(
                            abund = mean(Abund)
                        ) %>% mutate(
                            rel_abund = abund / sum(abund)
                        ) %>%
                        group_by(Season) %>%
                        arrange(desc(rel_abund)) %>%
                        mutate(cum_abund = cumsum(rel_abund)) %>%
                        dplyr::filter(cum_abund <= 0.99 | row_number() == 1) %>% 
                        mutate(
                            rel_abund = abund / sum(abund)
                        ) %>%
                        ungroup()
pie_data %>% arrange(desc(abund)) %>% dplyr::filter(Season == "Winter") %>% 
mutate(c2 = lag(cum_abund, 1))
pie_data$Genus <- substr(pie_data$Genus, 1, 4)

season <- "Winter"

        names(pie_charts) <- unique(pie_data$Season)

season_data
label_pos
df <- data.frame(value = c(15, 25, 32, 28),
                 group = paste0("G", 1:4))
df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df, aes(x = "" , y = value, fill = group)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()

df
df2
