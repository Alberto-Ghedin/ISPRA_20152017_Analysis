library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
library(readxl)

IMAGE_FORMAT <- "svg"

ggplot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25, face = "bold"),
    ) 

ggplot_fill_scale <- function(fill_limits, fill_title) {
    z <- list(
        ggplot2::scale_fill_continuous(type = "viridis", limits = fill_limits, oob=scales::squish), 
        ggplot2::scale_colour_manual(values=c("white"="white", "black"="black")), 
        ggplot2::guides(colour = "none"), 
        ggplot2::guides(
        fill = guide_colourbar(
            title = fill_title, 
            title.position = "left", 
            title.theme = element_text(size = 25, face = "bold", margin = margin(r = 3, l = -1), vjust = 1), 
            label.theme = element_text(size = 20),
            barwidth = unit(20, "lines"),
            ticks.linewidth = 1,
            frame.linewidth = 1,
            ticks.colour = "black",
            frame.colour  ='black'
            )
        )
    )
    return(z)     
}

sciencific_notation <- function(x) {
    y <- gsub("\\+0", "", scales::scientific(x, digits = 2))
    y <- gsub("0\\.0e0", "0", y)
    return(y)
}

HOME_ <- "."
phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))
phyto_abund <- phyto_abund %>% mutate(
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
)
phyto_abund$New_basin <- factor(phyto_abund$New_basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)

from_region_to_abreviation <- c(
    "Friuli-Venezia-Giulia" = "FVG",
    "Veneto" = "VEN", 
    "Emilia-Romagna" = "EMR",
    "Marche" = "MAR",
    "Abruzzo" = "ABR",
    "Molise" = "MOL",
    "Puglia" = "PUG",
    "Basilicata" = "BAS",
    "Calabria" = "CAL",
    "Sicilia" = "SIC",
    "Campania" = "CAM", 
    "Lazio" = "LAZ",
    "Toscana" = "TOS",
    "Liguria" = "LIG",
    "Sardegna" = "SAR"
)
ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")

ordered_transect <- c(
    "SMTS", "SMLG", "VENEZIA", "ROSOLINA", "PORTO_GARIBALDI", "CESENATICO", "RIMINI", "Chienti", "Esino", "GU",
    "VA", "R14001_B2", "FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO", "PORTO_CESAREO", "PUNTA_RONDINELLA", "SINNI", "Villapiana",
    "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche", "Isole_Ciclopi", "Plemmirio", "Isola_Correnti", "San_Marco", "Isole_Egadi", "Capo_Gallo","Vibo_marina", 
    "Cetraro", "Cilento", "Salerno", "Napoli", "Domizio", "m1lt01", "m1lt02",  "m1rm03",  "m1vt04", "Collelungo","Carbonifera",
    "Donoratico", "Fiume_Morto", "Mesco" , "Portofino"  ,  "Voltri"  , "Quiliano" , 
    "Olbia", "Arbatax", "Villasimius", "Cagliari", "Oristano", "Alghero", "Porto_Torres"

)
sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
library(rjson)
params <- fromJSON(file = file.path(HOME_, "params.json"))

abund <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Abund = sum(Num_cell_l),
        Basin = first(Basin),
        New_basin = first(New_basin),
        Region = first(Region), 
        Season = first(Season), 
        Closest_coast = first(Closest_coast),
        Longitude = first(Longitude), 
        Latitude = first(Latitude), 
        .groups = "drop"
    )
abund$Region <- factor(abund$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
abund$Basin <- factor(abund$Basin, levels = ordered_basins, ordered = TRUE)
abund$id <- factor(abund$id, levels = params$ordered_id, ordered = TRUE)
abund <- merge(
    abund, 
    sea_depth %>% select(id, SeaDepth, Transect)
)
abund$Transect <- factor(abund$Transect, levels = ordered_transect, ordered = TRUE)
abund$New_basin <- factor(abund$New_basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)
abund <- abund %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) 


one_every_n_item <- function(list, n) {
    return(
        list[seq(1, length(list), n)]
    )
}


ordered_latitude <- abund %>% dplyr::select(id, Latitude) %>% distinct() %>% 
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Latitude) %>% round(2)
ordered_longitude <- abund %>% dplyr::select(id, Longitude) %>% distinct() %>%
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Longitude) %>% round(2)
colors <- scales::hue_pal()(length(unique(abund$New_basin)))
#sorting to be cosistent with color legend of the map 
palette <- setNames(colors, sort(as.character(phyto_abund$New_basin %>% unique())))
p <- abund %>% 
mutate(index_id = match(id, params$ordered_id)) %>% 
ggplot(aes(x = index_id, y = log10(Abund +1))) + 
geom_boxplot(aes(group = id, fill = New_basin)) +
scale_fill_manual(values = palette) + 
labs(title = "Sample abundance across all stations", y = "Abundance [cells/L] (log scale)") +
scale_x_continuous(
    name = "Latitude",
    breaks =  seq(1, length(params$ordered_id), 3), 
    labels = one_every_n_item(ordered_latitude, 3),
    limits = c(-0.01, 162.01), 
    sec.axis = sec_axis(
      ~., 
      name = "Longitude",
      breaks = seq(1, length(params$ordered_id), 3),
      labels = one_every_n_item(ordered_longitude, 3)
    )
) + 
ggplot2::theme(
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 17),
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 17),
        axis.text.y = element_text(angle = 0, hjust = 0, size = 17),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18, face = "bold"),
    ) 
ggsave(
    file.path(HOME_, "abundance_per_basin.svg"), 
    p, 
    width = 18, 
    height = 8.5, 
    dpi = 300
)

kruskal_abund_test <- sapply(
    abund %>% pull(New_basin) %>% unique(),
    function(x) {
        test <- kruskal.test(Abund ~ Season, abund %>% dplyr::filter(New_basin == x))
        return(test)
    }, 
    simplify = FALSE
)
names(kruskal_abund_test) <- abund %>% pull(New_basin) %>% unique()

data.frame(
    New_basin = names(kruskal_abund_test),
    p.value = sapply(kruskal_abund_test, function(x) x$p.value)
) %>% write.csv(file.path(HOME_, "kruskal_over_season_abund_test.csv"), row.names = FALSE)

basin <- "SIC"
abund %>% dplyr::filter(New_basin == basin) %>%
ggplot(aes(x = id, y = log10(Abund + 1))) +
geom_boxplot(aes(fill = Season)) 

pairwise.wilcox.test(abund %>% dplyr::filter(New_basin == basin) %>% pull(Abund), 
abund %>% dplyr::filter(New_basin == basin) %>% pull(Season), p.adjust.method = "BH")



abund_groups <- phyto_abund %>% mutate(
    higher_group = case_when(
        Class == "Dinoflagellata incertae sedis" ~ "Dinoflagellata",
        Taxon == "Noctilucea" ~ "Dinoflagellata",
        Class == "nan" ~ "Unknown", 
        Class == "Dinophyceae" ~ "Dinoflagellata", 
        TRUE ~ Class
    )
) %>% group_by(Date, id, higher_group) %>%
summarise(
    Abund = sum(Num_cell_l),
    Region = first(Region),
    Season = first(Season), 
    Basin = first(Basin), 
    .groups = "drop"
)
abund_groups <- abund_groups %>% mutate(
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
)
abund_groups$Region <- factor(abund_groups$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
abund_groups$Basin <- factor(abund_groups$Basin, levels = ordered_basins, ordered = TRUE)
abund_groups$Season <- factor(abund_groups$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
abund_groups$New_basin <- factor(abund_groups$New_basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)
abund_groups <- abund_groups %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30"))

library(ggh4x)
library(legendry)
p <- abund_groups %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(Season, New_basin, Group) %>%
summarise(
    Abund = quantile(Abund, 0.75),
    .groups = "drop"
) %>% complete(Group, nesting(Season, New_basin), fill = list(Abund = 0)) %>% 
ggplot(aes(x = Group, y = Season, fill = log10(Abund +1))) +
geom_tile(color = "gray", linewidth = 0.5) +
geom_text(aes(label = sciencific_notation(Abund), colour = ifelse(Abund > 1e4, "black", "white")),
size = 8) +
scale_y_discrete(limits = rev) + 
facet_grid(New_basin ~ ., scale = "free_y") +
labs(title = "Distribution of abundance of main phytoplankton groups", y = "Season", x = "Group") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(c(2, 6), "Abundance [cell/L] (log scale)") 
p
ggsave(
    file.path(HOME_, "abundance_per_group_heatmap.svg"), 
    p, 
    width = 9, 
    height = 16, 
    dpi = 300
)

p <- abund_groups %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(Season, New_basin, Group) %>%
summarise(
    Abund = quantile(Abund, 0.75),
    .groups = "drop"
) %>% ggplot(aes(x = Group, y = Season, size = log10(Abund +1), color = log10(Abund +1))) +
geom_point()+
scale_size(range = c(1, 6), "Abundance [cell/L] (log scale)") + 
scale_y_discrete(limits = rev) + 
facet_wrap(~New_basin, ncol = 1) +
labs(title = "Distribution of abundance of main phytoplankton groups", y = "Season", x = "Group") +
ggplot_theme + theme(legend.position = "bottom")# +
p
ggsave(
    file.path(HOME_, "abundance_per_group_bubbleplot.svg"), 
    p, 
    width = 9, 
    height = 16, 
    dpi = 300
)


p <- abund_groups %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(Season, New_basin, Group) %>%
summarise(
    Abund = quantile(Abund, 0.75),
    .groups = "drop"
) %>% group_by(Season, New_basin) %>% 
mutate(Rel_abund = Abund / sum(Abund)) %>% 
ggplot(aes(y = Season, x = Rel_abund, fill = Group)) +
geom_bar(stat = "identity") +
facet_wrap(~New_basin, ncol = 1) +
scale_y_discrete(limits = rev) + 
labs(title = "Distribution of abundance of main phytoplankton groups", y = "Season", x = "Relative abundance") +
ggplot_theme 
ggsave(
    file.path(HOME_, "abundance_per_group_barplot.svg"), 
    p, 
    width = 9, 
    height = 16, 
    dpi = 300
)


p <- ggplot(abund) + 
geom_point(aes(x = Longitude, y = Num_cell_l, color = Region)) +
    scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
ggsave(file.path(HOME_, "abundance_per_longitude.pdf"), p, width = 22, height = 13, dpi = 300)

p <- ggplot(abund) + 
geom_point(aes(x = Latitude, y = Num_cell_l, color = Region)) +
    scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
ggsave(file.path(HOME_, "abundance_per_latitude.pdf"), p, width = 22, height = 13, dpi = 300)




abund_only_genera <- phyto_abund %>% dplyr::mutate(
    Genus = case_when(
        Class == "nan" ~ Taxon,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
    ) %>% group_by(Date, id, Genus) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    Basin = first(New_basin),
    Season = first(Season),
    Class = first(Class),
    .groups = "drop"
) 



sheets <- excel_sheets("./indval_only_genera_per_basin.xlsx")
all_data <- sapply(sheets, function(sheet) {
    data <- read_excel("./indval_only_genera_per_basin.xlsx", sheet = sheet, col_names = TRUE)
    colnames(data)[1] <- "Taxon"
    data <- data %>% dplyr::filter(Taxon != "Other phytoplankton")
    return(data)
}, simplify = FALSE)
names(all_data) <- sheets


quantile_abund_genera <- abund_only_genera %>% group_by(Date, id, Basin, Class) %>% 
mutate(sample_abund = sum(Abund)) %>% group_by(Date, id, Genus) %>% #mutate(rel_abund = Abund / sample_abund) %>% 
group_by(Season, Basin, Genus) %>%
summarise(
    Abund =quantile(Abund, 0.75),
    Class = first(Class),
    .groups = "drop"
)



basin <- "NA"
selected_genera <- all_data[[basin]] %>% rowwise() %>%
  dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
  ungroup() %>% pull(Taxon)

quantile_abund_genera %>% head()
quantile_abund_genera %>% dplyr::filter(Genus %in% selected_genera & Basin == basin) %>% head() %>% 
complete(Genus, nesting(Season, Basin)) # %>% 
ggplot(aes(x = Season, y = Genus, fill = log10(Abund +1))) +
geom_tile() + 
geom_text(aes(label = sciencific_notation(Abund), colour = ifelse(Abund > 1e4, "black", "white")),
size = 8) +
#facet_grid(Class ~ ., scale = "free_y") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(NULL, "Abundance [cell/L] (log scale)") 


basin <- "CA"
dominant_species_path <- paste(HOME_, "Dominant_species", sep = "/")
dir.create(dominant_species_path, showWarnings = FALSE)

sapply(
    names(all_data),
    function(basin) {
        selected_genera <- all_data[[basin]] %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
        ungroup() %>% pull(Taxon)
        dominant_genera <- quantile_abund_genera %>% dplyr::filter(Genus %in% selected_genera & Basin == basin) %>% 
            group_by(Season, Basin) %>% 
            mutate(
                rel_abund = Abund / sum(Abund), 
                .groups = "drop"
            ) %>% dplyr::filter(rel_abund > 0.02) %>% pull(Genus) %>% unique()
        p <-  quantile_abund_genera %>% dplyr::filter(Genus %in% dominant_genera & Basin == basin) %>% 
            mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
            group_by(Season, Basin) %>% 
            mutate(
                rel_abund = Abund / sum(Abund), 
                .groups = "drop"
            ) %>% 
            ggplot(aes(x = Season, fill = Genus, y = rel_abund)) +
            geom_bar(stat = "identity", linewidth = 0.5, color = "black") + 
            labs(y = "Relative abundance", x = "Season") +
            ggplot_theme + theme(legend.position = "right") + 
            scale_fill_manual(values = unname(colorBlindness::paletteMartin))
        ggsave(
            file.path(dominant_species_path, basin),
            device = IMAGE_FORMAT, 
            p, 
            width = 12, 
            height = 8, 
            dpi = 300
        )
    }
)
selected_genera <- all_data[["SAR"]] %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
        ungroup() %>% pull(Taxon)
dominant_genera <- quantile_abund_genera %>% dplyr::filter(Genus %in% selected_genera & Basin == "SAR") %>% 
            group_by(Season, Basin) %>% 
            mutate(
                rel_abund = Abund / sum(Abund), 
                .groups = "drop"
            ) %>% dplyr::filter(rel_abund > 0.02) %>% pull(Genus) %>% unique()
dominant_genera
ibrary(colorBlindness)



unname(colorBlindness::paletteMartin)

fold_change <- abund_groups %>% #dplyr::filter(Region %in% c("FVG", "VEN", "EMR")) %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(id, Date, Group) %>% summarise(Abund = sum(Abund), Season = first(Season), New_basin = first(New_basin), .groups = "drop") %>%
complete(Group, nesting(Season, New_basin), fill = list(Abund = 0)) %>% 
group_by(New_basin, Group) %>% mutate(annual_mean = mean(Abund), annual_median = median(Abund), annual_quantile = quantile(Abund, 0.75), .groups = "drop") %>% 
group_by(Season, New_basin, Group) %>% 
summarise(
    seasonal_mean = mean(Abund), seasonal_median = median(Abund), seasonal_quantile = quantile(Abund, 0.75),
    annual_mean = first(annual_mean), annual_median = first(annual_median), annual_quantile = first(annual_quantile), 
    .groups = "drop"
    ) %>% 
mutate(fold_mean = seasonal_mean / annual_mean, fold_median = seasonal_median / annual_median, fold_quantile = seasonal_quantile / annual_quantile)



p <- fold_change %>% ggplot(aes(x = Season, y = Group, fill = fold_mean)) +
geom_tile() +
geom_text(aes(label = round(fold_mean, 2), colour = ifelse(fold_mean > 1, "black", "white")), 
size = 8) + 
facet_wrap(~New_basin, ncol = 2) + 
labs(title = "Distribution of abundance of main phytoplankton groups", x = "Season", y = "Group") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(NULL, "Abundance [cell/L] (log scale)") 
ggsave(
    file.path(HOME_, "abundance_per_group_heatmap_fold_mean.svg"), 
    p, 
    width = 12, 
    height = 13, 
    dpi = 300
)










library(patchwork)

species_abund_rich <- phyto_abund %>% mutate(
    Group = case_when(
        Class == "Dinoflagellata incertae sedis" ~ "DINO",
        Taxon == "Noctilucea" ~ "DINO",
        Class == "nan" ~ "UNK", 
        Class == "Dinophyceae" ~ "DINO", 
        Class == "Bacillariophyceae" ~ "DIAT",
        Class == "Coccolithophyceae" ~ "COCC",
        Class == "Cryptophyceae" ~ "CRYP",
        Class != "Unknown" ~ "ElSE",
        Class == "Unknown" ~ "UNK"
    )) %>% dplyr::filter(Group %in% c("DIAT", "DINO") & Det_level %in% c("Genus", "Species")) %>% 
    group_by(Date, id, Group) %>% 
    summarise(
    n_genus = n_distinct(Genus),
    n_species = n_distinct(Taxon[Det_level == "Species"]), 
    abund_genus = sum(Num_cell_l),
    abund_species = sum(Num_cell_l[Det_level == "Species"]),
    .groups = "drop"
  )


p1 <- species_abund_rich %>% 
ggplot() + 
geom_point(aes(y = n_genus, x = log10(abund_genus + 1), fill = Group),color = "black", shape = 21) +
geom_smooth(aes(y = n_genus, x = log10(abund_genus + 1), color = Group), method = "gam") + 
labs(y = "Genus richness", x = "Abundance [cells/L] (log scale)") +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    plot.title = element_text(size = 25, hjust = 0.5),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18)
)
p2 <- species_abund_rich %>%
ggplot() +
geom_point(aes(y = n_species, x = log10(abund_species  +1), fill = Group), color = "black", shape = 21) + 
geom_smooth(aes(y = n_species, x = log10(abund_species + 1), color = Group), method = "gam") + 
labs(y = "Species richness", x = "") +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    plot.title = element_text(size = 25, hjust = 0.5),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18)
) 
ylims <- range(
    ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,
    ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range
    )
p <- p1 + ylim(ylims) + 
p2 + ylim(ylims) + 
plot_annotation(title = 'Richness-abundance curve for diatoms and Dinoflagellata') + 
plot_layout(guides = "collect") & theme(
    legend.position = "right", 
    plot.title = element_text(size = 25, hjust = 0.5)
    )

ggsave(
    file.path(HOME_, "richness_abundance_curve.svg"), 
    p, 
    width = 12, 
    height = 6, 
    dpi = 300
)



