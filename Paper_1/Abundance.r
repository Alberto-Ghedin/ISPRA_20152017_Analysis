library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

HOME_ <- "."
phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))

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


abund <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Num_cell_l = sum(Num_cell_l),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season), 
        Longitude = first(Longitude), 
        Latitude = first(Latitude), 
        .groups = "drop"
    )
abund$Region <- factor(abund$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
abund$Basin <- factor(abund$Basin, levels = ordered_basins, ordered = TRUE)

abund %>% dplyr::filter(Region == "SIC") %>% pull(id) %>% unique() %>% sort()
library(rjson)
params <- fromJSON(file = file.path(HOME_, "params.json"))
abund$id <- factor(abund$id, levels = params$ordered_id, ordered = TRUE)


chem_phys <- read.csv("./df_chem_phys.csv")
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))
chem_phys$id <- factor(chem_phys$id, levels = params$ordered_id, ordered = TRUE)


colors <- scales::hue_pal()(length(unique(abund$Region)))
palette <- setNames(colors, unique(abund$Region))
abund %>% mutate(index_id = match(id, params$ordered_id))  %>% 
ggplot(aes(x = id, y = log10(Num_cell_l +1), fill = Region)) +
geom_boxplot() + 
facet_wrap(~Season) +
scale_fill_manual(values = palette)

abund %>% mutate(index_id = match(id, params$ordered_id))  %>% 
ggplot(aes(x = index_id, y = log10(Num_cell_l +1))) +
geom_point(aes(color = Region), shape = 20) + 
scale_color_manual(values = palette) +
geom_smooth(formula = y ~ s(x, bs = "cs", fx = TRUE, k = 30))  


colors <- scales::hue_pal()(length(unique(abund$Region)))
palette <- setNames(colors, unique(abund$Region))
chem_phys %>% ggplot() +
geom_boxplot(aes(x = id, y = Salinity, fill = Region)) +
scale_fill_manual(values = palette) + 
scale_y_reverse()

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
abund_groups$Region <- factor(abund_groups$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
abund_groups$Basin <- factor(abund_groups$Basin, levels = ordered_basins, ordered = TRUE)
abund_groups$Season <- factor(abund_groups$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)


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
        ggplot2::scale_fill_continuous(type = "viridis", limits = fill_limits), 
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
}
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
) %>% group_by(Season, Basin, Group) %>%
summarise(
    Abund = mean(Abund),
    .groups = "drop"
) %>% ggplot(aes(x = Season, y = Group, fill = log10(Abund + 1))) +
geom_tile() +
geom_text(aes(label = sciencific_notation(Abund), colour = ifelse(Abund > 1e4, "black", "white")), size = 8) + 
facet_wrap(~Basin, ncol = 2) + 
labs(title = "Distribution of abundance of main phytoplankton groups", x = "Season", y = "Group") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(NULL, "Abundance [cell/L] (log scale)") 
ggsave(
    file.path(HOME_, "abundance_per_group_heatmap.svg"), 
    p, 
    width = 12, 
    height = 13, 
    dpi = 300
)

fold_change <- abund_groups %>% dplyr::filter(Region %in% c("FVG", "VEN", "EMR")) %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(id, Date, Group) %>% summarise(Abund = sum(Abund), Season = first(Season), Basin = first(Basin), .groups = "drop") %>%
group_by(Season, Basin, Group) %>% 
mutate(seasonal_abund = median(Abund)) %>% 
group_by(Basin, Group) %>% mutate(mean_abund = median(Abund)) %>% 
mutate(fold_change = seasonal_abund / mean_abund) 

abund_percetile <- abund_groups %>% 
mutate(
    Group = case_when(
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK"
    )
) %>% group_by(id, Date, Group) %>% summarise(Abund = sum(Abund), Season = first(Season), Basin = first(Basin), .groups = "drop") %>% 
group_by(Season, Basin, Group) %>% summarise(percetile = quantile(Abund, 0.95), .groups = "drop")

p <- fold_change %>% ggplot(aes(x = Season, y = Group, fill = fold_change)) +
geom_tile() +
geom_text(aes(label = round(fold_change, 2)), #colour = ifelse(Abund > 1e4, "black", "white")), 
size = 8) + 
#facet_wrap(~Basin, ncol = 2) + 
labs(title = "Distribution of abundance of main phytoplankton groups", x = "Season", y = "Group") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(c(0, 5), "Abundance [cell/L] (log scale)") 
p

p <- abund_percetile %>% ggplot(aes(x = Season, y = Group, fill = log10(percetile +1))) +
geom_tile() +
geom_text(aes(label = sciencific_notation(percetile)), #colour = ifelse(Abund > 1e4, "black", "white")),
size = 8) +
facet_wrap(~Basin, ncol = 2) +
labs(title = "Distribution of abundance of main phytoplankton groups", x = "Season", y = "Group") +
ggplot_theme + theme(legend.position = "bottom") +
ggplot_fill_scale(NULL, "Abundance [cell/L] (log scale)")
p

select_group <-  c("DIAT", "DINO", "COCC", "CRYP")
group <- "CRYP"
selected_genera <- phyto_abund %>% mutate(
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
    )) %>% dplyr::filter(Group == group) %>% group_by(Season, Basin, Genus) %>%
summarise(
    Abund = mean(Num_cell_l),
    Basin = first(Basin),
    Season = first(Season),
    .groups = "drop"
) %>% group_by(Basin, Season) %>% 
mutate(
    rel = Abund / sum(Abund)
) %>% mutate(Genus = case_when(Genus == "" ~ "Undet", TRUE ~ Genus)) %>% dplyr::filter(rel > 0.10) %>% pull(Genus) %>% unique()

selected_genera

data_plot <- phyto_abund %>% mutate(
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
    ))%>% dplyr::filter(Group == group) %>% 
mutate(Genus = case_when(Genus == "" ~ "Undet", TRUE ~ Genus)) %>% 
dplyr::filter(Genus %in% selected_genera) %>% group_by(Season, Basin, Genus) %>%
summarise(
    Abund = mean(Num_cell_l),
    Basin = first(Basin),
    Season = first(Season),
    .groups = "drop"
) %>% complete(Genus, nesting(Season, Basin), fill = list(Abund = 0))
data_plot$Season <- factor(data_plot$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
data_plot$Basin <- factor(data_plot$Basin, levels = ordered_basins, ordered = TRUE) 
p <- data_plot %>% ggplot(aes(x = Season, y = Genus, fill = log10(Abund + 1))) +
geom_tile() +
geom_text(aes(label = sciencific_notation(Abund), colour = ifelse(Abund > 1e4, "black", "white")), size = 6) + 
facet_wrap(~Basin, ncol = 2) + 
labs(title = paste("Distribution of abundance of", group, sep = " "), x = "Season", y = "Group") +
ggplot_theme + 
ggplot_fill_scale(NULL, "Abundance [cell/L] (log scale)") + 
theme(
        legend.position = "bottom",       # Moves legend to the bottom
        legend.justification = c(-1.5, 0),   # Aligns legend to the left
        legend.box.just = "left",         # Ensures it stays left-aligned
        legend.margin = margin(t = 10)    # Adds space between plot and legend
    )
p
ggsave(
    file.path(HOME_, paste("abundance_per_group_heatmap_",group,".svg", sep = "")), 
    p, 
    width = 12, 
    height = 8.5, 
    dpi = 300
)



colors <- scales::hue_pal()(length(unique(abund$Region)))
palette <- setNames(colors, unique(abund$Region))
p <- ggplot(abund, aes(x = Region, y = Num_cell_l, fill = Region)) +
    geom_hline(yintercept = 1e6, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.9) + 
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
    #scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = palette) +
    facet_wrap(~Basin, scales = "free_x", ncol = 2) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        strip.text = element_text(size = 23),
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "none", 
        strip.text.x = element_text(size = 15, face = "bold"), 
        strip.text.y = element_text(size = 15, face = "bold"), 
        panel.spacing = unit(1, "lines")
    ) +
    scale_y_continuous(trans="log10", breaks = trans_breaks("log10", function(x) 10^x, breaks = breaks_extended(4)),
                labels = trans_format("log10", math_format(10^.x))) + 
    labs(
        y = "Abundance [cells/L]",
        title = "Sample abundance in each region"
    ) 
p
ggsave(file.path(HOME_, "abundance_per_basin.pdf"), p, width = 12, height = 16, dpi = 300)




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

library(patchwork)
p1 <- species_abund_rich %>% 
ggplot() + 
geom_point(aes(x = n_genus, y = log10(abund_genus + 1), fill = Group),color = "black", shape = 21) +
geom_smooth(aes(x = n_genus, y = log10(abund_genus + 1), color = Group), method = "gam") + 
labs(x = "Genus richness", y = "Abundance [cells/L] (log scale)") +
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
geom_point(aes(x = n_species, y = log10(abund_species  +1), fill = Group), color = "black", shape = 21) + 
geom_smooth(aes(x = n_species, y = log10(abund_species + 1), color = Group), method = "gam") + 
labs(x = "Species richness", y = "") +
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



