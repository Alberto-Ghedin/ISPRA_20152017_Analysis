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



# Create the plot
bloom <- data.frame(
    Basin = rep(ordered_basins, each = 4),
    Season = rep(c("Winter", "Spring", "Summer", "Autumn"), times = 6),
    Num_cell_l = rep(c(1e4), 24), 
    Region = rep("BAS", times = 24), 
    label = rep("", times = 24)
)
bloom[which(bloom$Season == "Winter" & bloom$Basin == "NorthAdr"), "label"] <- "Bloom"
bloom$Region <- factor(bloom$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
bloom$Season <- factor(bloom$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)

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

top_classes <- phyto_abund %>% group_by(Class) %>% summarise(Frequency = n_distinct(Date, id)) %>% arrange(desc(Frequency)) %>% dplyr::filter(Class != "nan")

c(top_classes$Class[c(1:4)],"Pyramimonadophyceae", "Dinoflagellata")
abund_groups %>% dplyr::filter(higher_group %in% c(top_classes$Class[c(1:4)],"Pyramimonadophyceae", "Dinoflagellata")) %>% 
mutate(
    Abbr = case_when(
        higher_group == "Pyramimonadophyceae" ~ "PYR",
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
    )
) %>% 
ggplot() + 
geom_boxplot(aes(x = Abbr, y = Abund), width = 0.5) + 
facet_grid(Basin ~ Season) + 
theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "none", 
        strip.text.x = element_text(size = 15, face = "bold"), 
        strip.text.y = element_text(size = 15, face = "bold"), 
        panel.spacing = unit(1, "lines")
    ) +
    #ylim(100, 1e6) + 
    scale_y_continuous(trans="log10", breaks = trans_breaks("log10", function(x) 10^x, breaks = breaks_extended(4)),
                labels = trans_format("log10", math_format(10^.x))) + 
    labs(
        y = "Abundance [cells/L]",
        x = "Phytoplankton groups",
        title = "Sample abundance per basin and season"
    )

bloom <- data.frame(
    Basin = rep(ordered_basins, each = 4),
    Season = rep(c("Winter", "Spring", "Summer", "Autumn"), times = 6),
    Num_cell_l = rep(c(1e4), 24), 
    Region = rep("BAS", times = 24), 
    label = rep("", times = 24)
)
bloom[which(bloom$Season == "Winter" & bloom$Basin == "NorthAdr"), "label"] <- "Bloom"
bloom$Region <- factor(bloom$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
bloom$Season <- factor(bloom$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
group_plot_data <- abund_groups %>% 
mutate(
    Abbr = case_when(
        #higher_group == "Pyramimonadophyceae" ~ "PYR",
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        higher_group == "Unknown" ~ "UNK",
    )
) %>% group_by(Date, id, Abbr) %>%
summarise(
    Abund = sum(Abund),
    Region = first(Region),
    Season = first(Season), 
    Basin = first(Basin), 
    .groups = "drop"
)
p <- group_plot_data %>% 
#dplyr::filter(Abbr != "Unknown") %>%
ggplot() + 
geom_boxplot(aes(x = Abbr, y = Abund), width = 0.5) + 
geom_hline(yintercept = 1e6, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.9) + 
facet_wrap(~Basin, ncol = 2) + 
theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        strip.text = element_text(size = 20),
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
        x = "Phytoplankton groups",
        title = "Abundance of phytoplankton groups in each basin"
    )
p
ggsave(file.path(HOME_, "abundance_per_group.pdf"), p, width = 12, height = 16, dpi = 300)

abund_groups %>% dplyr::filter(higher_group %in% c(top_classes$Class[c(1:4)],"Pyramimonadophyceae", "Dinoflagellata")) %>% 
mutate(
    Abbr = case_when(
        #higher_group == "Pyramimonadophyceae" ~ "PYR",
        higher_group == "Dinoflagellata" ~ "DIN",
        higher_group == "Bacillariophyceae" ~ "DIA",
        higher_group == "Coccolithophyceae" ~ "COC",
        higher_group == "Cryptophyceae" ~ "CRY",
        higher_group != "Unknown" ~ "Else",
        TRUE ~ higher_group
    )
) %>% group_by(Date, id, Abbr) %>%
summarise(
    Abund = sum(Abund),
    Region = first(Region),
    Season = first(Season), 
    Basin = first(Basin), 
    .groups = "drop"
) %>% head()
