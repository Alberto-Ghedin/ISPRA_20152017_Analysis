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

phyto_abund$Basin <- factor(phyto_abund$Basin, levels = ordered_basins, ordered = TRUE)
phyto_abund$Region <- factor(phyto_abund$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
phyto_abund$Season <- factor(phyto_abund$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
abund <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Num_cell_l = sum(Num_cell_l),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season)
    )

# Define colors for regions
colors <- scales::hue_pal()(length(unique(abund$Region)))
palette <- setNames(colors, unique(abund$Region))


# Create the plot
bloom <- data.frame(
    Basin = rep(ordered_basins, each = 4),
    Season = rep(c("Winter", "Spring", "Summer", "Autumn"), times = 6),
    Num_cell_l = rep(c(1e4), 24), 
    Region = rep("BAS", times = 24), 
    label = rep("", times = 24)
)

bloom[which(bloom$Season == "Winter" & bloom$Basin == "NorthAdr"), "label"] <- "Bloom"
bloom
bloom$Region <- factor(bloom$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
bloom$Season <- factor(bloom$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
p <- ggplot(abund, aes(x = Region, y = Num_cell_l, fill = Region)) +
    geom_hline(yintercept = 1e6, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.9) + 
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
    #scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = palette) +
    facet_grid(Season ~ Basin, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "none", 
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        panel.spacing = unit(1, "lines")
    ) +
    scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
    labs(
        y = "Abundance [cells/L]",
        title = "Sample abundance per basin and season"
    ) 
p    
ggsave(file.path(HOME_, "abundance_per_basin_season.pdf"), p, width = 22, height = 13, dpi = 300)

richness <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Taxon = n(),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season)
    )

p <- ggplot(richness, aes(x = Region, y = Taxon, fill = Region)) +
    #geom_hline(yintercept = 1e6, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.9) + 
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
    #scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = palette) +
    facet_grid(Season ~ Basin, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "none", 
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        panel.spacing = unit(1, "lines")
    ) +
    #scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
    #            labels = trans_format("log10", math_format(10^.x))) + 
    labs(
        y = "Abundance [cells/L]",
        title = "Sample abundance per basin and season"
    ) 
p