library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ComplexHeatmap)
library(readxl)
library(tibble)
library(gridExtra)
library(cowplot)
library(tidytext)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(grid)
library(tidyverse)
library(paletteer)

HOME_ <- "./Paper_1"
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




samples <- phyto_abund %>%
    select(Region, Date, id) %>%
    distinct()


data <- samples %>%
    group_by(Region, Date) %>%
    summarise(id = n_distinct(id)) %>%
    ungroup() %>%
    complete(Region, Date, fill = list(id = 0)) %>%
    group_by(Region) %>%
    mutate(id = id / n_distinct(samples$id[samples$Region == Region])) %>%
    spread(Date, id, fill = 0)

data <- data %>%
    column_to_rownames("Region") %>%
    as.matrix()

annot_label <- ifelse(data == 0 | data == 1, "", sprintf("%.2f", data))
data[data > 0 & data < 1] <- 0.5
colnames(data) <- sapply(strsplit(colnames(data), "-"), function(x) paste(rev(x[1:2]), collapse = "-"))
cmap <- colorRamp2(c(1,0.5,0), c("chartreuse4", "#FFDE3B", "red"))

pdf(file.path(HOME_, "heatmap_samples_per_region.pdf"), width = 22, height = 13)
lgd <- Legend(at = 1:3, 
labels = c("Sampled", "Partially sampled", "Not sampled"), 
title = "", 
legend_gp = gpar(fill = cmap(c(1,0.5,0))),
labels_gp = gpar(fontsize = 18),
ncol = 3, 
gap = unit(1, "cm"), 
border = "black"
)

heatmap <- Heatmap(data, 
                   column_title = "Sampling effort per region",
                   column_title_gp = gpar(fontsize = 25, fontface = "bold"),
                   col = cmap, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   show_column_names = TRUE, 
                   show_row_names = TRUE, 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x, y, width, height, gp = gpar(col = "black", lty = "solid", lwd = 2))
                       grid.text(annot_label[i, j], x, y, gp = gpar(col = "black", fontsize = 18))
                   },
                   border_gp = gpar(col = "black", lty = "solid", linewidth = 2),
                   row_names_centered = FALSE,
                   column_names_rot = -45, 
                   row_names_gp = gpar(fontsize = 15),
                   column_names_gp = gpar(fontsize = 15), 
                   show_heatmap_legend = FALSE
                   )
draw(heatmap, annotation_legend_list = list(lgd), annotation_legend_side = "bottom")
dev.off()





richness <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Taxon = n(),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season)
    )

p <- ggplot(richness, aes(x = Region, y = Taxon, fill = Region)) +
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
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
    labs(
        y = "Taxa richness",
        title = "Sample richness per basin and season"
    ) 
ggsave(file.path(HOME_, "richness_per_basin_season.pdf"), p, width = 22, height = 13, dpi = 300)

richness_species <- phyto_abund %>% dplyr::filter(Det_level == "Species") %>%
    group_by(Date, id) %>%
    summarise(
        Taxon = n(),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season)
    )

richness_species
richness_genera <- phyto_abund %>% dplyr::filter(Det_level == "Genus") %>%
    group_by(Date, id) %>%
    summarise(
        Taxon = n(),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season)
    )





Genus_abund <- phyto_abund %>% group_by(id, Date, Genus) %>% 
    summarise(Num_cell_l = sum(Num_cell_l), 
                Basin = first(Basin), 
                Region = first(Region), 
                Season = first(Season)
                )

top_genera <- top_genera[c(1:10),] %>% pull(Genus)
Genus_abund %>% filter(Genus %in% top_genera) %>% 
    group_by(Genus, Basin, Season) %>% 
    summarise(Abund = median(Num_cell_l)) %>%
    ggplot(aes(x = Season, y = Abund, group = Genus, fill = Genus, color = Genus)) +
    facet_wrap(~Basin) +
    scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
               labels = trans_format("log10", math_format(10^.x)), limits = c(10, 4 * 10^4)) + 
    geom_point(size = 4) + 
    geom_line(aes(linetype = Genus), size = 1.5) + 
    scale_x_discrete(expand = c(0,0.2))




top_taxon <- top_taxa[c(1:10),] %>% pull(Taxon)
Taxon_abund <- phyto_abund %>% group_by(id, Date, Taxon) %>% 
    summarise(Num_cell_l = sum(Num_cell_l), 
                Basin = first(Basin), 
                Region = first(Region), 
                Season = first(Season)
                )
Taxon_abund %>% filter(Taxon %in% top_taxon) %>% 
    group_by(Taxon, Basin, Season) %>% 
    summarise(Abund = median(Num_cell_l)) %>%
    ggplot(aes(x = Season, y = Abund, group = Taxon, fill = Taxon, color = Taxon)) +
    facet_wrap(~Basin) +
    scale_y_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),
               labels = trans_format("log10", math_format(10^.x))) + 
    geom_point(size = 4) + 
    geom_line(aes(linetype = Taxon), size = 1.5) #+ 
    #scale_x_discrete(expand = c(0,0.2))



