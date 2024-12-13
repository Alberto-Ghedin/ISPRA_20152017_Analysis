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

top_taxa <- read.csv("./Taxa_freq_95.csv")
top_species <- read.csv("./Species_freq_top.csv")
top_genera <- read.csv("./Genera_freq_top.csv")
top_classes <- read.csv("./Classes_freq_top.csv") %>% dplyr::filter(Class != 'nan')

select_and_order <- function(data, title, n_otu = 10, n_samples = 2220) {
    data <- data[1:n_otu, ]
    data$Frequency <- round(data[, 2] / n_samples, 2)
    data$Taxon <- factor(data[, 1], levels = data[, 1][order(data$Frequency, decreasing = FALSE)])
    data$type <- title
    return(data %>% dplyr::select(Taxon, Frequency, type))
}

data <- list(
    select_and_order(top_taxa, "Most common taxa", n_otu = 10, n_samples = 2220),
    select_and_order(top_classes, "Most common Classes", n_otu = 10, n_samples = 2220),
    select_and_order(top_genera, "Most common Genera", n_otu = 10, n_samples = 2220),
    select_and_order(top_species, "Most common Species", n_otu = 10, n_samples = 2220)
) %>% bind_rows()

data$type <- factor(data$type, levels = c("Most common taxa", "Most common Classes", "Most common Genera", "Most common Species"))
data <- data %>% arrange(type, desc(Frequency))
data %>%  mutate(Taxon = reorder_within(Taxon, Frequency, within = type))
plot_most_common <- function(data) {    
    p <- ggplot(data, aes(x = Frequency, y = Taxon)) +
        geom_point(size = 5) +
        labs(x = "Frequency of occurrence", y = "") +
        theme(
            axis.text.x = element_text(size = 20, color = "black"),
            axis.text.y = element_text(size = 20, face = "italic", color = "black"),
            axis.title.x = element_text(size = 20), 
            axis.title.y = element_text(size = 20),
            strip.text = element_text(size=22, face = "bold"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
            strip.background = element_rect(color = "black", linewidth = 1)
        ) +
        scale_x_continuous(labels = scales::percent) + 
        tidytext::scale_y_reordered() +
        facet_wrap(~type, ncol = 2, scales = "free", labeller = label_wrap_gen(multi_line = FALSE))
    
    return(p)
}

pdf(file.path(HOME_, "most_common_taxa_classes_genera_species.pdf"), width = 24, height = 14)
plot_most_common(data %>%  mutate(Taxon = reorder_within(Taxon, Frequency, within = type)))
grid::grid.text(label = "(a)",
          x=0.13, 
          y=0.98,
          gp=gpar(fontsize=22, col="black"))
grid::grid.text(label = "(c)",
          x=0.13, 
          y=0.50,
          gp=gpar(fontsize=22, col="black"))
grid::grid.text(label = "(b)",
          x=0.62, 
          y=0.98,
          gp=gpar(fontsize=22, col="black"))
grid::grid.text(label = "(d)",
          x=0.62, 
          y=0.50,
          gp=gpar(fontsize=22, col="black"))
dev.off()

rich_classes <- phyto_abund %>%
    filter(Taxon != "Other phytoplankton" & Det_level == "Species") %>%
    group_by(Class) %>%
    summarise(n_distinct = n_distinct(Taxon)) %>%
    arrange(desc(n_distinct)) %>% rename(Taxon = Class)

rich_genera <- phyto_abund %>%
    filter(Taxon != "Other phytoplankton" & Det_level == "Species") %>%
    group_by(Genus) %>%
    summarise(n_distinct = n_distinct(Taxon)) %>%
    arrange(desc(n_distinct)) %>% rename(Taxon = Genus)

select_and_order <- function(data, title, n_otu = 10) {
    data <- data[1:n_otu, ]
    data$type <- title
    return(data %>% dplyr::select(Taxon, n_distinct, type))
}

data <- list(
    select_and_order(rich_classes, "Top species-rich classes", n_otu = 10), 
    select_and_order(rich_genera, "Top species-rich genera", n_otu = 10)
) %>% bind_rows()
data <- data %>%  mutate(Taxon = reorder_within(Taxon, n_distinct, within = type))

plot_most_common <- function(data) {    
    p <- ggplot(data, aes(x = n_distinct, y = Taxon)) +
        geom_point(size = 5) +
        labs(x = "Number of species", y = "") +
        theme(
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, face = "italic", color = "black"),
            axis.title.x = element_text(size = 15), 
            axis.title.y = element_text(size = 15),
            strip.text = element_text(size=18, face = "bold"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
            strip.background = element_rect(color = "black", linewidth = 1)
        ) +
        tidytext::scale_y_reordered() +
        facet_wrap(~type, ncol = 2, scales = "free", labeller = label_wrap_gen(multi_line = FALSE))
    
    return(p)
}

pdf(file.path(HOME_, "top_species_rich_classes_genera.pdf"), width = 12, height = 6)
plot_most_common(data %>%  mutate(Taxon = reorder_within(Taxon, n_distinct, within = type)))
grid::grid.text(label = "(a)",
          x=0.11, 
          y=0.98,
          default.unit = "npc",
          gp=gpar(fontsize=20, col="black"))
grid::grid.text(label = "(b)",
          x=0.59, 
          y=0.98,
          gp=gpar(fontsize=20, col="black"))
dev.off()

abund <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Num_cell_l = sum(Num_cell_l),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season), 
        Longitude = first(Longitude), 
        Latitude = first(Latitude)
    )

cat_contribution <- phyto_abund %>% 
    group_by(Region, Det_level) %>% 
    summarise(
        Abund = sum(Num_cell_l)
    ) %>% group_by(Region) %>%
    mutate(
        Rel_Abund = Abund / sum(Abund)
    )

cat_contribution <- phyto_abund %>% 
    group_by(Region, Date, id, Det_level) %>% 
    summarise(Cat_abund = sum(Num_cell_l)) %>%
    group_by(Region,Det_level) %>%
    summarise(mean_Abund =mean(Cat_abund)) %>%  group_by(Region, Det_level) %>% 
    summarise(
        tot_abund = sum(mean_Abund)
    ) %>% group_by(Region) %>%
    mutate(
        rel_cont = tot_abund / sum(tot_abund)
    )

cat_contribution$Det_level <- factor(cat_contribution$Det_level, levels = c("Species", "Genus", "Higher cat.", "Unknown"))
cat_contribution$Region <- factor(cat_contribution$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
p <- ggplot(cat_contribution, aes(x = Region, y = rel_cont, fill = Det_level)) +
    #scale_y_continuous(labels = scales::percent) + 
    scale_fill_manual(values = c("Species" = "#1B9E77", "Genus" = "#D95F02", "Higher cat." = "#7570B3", "Unknown" = "black")) +
    geom_bar(stat = "identity", color = "black", linewidth = 1) +
    labs(x = "Region", y = "Proportion of abundance", fill = "Identification \n level") +
    theme_minimal() +
    ggtitle("Average contribuion of each identification level to the sample abundance in each region") +
    #scale_fill_manual()
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25, face = "bold"),
        #legend.position = "none", 
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    ) 
p
ggsave(file.path(HOME_, "lebundance_per_region.pdf"), p, width = 22, height = 13, dpi = 300)


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
bloom
bloom$Region <- factor(bloom$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
bloom$Season <- factor(bloom$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
colors <- scales::hue_pal()(length(unique(phyto_abund$Region)))
palette <- setNames(colors, unique(phyto_abund$Region))
p <- ggplot(abund, aes(x = Region, y = Num_cell_l, fill = Region)) +
    geom_hline(yintercept = 1e6, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.9) + 
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
    #scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = palette) +
    facet_grid(Season ~ Basin, scales = "free_x", space = "free_x") +
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
        title = "Sample abundance per basin and season"
    ) 
p    
ggsave(file.path(HOME_, "abundance_per_basin_season.pdf"), p, width = 25, height = 20, dpi = 300)

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
        y = "Taxa richness",
        title = "Sample richness per basin and season"
    ) 
p
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


IndVal <- read_excel("./indval_per_basin.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
IndVal <- IndVal %>% dplyr::filter(Taxon != "Other phytoplankton")

order_species <- function(df, basins, threshold = 0.5) {
    characteristic_species <- df$Taxon[apply(df[, basins], 1, function(x) any(x >= threshold))]
    df_long <- reshape2::melt(df[, c("Taxon", basins)] %>% dplyr::filter(Taxon %in% characteristic_species), id.vars = "Taxon", variable.name = "Basin", value.name = "Value")
    df_long$Basin <- factor(df_long$Basin, levels = basins, ordered = TRUE)
    ordered_ids <- df_long %>%
    group_by(Taxon) %>%
    summarise(Max_Basin = Basin[which.max(Value)], Max_Value = max(Value)) %>%
    arrange(Max_Basin, desc(Max_Value)) %>% pull(Taxon)
    return(ordered_ids)
}



maps <- vector("list", 2)
taxa_list <- order_species(IndVal, ordered_basins, threshold = 0)
complete_indval <- IndVal %>% column_to_rownames(var = "Taxon") %>% .[taxa_list, ordered_basins] %>% as.matrix()
max_indval <- max(complete_indval)
col_fun = colorRamp2(c(0, max_indval / 2, max_indval), c("darkgreen", "white", "brown4"))
maps[[1]] <- grid.grabExpr(
    draw(
        Heatmap(complete_indval, 
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_centered = FALSE,
        show_column_names = TRUE,
        show_row_names = FALSE,
        column_names_rot = 45, 
        row_title = "Taxon",
        row_title_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 20),
        row_title_rot = 90,
        show_heatmap_legend = FALSE
        )
    )
)
taxa_list <- order_species(IndVal, ordered_basins, threshold = 0.5)
partial_indval <- IndVal %>% column_to_rownames(var = "Taxon") %>% .[taxa_list, ordered_basins] %>% as.matrix()
maps[[2]] <- grid.grabExpr(
    draw(
        Heatmap(partial_indval, 
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_centered = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_names_rot = 45, 
        row_names_gp = gpar(fontsize = 20, fontface = "italic"),
        column_names_gp = gpar(fontsize = 20),
        row_title_rot = 90,
        row_names_side = "left",
        heatmap_legend_param = list(
            title = "IndVal", 
            #at = seq(0, max(partial_indval), 0.2), 
            legend_width = unit(0.5, "inch"), 
            legend_height = unit(4, "inch"), 
            title_position = "topcenter",
            title_gp = gpar(fontsize = 25, fontface = "bold", lheight = 5),
            labels_gp = gpar(fontsize = 20), 
            grid_width = unit(0.5, "inch")
            ),
        width = unit(8, "inch"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", partial_indval[i, j]), x, y, gp = gpar(fontsize = 18))
        }
        ), 
        padding = unit(c(0.2, -2, 0.2, -7), "cm")
    )
)


combined_plot <- grid.arrange(grobs = maps, ncol=2, clip=FALSE, widths = c(1, 2)) 
text_grob_a <- grid::textGrob(label = "(a)",
          x=0.01, 
          y=0.98,
          gp=gpar(fontsize=20, col="black"))
text_grob_b <- grid::textGrob(label = "(b)",
            x=0.425, 
            y=0.98,
            gp=gpar(fontsize=20, col="black"))
grid.newpage()
grid.draw(combined_plot)
grid.draw(text_grob_a)
grid.draw(text_grob_b)
ggsave(file.path(HOME_, "IndVal_per_basin.pdf"), plot = grid.grab(), width = 22, height = 13, dpi = 600)



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


