library(dplyr)
library(tidyr)
library(tidytext)
library(ggplot2)

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

select_and_order <- function(data, title, n_otu = 10, n_samples = 2220) {
    data <- data[1:n_otu, ]
    data$Frequency <- data$Frequency <- round(data$Frequency / n_samples, 2)
    data$Taxon <- factor(data[[1]], levels = data[[1]])
    data$type <- title
    return(data %>% dplyr::select(Taxon, Frequency, type))
}
[order(data$Frequency, decreasing = FALSE)]
top_taxa <- read.csv("./Taxa_freq_95.csv")
top_species <- read.csv("./Species_freq_top.csv")
top_genera <- read.csv("./Genera_freq_top.csv")
top_classes <- read.csv("./Classes_freq_top.csv") %>% dplyr::filter(Class != 'nan')

top_taxa <- phyto_abund %>% group_by(Taxon) %>% summarise(Frequency = n_distinct(Date, id)) %>% dplyr::filter(Taxon != "Other phytoplankton") %>% arrange(desc(Frequency)) 
top_classes <- phyto_abund %>% group_by(Class) %>% summarise(Frequency = n_distinct(Date, id)) %>% arrange(desc(Frequency)) %>% dplyr::filter(Class != "nan")
top_genera <- phyto_abund %>% group_by(Genus) %>% summarise(Frequency = n_distinct(Date, id)) %>% arrange(desc(Frequency)) %>% dplyr::filter(Genus != "")
top_species <- phyto_abund %>% dplyr::filter(Det_level == "Species") %>% group_by(Taxon) %>% summarise(Frequency = n_distinct(Date, id)) %>% arrange(desc(Frequency)) %>% na.omit()


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
ggsave(file.path(HOME_, "relative_abundance_per_region.pdf"), p, width = 22, height = 13, dpi = 300)
