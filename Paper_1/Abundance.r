library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
#library(ggh4x)
#library(legendry)
library(rjson)
library(openxlsx)
library(colorBlindness)

HOME_ <- "./Paper_1"
IMAGE_FORMAT <- "svg"
source(file.path(HOME_, "utils.r"))


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

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))

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


abund_groups <- process_abund_groups(phyto_abund)


ordered_latitude <- abund %>% dplyr::select(id, Latitude) %>% distinct() %>% 
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Latitude) %>% round(2)
ordered_longitude <- abund %>% dplyr::select(id, Longitude) %>% distinct() %>%
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Longitude) %>% round(2)
p <- plot_variable_along_coast(
    abund %>% mutate(
        Abund = log10(Abund + 1)
    ), 
    var = "Abund", 
    group = "New_basin", 
    title = "Sample abundance across all stations", 
    ylab = "Abundance [cells/L] (log scale)", 
    ordered_latitude = ordered_latitude, 
    ordered_longitude = ordered_longitude
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




abund_groups %>% 
mutate(Month = format(as.Date(Date), "%m"), 
        Year_month = format(as.Date(Date), "%Y-%m")) %>%
pivot_wider(
    names_from = Group, 
    values_from = Abund, 
    values_fill = 0
) %>% mutate(
    DIA_DIN = DIA / DIN
) %>% ggplot() + 
geom_boxplot(aes(x = Year_month, y = log10(DIA_DIN), fill = Season))# + 
facet_wrap(~New_basin, ncol = 2)

p <- abund_groups %>% group_by(Season, New_basin, Group) %>%
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
    Abund = mean(Abund),
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













sheets <- getSheetNames(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"))
all_data <- sapply(sheets, function(sheet) {
    data <- openxlsx::read.xlsx(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), sheet = sheet, colNames = TRUE)
    colnames(data)[1] <- "Taxon"
    data <- data %>% dplyr::filter(Taxon != "Other phytoplankton")
    return(data)
}, simplify = FALSE)
names(all_data) <- sheets




quantile_abund_genera <- abund_only_genera %>% 
group_by(Date, id, Basin, Genus) %>% 
mutate(sample_abund = sum(Abund)) %>%
group_by(Season, Basin, Genus) %>%
summarise(
    Abund =quantile(Abund, 0.5),
    .groups = "drop"
)


dominant_species_path <- paste(HOME_, "Dominant_species", sep = "/")
dir.create(dominant_species_path, showWarnings = FALSE)

library(patchwork)
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
            ) %>% dplyr::filter(rel_abund > 0.03) %>% pull(Genus) %>% unique()
        barplot <-  quantile_abund_genera %>% dplyr::filter(Genus %in% dominant_genera & Basin == basin) %>% 
            mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
            group_by(Season, Basin) %>% 
            mutate(
                rel_abund = Abund / sum(Abund), 
                .groups = "drop"
            ) %>% 
            ggplot(aes(y = Season, fill = Genus, x = rel_abund)) +
            scale_y_discrete(limits = rev) +
            geom_bar(stat = "identity", linewidth = 0.5, color = "black") + 
            labs(x = "Relative abundance", y = "Season") +
            ggplot_theme + theme(legend.position = "bottom") + 
            guides(fill = guide_legend(ncol = 3)) + 
            scale_fill_manual(values = unname(colorBlindness::paletteMartin))
        
        selected_genera <- all_data[[basin]] %>% rowwise() %>%
        dplyr::filter(max(c_across(where(is.numeric))) > 0.3) %>%
        ungroup() %>% pull(Taxon)
        sites_genera <- abund_only_genera %>% dplyr::select(
            Date, id, Genus, Abund, Basin, Season
        ) %>% dplyr::filter(Genus %in% selected_genera & Basin == basin) %>% pivot_wider(
            names_from = Genus, 
            values_from = Abund, 
            values_fill = 0
        )
        Shannon_index <- sites_genera %>% dplyr::select(where(is.numeric)) %>% diversity()
        Pielou_even <- Shannon_index / log(specnumber(sites_genera %>% dplyr::select(where(is.numeric))))
        even <- cbind(
            sites_genera %>% dplyr::select(Date, id, Basin, Season),
            Pielou_even
        ) %>% mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE))
        evenness_plot <- even %>% 
            ggplot() + 
            geom_boxplot(aes(y = Season, x = Pielou_even)) +
            scale_y_discrete(limits = rev) +
            labs(x = "Pielou's evenness") +
            theme_minimal() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none",
                axis.text.x = element_text(angle = 0, hjust = 0, size = 17),
                axis.title.x = element_text(size = 22),
            )
        
        ggsave(
            file.path(dominant_species_path, basin),
            device = IMAGE_FORMAT, 
            barplot + evenness_plot + plot_layout(ncol = 2, widths = c(0.5, 0.5)), 
            width = 15, 
            height = 10, 
            dpi = 300
        )
    }
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



