library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(readxl)
library(reshape2)
library(grid)
library(tidyverse)
library(paletteer)

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
HOME_ <- "."
phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))
IndVal <- read_excel("./indval_only_genera.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
#IndVal_taxon <- c(IndVal %>% filter(apply(select(., -Taxon), 1, max) > 0.5) %>% pull(Taxon), "Unknown")
ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")
IndVal_taxon <- order_species(IndVal, ordered_basins, threshold = 0.5)
IndVal_taxon <- IndVal_taxon[! IndVal_taxon %in% c("Bacillariophyceae", "Dinoflagellata")]



### BIODIVERSITY ###
data_plot <- phyto_abund %>% dplyr::select(c(Basin, Date, id, Season, Genus, Num_cell_l)) %>% group_by(Basin, Season) %>% mutate(n_sample_basin = n_distinct(id, Date)) %>%
dplyr::filter(Genus %in% IndVal_taxon) %>% group_by(Basin, Season, Genus) %>% summarise(Abund = sum(Num_cell_l), n_sample_basin = first(n_sample_basin), .groups = "drop") %>% 
mutate(Abund = Abund / n_sample_basin) %>% complete(Basin, Season, Genus, fill = list(Abund = 0))
data_plot <- phyto_abund %>% dplyr::select(c(Basin, Date, id, Season, Genus, Num_cell_l)) %>% group_by(Basin, Season) %>% mutate(n_sample_basin = n_distinct(id, Date)) %>%
dplyr::filter(Genus %in% IndVal_taxon) %>% group_by(Basin, Season, Genus) %>% summarise(Abund = median(Num_cell_l), n_sample_basin = first(n_sample_basin), .groups = "drop") %>% 
complete(Basin, Season, Genus, fill = list(Abund = 0))
data_plot$Genus <- factor(data_plot$Genus, level = IndVal_taxon)
data_plot$Season <- factor(data_plot$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
data_plot$Basin <- factor(data_plot$Basin, levels = ordered_basins)
colors <- paletteer::paletteer_d("ggsci::category20_d3")
custom_palette <- setNames(
    as.character(colors[1:length(unique(data_plot$Genus))]), 
    unique(data_plot$Genus)
    )



data_plot %>% mutate(Abund = log10(Abund +1)) %>% 
ggplot(aes(x = Season, y = Basin, fill = Abund)) +
facet_wrap(~Genus) +
geom_tile() + 
geom_text(aes(label = round(Abund, 2), colour = ifelse(Abund > 3, "black", "white")), size = 6) +
theme_bw() +
scale_fill_continuous(type = "viridis") + 
scale_colour_manual(values=c("white"="white", "black"="black")) +
scale_y_discrete(limits = rev) + 
theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    ) + 
    guides(
        color = "none", 
        fill = guide_colourbar(
            title = "IndVal", 
            title.position = "left", 
            title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
            label.theme = element_text(size = 20),
            barwidth = unit(20, "lines"),
            ticks.linewidth = 1,
            frame.linewidth = 1,
            ticks.colour = "black",
            frame.colour  ='black'
            )
    )

data_plot %>% mutate(Abund = log10(Abund +1)) %>% 
ggplot(aes(x = Season, y = Genus, fill = Abund)) +
facet_wrap(~Basin, nrow = 1) +
geom_tile() + 
geom_text(aes(label = round(Abund, 2), colour = ifelse(Abund > 3, "black", "white")), size = 6) +
theme_bw() +
scale_fill_continuous(type = "viridis") + 
scale_y_discrete(limits = rev) + 
scale_colour_manual(values=c("white"="white", "black"="black")) +
theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    )+ 
    guides(
        color = "none", 
        fill = guide_colourbar(
            title = "IndVal", 
            title.position = "left", 
            title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
            label.theme = element_text(size = 20),
            barwidth = unit(20, "lines"),
            ticks.linewidth = 1,
            frame.linewidth = 1,
            ticks.colour = "black",
            frame.colour  ='black'
            )
    )


## TESTS ##
mean_abund_per_season_basin <- phyto_abund %>% mutate(
    Det_level = case_when(
        Det_level == "Higher cat." ~ Taxon,
        Class == "nan" ~ Det_level,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
    )  %>% dplyr::filter(Det_level != "Unknown")  %>% group_by(Basin, Season, Date, id, Det_level) %>%
    summarise(Abund = sum(Num_cell_l), .groups = "drop") %>% group_by(Basin, Season, Det_level) %>%
    summarise(mean_abund = mean(Abund),.groups = "drop") %>% group_by(Basin, Season) %>% 
    mutate(tot_abund = sum(mean_abund), rel_abund = mean_abund / tot_abund)


colors <- paletteer::paletteer_d("ggsci::category20_d3")
custom_palette <- setNames(
    as.character(colors[1:length(unique(IndVal_taxon))]), 
    unique(IndVal_taxon)
    )
custom_palette["Else"] <- "darkgrey"
abund_most_important_taxa <- phyto_abund %>% dplyr::select(c(Basin, Date, id, Season, Genus, Num_cell_l)) %>% 
mutate(Genus = ifelse(Genus %in% IndVal_taxon, Genus, "Else")) %>% 
group_by(Basin, Season, Genus) %>% summarise(Abund = sum(Num_cell_l), .groups = "drop") %>% 
group_by(Basin, Season) %>% mutate(rel_abund = Abund / sum(Abund)) 
abund_most_important_taxa$Genus <- factor(abund_most_important_taxa$Genus, level = c(IndVal_taxon, "Else"))
## NO ELSE ###
abund_most_important_taxa <- phyto_abund %>% dplyr::select(c(Basin, Date, id, Season, Genus, Num_cell_l)) %>% 
dplyr::filter(Genus %in% IndVal_taxon) %>% 
group_by(Basin, Season, Genus) %>% summarise(Abund = sum(Num_cell_l), .groups = "drop") %>% 
group_by(Basin, Season) %>% mutate(rel_abund = Abund / sum(Abund)) 
abund_most_important_taxa$Genus <- factor(abund_most_important_taxa$Genus, level = IndVal_taxon)
## NO ELSE ###
abund_most_important_taxa$Season <- factor(abund_most_important_taxa$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
abund_most_important_taxa$Basin <- factor(abund_most_important_taxa$Basin, levels = ordered_basins)

abund_most_important_taxa %>% 
ggplot() +
geom_col(aes(x = Season, y = rel_abund, group = Genus, fill = Genus), color = "black", linewidth = 0.5) +
facet_wrap(~Basin) +
theme_bw() +
labs(y = "Relative abundace") +
scale_fill_manual(values = custom_palette) + 
theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25, face = "bold"),
        #legend.position = "none", 
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    ) + 
guides(fill = guide_legend(ncol = 1))


freq_most_important_taxa <- phyto_abund %>% dplyr::select(c(Basin, Date, id, Season, Genus, Num_cell_l)) %>% group_by(Basin, Season) %>% mutate(n_sample_basin = n_distinct(id, Date)) %>%
dplyr::filter(Genus %in% IndVal_taxon) %>% distinct(Date, id, Genus, .keep_all = TRUE) %>% group_by(Basin, Season, Genus) %>% 
summarise(Freq = n(), n_sample_basin = first(n_sample_basin), .groups = "drop") %>%
mutate(Freq = Freq / n_sample_basin) %>% complete(Basin, Season, Genus, fill = list(Freq = 0)) 
freq_most_important_taxa$Genus <- factor(freq_most_important_taxa$Genus, level = c(IndVal_taxon, "Else"))
freq_most_important_taxa$Basin <- factor(freq_most_important_taxa$Basin, levels = ordered_basins)
freq_most_important_taxa$Season <- factor(freq_most_important_taxa$Season, levels = c("Winter", "Spring", "Summer", "Autumn")) 

freq_most_important_taxa %>% ggplot(aes(x = Season, y = Genus, fill = Freq)) +
geom_tile() + 
geom_text(aes(label = round(Freq, 2), colour = ifelse(Freq > 0.5, "black", "white")), size = 5.5) +
facet_wrap(~Basin) +
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22, face = "italic"),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    legend.position = "bottom", 
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
scale_colour_manual(values=c("white"="white", "black"="black")) + 
guides(colour = "none") + guides(
    fill = guide_colourbar(
        title = "Frequency", 
        title.position = "left", 
        title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
        label.theme = element_text(size = 20),
        barwidth = unit(20, "lines"),
        ticks.linewidth = 1,
        frame.linewidth = 1,
        ticks.colour = "black",
        frame.colour  ='black'
        )
    )


average_genera_cont <- mean_abund_per_season_basin %>% dplyr::filter(Det_level %in% IndVal_taxon) %>% group_by(Basin, Season, Det_level) %>% 
summarise(mean_abund = sum(mean_abund), .groups = "drop") %>% group_by(Basin, Season) %>%  
mutate(tot_abund = sum(mean_abund), rel_abund = mean_abund / tot_abund) %>% rename(Taxa = Det_level)
average_genera_cont$Basin <- factor(average_genera_cont$Basin, levels = ordered_basins)
average_genera_cont$Season <- factor(average_genera_cont$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))


### BARPLOT FOR BASIN ###   (TO CONTROL)
data_plot %>%
ggplot(aes(x = Season, y = Abund +1, group = Genus, fill = Genus, color = Genus)) +
facet_wrap(~Basin, ncol = 2) +
geom_bar(stat="identity", position=position_dodge(), color = "black", linewidth = 0.5) + 
labs(y = "Propostion of abundace") + 
theme_bw() +
scale_y_log10() +
scale_fill_manual(values = custom_palette) + 
theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    ) + 
guides(fill = guide_legend(ncol = 7))

#guides(fill = guide_legend(ncol = 7))


### BOXPLOT ###
abund_most_important_taxa <-  phyto_abund %>% mutate(
    Taxa = case_when(
    Det_level == "Higher cat." ~ Taxon,
    Class == "nan" ~ Det_level,
    Genus == "" ~ Class,
    TRUE ~ Genus
    )
    ) %>% dplyr::filter(Taxa %in% IndVal_taxon)  %>% group_by(Basin, Season, Date, id, Taxa) %>%
    summarise(Abund = sum(Num_cell_l), .groups = "drop")
abund_most_important_taxa$Season <- factor(abund_most_important_taxa$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
abund_most_important_taxa$Basin <- factor(abund_most_important_taxa$Basin, levels = ordered_basins)
abund_most_important_taxa$Taxa <- factor(abund_most_important_taxa$Taxa, level = IndVal_taxon)

abund_most_important_taxa %>% 
ggplot(aes(x = Taxa, y = Abund, group = Taxa, fill = Taxa, color = Taxa)) +
facet_grid(Basin ~ Season, scale = "free_x") +
geom_boxplot(color = "black", linewidth = 0.2) + 
labs(y = "Propostion of abundace") + 
theme_bw() +
scale_fill_paletteer_d("ggsci::category20_d3") +
scale_y_log10() + 
theme(
        #axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        #axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        #strip.text.x = element_text(size = 16), 
        #strip.text.y = element_text(size = 16), 
        #panel.spacing = unit(1, "lines")
    ) + 
guides(fill = guide_legend(ncol = 8))


### BARPLOT DIVIDED BY SPECIES ###
p <- average_genera_cont %>% ungroup()  %>% complete(Basin, Taxa, Season) %>% replace(., is.na(.), 1) %>% 
ggplot(aes(x = mean_abund, y = Taxa, group = Season, fill = Season, color = Season)) +
facet_wrap(~Basin, scale = "free_x") +
geom_col(position = position_dodge2(width = 0.5, preserve = "single"), color = "black", linewidth = 1) + 
labs(y = "Propostion of abundace") + 
theme_bw() +
scale_x_log10() + 
theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
    ) + 
guides(fill = guide_legend(ncol = 8))

ggsave(file.path(HOME_, "average_genera_cont.png"), p, width = 13, height = 22, dpi = 300)


IndVal_A <- read_excel("./indval_only_genera.xlsx", sheet = "Indval_A", col_names = TRUE) 
colnames(IndVal_A)[1] <- "Taxon"
IndVal_B <- read_excel("./indval_only_genera.xlsx", sheet = "Indval_B", col_names = TRUE) 
colnames(IndVal_B)[1] <- "Taxon"

IndVal <- read_excel("./indval_only_genera.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
#IndVal_taxon <- c(IndVal %>% filter(apply(select(., -Taxon), 1, max) > 0.5) %>% pull(Taxon), "Unknown")
IndVal_taxon <- order_species(IndVal, ordered_basins, threshold = 0.5)
IndVal_taxon <- IndVal_taxon[! IndVal_taxon %in% c("Bacillariophyceae", "Dinoflagellata")]


indval_components <- rbind(
    IndVal_A %>% dplyr::filter(Taxon %in% IndVal_taxon) %>% pivot_longer(cols = -Taxon, names_to = "Basin", values_to = "IndVal") %>% mutate(Component = "A"),
    IndVal_B %>% dplyr::filter(Taxon %in% IndVal_taxon) %>% pivot_longer(cols = -Taxon, names_to = "Basin", values_to = "IndVal") %>% mutate(Component = "B")
) 
indval_components$Basin <- factor(indval_components$Basin, levels = ordered_basins)
indval_components$Taxon <- factor(indval_components$Taxon, level = IndVal_taxon)

indval_components %>% ggplot(aes(x = Basin, y = Taxon)) +
geom_tile(aes(fill = IndVal)) +
geom_text(aes(label = round(IndVal, 2), colour = ifelse(IndVal > 0.5, "black", "white")), size = 5.5) +
facet_wrap(~Component) +
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22, face = "italic"),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    legend.position = "bottom", 
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
scale_colour_manual(values=c("white"="white", "black"="black")) + 
guides(colour = "none") + guides(
    fill = guide_colourbar(
        title = "Frequency", 
        title.position = "left", 
        title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
        label.theme = element_text(size = 20),
        barwidth = unit(20, "lines"),
        ticks.linewidth = 1,
        frame.linewidth = 1,
        ticks.colour = "black",
        frame.colour  ='black'
        )
    )
