library(iNEXT)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)
library(pals)
library(jsonlite)
library(gridExtra)
library(scales)

capitalize_first <- function(string) {
  paste0(toupper(substring(string, 1, 1)), substring(string, 2))
}

make_raref_plot_richness <- function(acc_curve, plot_title, taxon_level = "Taxa", plot_path = ".") {
  acculumated_richness <- acc_curve %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1)
  acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG",  "SAR"))

  custom_palette <- scales::hue_pal()(length(unique(acculumated_richness$Region)))
  names(custom_palette) <- unique(acculumated_richness$Region)
  line_types <- c("Rarefaction" = "solid", "Extrapolation" = "dotted")
  p1 <- ggplot() +
    geom_line(data = acc_curve %>% filter(Method != "Observed"), aes(x = t, y = qD, color = Region, linetype = Method), linewidth = 1) +
    geom_ribbon(data = acc_curve, aes(x = t, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = Region), alpha = 0.2) +
    scale_linetype_manual(values = line_types) +
    scale_color_manual(values = custom_palette) +
    scale_fill_manual(values = custom_palette) +
    geom_point(data = acc_curve %>% filter(Method == "Observed"), aes(x = t, y = qD, color = Region), size = 3) +
    theme_bw(base_size = 18) + 
    geom_vline(xintercept = 50, linetype = "dashed", color = "black") + 
    labs(x = "Number of samples", y = paste(capitalize_first(taxon_level), "richness", sep = " ")) +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
    
  p2 <- ggplot(acculumated_richness) + 
    geom_point(aes(x = Region, y = qD, color = Region), size = 7) +
    geom_errorbar(aes(x = Region, ymin = qD.LCL, ymax = qD.UCL, color = Region), width = 0.5)  +
    theme_bw(base_size = 18) +
    scale_color_manual(values = custom_palette) +
    labs(x = "Region", y = paste("Estimated", tolower(taxon_level), "richness", sep = " ")) +
    theme(legend.position = "none") +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
  ggsave(paste(plot_path, plot_title, sep = "/"), grid.arrange(p1, p2, ncol = 2, widths = c(1, 1), heights = c(1)), width = 25, height = 10)
}

# Load data
HOME_ <- "."

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))

plot_path <- paste(HOME_, "ISPRA_20152017_Analysis/Plots/Rich_levels/Acc_curve", sep = "/")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
for (region in names(data)) {
        data[[region]][, c(1:3)] <- data[[region]][, c(1:3)] %>% fill(Season, id)
        incidence_freq_list <- list()
        incidence_freq_dataframe <- merge(
        data[[region]] %>% group_by(Season) %>% summarise(n_samples = n_distinct(id, Date)), 
        data[[region]] %>% select(-c(id, Date)) %>% group_by(Season) %>% summarise_all(~sum(. != 0)), 
        by = "Season"
        )
        for (season in unique(incidence_freq_dataframe$Season)) {
            incidence_freq_list[[season]] <- incidence_freq_dataframe[incidence_freq_dataframe$Season == season, -1] %>% .[. != 0]
        }
        endpoint <- max(incidence_freq_dataframe$n_samples) + 10
        res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
        ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18)
        ggsave(paste(plot_path, "/acc_curve_", region, ".png", sep = ""), width = 10, height = 10)
}

make_raref_plot_richness(acc_curve.size_based, "acc_curve_all_regions_only_species.pdf", "Species")



#SAC per region
species_observation <- phyto_abund %>%
filter(Det_level == "Species") %>%
group_by(Region, Date, id, Taxon) %>% summarise(Abund = sum(Num_cell_l)) %>%
pivot_wider( names_from = Taxon, values_from = Abund, values_fill = 0) %>% ungroup()

genus_observation <- phyto_abund %>% 
filter(Det_level == "Genus" | Det_level == "Species") %>% 
group_by(Region, Date, id, Genus) %>% summarise(Abund = sum(Num_cell_l)) %>% 
pivot_wider( names_from = Genus, values_from = Abund, values_fill = 0) %>% ungroup()

common_species <- phyto_abund %>% filter(Det_level == "Species") %>% group_by(Taxon) %>% 
summarise(Frequency = n_distinct(Date, id) / 2220) %>% arrange(desc(Frequency)) %>% filter(Frequency > 0.1) %>% pull(Taxon) 
common_species_observation <- phyto_abund %>% filter(Taxon %in% common_species) %>% 
group_by(Region, Date, id, Taxon) %>% summarise(Abund = sum(Num_cell_l), .groups = "drop") %>%
pivot_wider( names_from = Taxon, values_from = Abund, values_fill = 0) 

common_genera <- phyto_abund %>% filter(Det_level == "Genus" | Det_level == "Species") %>% group_by(Genus) %>%
summarise(Frequency = n_distinct(Date, id) / 2220) %>% arrange(desc(Frequency)) %>% filter(Frequency > 0.1) %>% pull(Genus)
common_genera_observation <- phyto_abund %>% filter(Genus %in% common_genera) %>%
group_by(Region, Date, id, Genus) %>% summarise(Abund = sum(Num_cell_l), .groups = "drop") %>%
pivot_wider( names_from = Genus, values_from = Abund, values_fill = 0)

incidence_freq_list <- list()
for (region in unique(species_observation$Region)) {
  region_df <- species_observation %>% filter(Region == region)
  incidence_freq_dataframe <- merge(
      region_df %>% summarise(n_samples = n_distinct(id, Date)),
      region_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[region]] <- incidence_freq_dataframe
}
endpoint <- species_observation %>% group_by(Region) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_all_regions_only_species.csv", sep = "/"))


incidence_freq_list <- list()
for (region in unique(genus_observation$Region)) {
  regional_df <- genus_observation %>% filter(Region == region)
  incidence_freq_dataframe <- merge(
      regional_df %>% summarise(n_samples = n_distinct(id, Date)),
      regional_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[region]] <- incidence_freq_dataframe
}
endpoint <- genus_observation %>% group_by(Region) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_all_regions_only_genera.csv", sep = "/"))


incidence_freq_list <- list()
for (region in unique(common_species_observation$Region)) {
  region_df <- common_species_observation %>% filter(Region == region)
  incidence_freq_dataframe <- merge(
      region_df %>% summarise(n_samples = n_distinct(id, Date)),
      region_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[region]] <- incidence_freq_dataframe
}
endpoint <- common_species_observation %>% group_by(Region) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_all_regions_only_common_species.csv", sep = "/"))

incidence_freq_list <- list()
for (region in unique(common_genera_observation$Region)) {
  region_df <- common_genera_observation %>% filter(Region == region)
  incidence_freq_dataframe <- merge(
      region_df %>% summarise(n_samples = n_distinct(id, Date)),
      region_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[region]] <- incidence_freq_dataframe
}
endpoint <- common_genera_observation %>% group_by(Region) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_all_regions_only_common_genera.csv", sep = "/"))


acc_curve_genera <- read.csv(paste(HOME_, "acc_curve_all_regions_only_genera.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve_species <- read.csv(paste(HOME_, "acc_curve_all_regions_only_species.csv", sep = "/")) %>% rename(Region = Assemblage) 
acculumated_richness <- rbind(
  acc_curve_genera %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1) %>% 
  mutate(Level = "Genera") , 
  acc_curve_species %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1) %>% 
  mutate(Level = "Species") 
)
acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"))


custom_palette <- rep("black", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness %>% arrange(desc(Level)), aes(x = Region)) + 
    geom_bar(stat = "identity", position = position_dodge(), aes(y = qD, fill = Level), colour = "black", size = 1, alpha = 1) +
    geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL, group = Level), width = 0.5, position = position_dodge(width = 1))  +
    theme_bw(base_size = 18) +
    labs(x = "Region", y = "Estimated richness") +
    ggtitle("Estimated richness of species and genera") +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(t = -10)
    )
plot_title <- "acc_curve_all_regions_genera_species_bar.pdf"
ggsave(paste(HOME_, plot_title, sep = "/"), p, width = 15, height = 10)


acc_curve_genera <- read.csv(paste(HOME_, "acc_curve_all_regions_only_common_genera.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve_species <- read.csv(paste(HOME_, "acc_curve_all_regions_only_common_species.csv", sep = "/")) %>% rename(Region = Assemblage) 
acculumated_richness <- rbind(
  acc_curve_genera %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1) %>% 
  mutate(Level = "Genera") , 
  acc_curve_species %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1) %>% 
  mutate(Level = "Species") 
)
acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"))

custom_palette <- rep("black", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness %>% arrange(desc(Level)), aes(x = Region)) + 
    geom_bar(stat = "identity", position = position_dodge(), aes(y = qD, fill = Level), colour = "black", linewidth = 1, alpha = 1) +
    geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL, group = Level), width = 0.5, position = position_dodge(width = 1))  +
    theme_bw(base_size = 18) +
    labs(x = "Region", y = "Estimated richness") +
    ggtitle("Estimated richness of species and genera") +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(t = -10)
    )
p
plot_title <- "acc_curve_all_regions_genera_species_bar.pdf"
ggsave(paste(HOME_, plot_title, sep = "/"), p, width = 15, height = 10)

phyto_abund %>% dplyr::filter(Det_level == "Unknown") %>% group_by(Region) %>% 
ggplot() + 
geom_boxplot(aes(x = Region, y = log10(Num_cell_l)))

phyto_abund %>% group_by(Region, Date, id) %>% summarise(
  Abund = sum(Num_cell_l[Det_level == "Unknown"]),
  n_species = n_distinct(Taxon[Det_level == "Species"]),
  n_genera = n_distinct(Genus[Det_level == "Genus"])
) %>% ggplot() + 
geom_point(aes(x = log10(Abund + 1), y = n_genera, color = Region))
##PER BASIN##
phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))
genus_observation <- phyto_abund %>% 
filter(Det_level == "Genus" | Det_level == "Species") %>% 
group_by(Basin, Date, id, Genus) %>% summarise(Abund = sum(Num_cell_l)) %>% 
pivot_wider( names_from = Genus, values_from = Abund, values_fill = 0) %>% ungroup()
species_observation <- phyto_abund %>%
filter(Det_level == "Species") %>%
group_by(Basin, Date, id, Taxon) %>% summarise(Abund = sum(Num_cell_l)) %>%
pivot_wider( names_from = Taxon, values_from = Abund, values_fill = 0) %>% ungroup()

incidence_freq_list <- list()
for (basin in unique(genus_observation$Basin)) {
  basin_df <- genus_observation %>% filter(Basin == basin)
  incidence_freq_dataframe <- merge(
      basin_df %>% summarise(n_samples = n_distinct(id, Date)),
      basin_df %>% select(-c(id, Date, Basin)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[basin]] <- incidence_freq_dataframe
}
endpoint <- genus_observation %>% group_by(Basin) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_basins_only_genera.csv", sep = "/"))

incidence_freq_list <- list()
for (basin in unique(species_observation$Basin)) {
  basin_df <- species_observation %>% filter(Basin == basin)
  incidence_freq_dataframe <- merge(
      basin_df %>% summarise(n_samples = n_distinct(id, Date)),
      basin_df %>% select(-c(id, Date, Basin)) %>% summarise_all(~sum(. != 0))
      )
  incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
  incidence_freq_list[[basin]] <- incidence_freq_dataframe
}
endpoint <- species_observation %>% group_by(Basin) %>% summarise(n_samples = n_distinct(id, Date)) %>% pull(n_samples) %>% max() + 10
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_basins_only_species.csv", sep = "/"))


acc_curve_genera <- read.csv(paste(HOME_, "acc_curve_basins_only_genera.csv", sep = "/")) %>% rename(Basin = Assemblage) 
acc_curve_species <- read.csv(paste(HOME_, "acc_curve_basins_only_species.csv", sep = "/")) %>% rename(Basin = Assemblage) 
acculumated_richness <- rbind(
  acc_curve_genera %>%
  group_by(Basin) %>%
  filter(abs(t - 300) == min(abs(t - 300))) %>% slice(1) %>% 
  mutate(Level = "Genera") , 
  acc_curve_species %>%
  group_by(Basin) %>%
  filter(abs(t - 300) == min(abs(t - 300))) %>% slice(1) %>% 
  mutate(Level = "Species") 
)
acculumated_richness$Basin <- factor(acculumated_richness$Basin, levels = c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed"))

custom_palette <- rep("black", length(unique(acculumated_richness$Basin)))
names(custom_palette) <- unique(acculumated_richness$Basin)
p <- ggplot(acculumated_richness %>% arrange(desc(Level)), aes(x = Basin)) + 
    geom_bar(stat = "identity", position = position_dodge(), aes(y = qD, fill = Level), colour = "black", linewidth = 1, alpha = 1) +
    geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL, group = Level), width = 0.5, position = position_dodge(width = 1))  +
    theme_bw(base_size = 18) +
    #scale_color_manual(values = custom_palette) +
    labs(x = "Basin", y = "Estimated richness") +
    ggtitle("Estimated richness of species and genera") +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 15),
      #legend.title = element_text(size = 18, face = "bold"), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(t = -10)
    )
plot_title <- "acc_curve_all_basins_genera_species_bar.pdf"
ggsave(paste(HOME_, plot_title, sep = "/"), p, width = 15, height = 10)


##TESTING DEPENDENCY ON COAST LENGTH##

##values from https://www.isprambiente.gov.it/files/pubblicazioni/statoambiente/tematiche-2012/Cap.5_Mare_ambiente_costiero.pdf
coast_length <- data.frame(
  Region = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"),
  Coast_length = c(116, 216, 174, 176, 129, 37, 957, 66, 734, 1603, 502, 380, 646, 378, 2160)
)

acculumated_richness <- merge(acculumated_richness, coast_length, by = "Region")


p <- ggplot(acculumated_richness %>% arrange(desc(Level)), aes(x = Region)) + 
    geom_bar(stat = "identity", position = position_dodge(), aes(y = qD / Coast_length, fill = Level), colour = "black", size = 1, alpha = 1) +
    geom_errorbar(aes(ymin = qD.LCL / Coast_length, ymax = qD.UCL / Coast_length, group = Level), width = 0.5, position = position_dodge(width = 1))  +
    theme_bw(base_size = 18) +
    #scale_color_manual(values = custom_palette) +
    labs(x = "Region", y = "Estimated richness") +
    ggtitle("Estimated richness of species and genera") +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(t = -10)
    )