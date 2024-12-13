library(iNEXT)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)
library(pals)
library(jsonlite)
library(gridExtra)
library(scales)

# Load data
HOME_ <- "."

data <- read.csv(paste(HOME_, "sites_taxa.csv", sep = "/"))

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


data <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/site_taxa_matrix_w_det_level.csv", sep = "/"))
data[1, c(1,2,3)] <- data[2, c(1,2,3)]
data <- data[-2, ]
is_species <- grepl("Sp", colnames(data))
colnames(data) <- data[1, ]
species_list <- colnames(data)[is_species]
data <- data[-1, ] %>%
  mutate_at(vars(-c(1, 2, 3)), as.numeric)
regions <- unique(data$Region)

incidence_freq_list <- list()
for (region in regions) {
  regional_df <- data %>% filter(Region == region)
        incidence_freq_dataframe <- merge(
            regional_df %>% summarise(n_samples = n_distinct(id, Date)),
            regional_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
            )
        incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
        incidence_freq_list[[region]] <- incidence_freq_dataframe
}
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = 200)
res$iNextEst
write.csv(res$iNextEst$size_based, paste(HOME_, "ISPRA_20152017_Analysis/acc_curve_all_regions_all_taxa.csv", sep = "/"))



incidence_freq_list <- list()
only_species <- data %>% select(1:3, all_of(species_list))
for (region in regions) {
  regional_df <- only_species %>% filter(Region == region)
        incidence_freq_dataframe <- merge(
            regional_df %>% summarise(n_samples = n_distinct(id, Date)),
            regional_df %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
            )
        incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
        incidence_freq_list[[region]] <- incidence_freq_dataframe
}
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = 200)
write.csv(res$iNextEst$size_based, paste(HOME_, "ISPRA_20152017_Analysis/acc_curve_all_regions_only_species.csv", sep = "/"))


# SACs for genera
phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))
genus_observation <- phyto_abund %>% 
filter(Det_level == "Genus" | Det_level == "Species") %>% 
group_by(Region, Date, id, Genus) %>% summarise(Abund = sum(Num_cell_l)) %>% 
pivot_wider( names_from = Genus, values_from = Abund, values_fill = 0) %>% ungroup()
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
incidence_freq_list
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)
res$iNextEst
write.csv(res$iNextEst$size_based, paste(HOME_, "acc_curve_all_regions_only_genera.csv", sep = "/"))



acc_curve.size_based <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/acc_curve_all_regions_all_taxa.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve.size_based <- read.csv(paste(HOME_, "acc_curve_all_regions_only_species.csv", sep = "/")) %>% rename(Region = Assemblage) 

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


make_raref_plot_richness(acc_curve.size_based, "acc_curve_all_regions_only_species.pdf", "Species")

taxon_level <- "Species"
plot_path <- "."
acc_curve <- acc_curve.size_based
acculumated_richness <- acc_curve %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1)
acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"))
custom_palette <- rep("dodgerblue3", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness) + 
    geom_bar(stat = "identity", aes(x = Region, y = qD, fill = Region), colour = "black", linewidth = 1) +
    geom_errorbar(aes(x = Region, ymin = qD.LCL, ymax = qD.UCL), linewidth = 0.5, width = 0.5)  +
    theme_bw(base_size = 18) +
    scale_fill_manual(values = custom_palette) +
    labs(x = "Region", y = paste("Estimated", tolower(taxon_level), "richness", sep = " ")) +
    theme(legend.position = "none") +
    theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) 
p
plot_title <- "acc_curve_all_regions_only_species_bar.pdf"
ggsave(paste(plot_path, plot_title, sep = "/"), p, width = 15, height = 10)


acc_curve.size_based <- read.csv(paste(HOME_, "acc_curve_all_regions_only_genera.csv", sep = "/")) %>% rename(Region = Assemblage) 
taxon_level <- "Genera"
plot_path <- "."
acc_curve <- acc_curve.size_based
acculumated_richness <- acc_curve %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1)
acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"))
custom_palette <- rep("dodgerblue3", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness) + 
    geom_bar(stat = "identity", aes(x = Region, y = qD, fill = Region), colour = "black", linewidth = 1) +
    geom_errorbar(aes(x = Region, ymin = qD.LCL, ymax = qD.UCL), linewidth = 0.5, width = 0.5)  +
    theme_bw(base_size = 18) +
    scale_fill_manual(values = custom_palette) +
    labs(x = "Region", y = paste("Estimated", tolower(taxon_level), "richness", sep = " ")) +
    theme(legend.position = "none") +
    theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) 
p
plot_title <- "acc_curve_all_regions_only_genera_bar.pdf"
ggsave(paste(plot_path, plot_title, sep = "/"), p, width = 15, height = 10)


acc_curve.size_based <- read.csv(paste(HOME_, "acc_curve_all_regions_only_genera.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve_genera <- acc_curve.size_based
acc_curve.size_based <- read.csv(paste(HOME_, "acc_curve_all_regions_only_species.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve_species <- acc_curve.size_based
taxon_level <- "Genera"
plot_path <- "."

acc_curve_genera %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1) 
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


p <- ggplot(acculumated_richness) + 
    geom_point(aes(x = Region, y = qD, color = Level), size = 7) +
    geom_errorbar(aes(x = Region, ymin = qD.LCL, ymax = qD.UCL, color = Taxon), width = 0.5)  +
    theme_bw(base_size = 18) +
    #scale_color_manual(values = custom_palette) +
    labs(x = "Region", y = paste("Estimated", tolower(taxon_level), "richness", sep = " ")) +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
p
custom_palette <- rep("black", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness %>% arrange(desc(Level)), aes(x = Region)) + 
    geom_bar(stat = "identity", position = position_dodge(), aes(y = qD, fill = Level), colour = "black", size = 1, alpha = 1) +
    geom_errorbar(aes(ymin = qD.LCL, ymax = qD.UCL, group = Level), width = 0.5, position = position_dodge(width = 1))  +
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
      #legend.title = element_text(size = 18, face = "bold"), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(t = -10)
    )
p
plot_title <- "acc_curve_all_regions_genera_species_bar.pdf"
ggsave(paste(plot_path, plot_title, sep = "/"), p, width = 15, height = 10)




acculumated_richness
#original plots
N <- length(incidence_freq_list)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
#scale_fill_manual(values = unname(alphabet(N))) + 
#scale_colour_manual(values = unname(alphabet(N))) + 
 guides(shape=FALSE)
ggsave(paste(plot_path, "/acc_curve_all_regions.png", sep = ""), width = 20, height = 10)


acculumated_richness 

merge(genus_observation %>% group_by(Region) %>%
summarise(n_samples = n_distinct(id, Date), n_stations = n_distinct(id)), 
acculumated_richness  %>% select(Region, qD, Level)
) %>% ggplot() + 
geom_point(aes(x = n_samples, y = qD, shape = Level, color = Level), size = 5)

res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 1, conf = 0.95, endpoint = 200)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
scale_fill_manual(values = unname(alphabet(N))) + 
scale_colour_manual(values = unname(alphabet(N))) + 
guides(shape=FALSE)

res$iNextEst$size_based %>% filter(between(SC, 0.90,0.95))