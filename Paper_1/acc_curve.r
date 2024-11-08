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

sheet_names <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"))

data <- list()  # Create an empty list to store the data from each sheet
data <- lapply(sheet_names, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"), sheet = sheet))
names(data) <-  sheet_names

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


acc_curve.size_based <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/acc_curve_all_regions_all_taxa.csv", sep = "/")) %>% rename(Region = Assemblage) 
acc_curve.size_based <- read.csv(paste(HOME_, "acc_curve_all_regions_only_species.csv", sep = "/")) %>% rename(Region = Assemblage) 

capitalize_first <- function(string) {
  paste0(toupper(substring(string, 1, 1)), substring(string, 2))
}

make_raref_plot_richness <- function(acc_curve, plot_title, taxon_level = "Taxa", plot_path = ".") {
  acculumated_richness <- acc_curve %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1)
  acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "CAM","LAZ", "TOS", "LIG", "SIC", "SAR"))

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


#original plots
N <- length(incidence_freq_list)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
#scale_fill_manual(values = unname(alphabet(N))) + 
#scale_colour_manual(values = unname(alphabet(N))) + 
 guides(shape=FALSE)
ggsave(paste(plot_path, "/acc_curve_all_regions.png", sep = ""), width = 20, height = 10)

res
res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 1, conf = 0.95, endpoint = 200)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
scale_fill_manual(values = unname(alphabet(N))) + 
scale_colour_manual(values = unname(alphabet(N))) + 
guides(shape=FALSE)

res$iNextEst$size_based %>% filter(between(SC, 0.90,0.95))