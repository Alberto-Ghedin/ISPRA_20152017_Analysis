library(iNEXT)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)
library(pals)
library(jsonlite)

# Load data
HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

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
is_species
colnames(data) <- data[1, ]
data <- data %>%
  mutate_at(vars(-c(1, 2, 3)), as.numeric)
regions <- unique(data$Region)

incidence_freq_list <- list()
for (region in regions) {
        incidence_freq_dataframe <- merge(
            data %>% filter(Region == region) %>% summarise(n_samples = n_distinct(id, Date)),
            data %>% filter(Region == region) %>% select(-c(id, Date, Region)) %>% summarise_all(~sum(. != 0))
            )
        incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
        incidence_freq_list[[region]] <- incidence_freq_dataframe
}


res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = 200)
write.csv(res$iNextEst$size_based, paste(HOME_, "ISPRA_20152017_Analysis/acc_curve_all_regions_all_taxa.csv", sep = "/"))

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






