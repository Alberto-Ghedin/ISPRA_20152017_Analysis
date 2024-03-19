library(iNEXT)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(pals)

# Load data
HOME_ <- "/mnt/d" 
sheet_names <- excel_sheets("/mnt/d/PHD/ISPRA_20152017_Analysis/eco_matrix_region.xlsx")

data <- list()  # Create an empty list to store the data from each sheet

for (sheet in sheet_names) {
    data[[sheet]] <- read_excel("/mnt/d/PHD/ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sheet = sheet)
}

plot_path <- paste(HOME_, "/PHD/ISPRA_20152017_Analysis/Plots/Rich_levels/Acc_curve", sep = "")
dir.create(plot_path, showWarnings = FALSE)
for (region in names(data)) {
        data[[region]][, c(1:3)] <- data[[region]][, c(1:3)] %>% fill(Region, Season, id)
        incidence_freq_list <- list()
        incidence_freq_dataframe <- merge(
        data[[region]] %>% group_by(Season) %>% summarise(n_samples = n_distinct(id, Date)), 
        data[[region]] %>% select(-c(Region,id, Date)) %>% group_by(Season) %>% summarise_all(~sum(. != 0)), 
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

incidence_freq_list <- list()
for (region in names(data)) {
        data[[region]][, c(1:3)] <- data[[region]][, c(1:3)] %>% fill(Region, Season, id)
        incidence_freq_dataframe <- merge(
            data[[region]] %>% summarise(n_samples = n_distinct(id, Date)),
            data[[region]] %>% select(-c(Region,id, Date, Season)) %>% summarise_all(~sum(. != 0))
            )
        incidence_freq_dataframe <- incidence_freq_dataframe %>% .[. != 0]
        incidence_freq_list[[region]] <- incidence_freq_dataframe
}

res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = 200)
N <- length(incidence_freq_list)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
scale_fill_manual(values = unname(alphabet(N))) + 
scale_colour_manual(values = unname(alphabet(N))) + 
guides(shape=FALSE)

res <- iNEXT(incidence_freq_list, datatype = "incidence_freq", q = 1, conf = 0.95, endpoint = 200)
ggiNEXT(res, type=1, color.var="Assemblage") + theme_bw(base_size = 18) + 
scale_shape_manual(values=rep(19,length(incidence_freq_list))) + 
scale_fill_manual(values = unname(alphabet(N))) + 
scale_colour_manual(values = unname(alphabet(N))) + 
guides(shape=FALSE)

res$iNextEst$size_based %>% filter(between(SC, 0.90,0.95))






