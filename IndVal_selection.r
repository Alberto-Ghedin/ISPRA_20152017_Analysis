library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(labdsv)
library(openxlsx)
library(readxl)
library(tibble)
library(jsonlite)

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")

sheet_names <- excel_sheets(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"))

data <- list()  # Create an empty list to store the data from each sheet
for (sheet in sheet_names) {
    data[[sheet]] <- read_excel(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"), sheet = sheet)
}

# add season column 
seasons <- c(rep("Winter", 3), rep("Spring", 3), rep("Summer", 3), rep("Autumn", 3))
names(seasons) <- c(1:12)

for (region in names(data)) {
    data[[region]]$Season <- sapply(data[[region]]$Date, FUN = function(x) seasons[month(x)])
}


ind_val <- list()
for (region in names(data)){
    regional_matrix <- data[[region]]
    ind_val[[region]] <- labdsv::indval(regional_matrix %>% select(-c(Season, id, Date)), regional_matrix$Season)
}

list_dfs <- list()
for (region in names(ind_val)){
    df <- as.data.frame(ind_val[[region]]$indval)
    df["pval"] <- ind_val[[region]]$pval
    df["max_val"] <- ind_val[[region]]$indcls
    df <- tibble::rownames_to_column(df, var = "Taxa")
    list_dfs[[region]] <- df
}

file_path <- "ISPRA_20152017_Analysis"
openxlsx::write.xlsx(list_dfs, file = paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"))



# read excel file
sheet_names <- excel_sheets(paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"))
list_dfs <- lapply(sheet_names, function(sheet) {
  read_excel(paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"), sheet = sheet)
})
names(list_dfs) <- sheet_names

results <- lapply(list_dfs, function(x) {
  data.frame(
    pval_less_than_0_05 = sum(x$pval < 0.05),
    indval_greater_than_0_25 = sum(x$max_val > 0.25),
    both_conditions = sum(x$pval < 0.05 & x$max_val > 0.25)
  )
})

results_df <- do.call(rbind, results)


selected_species <- lapply(list_dfs, function(df) {
  df %>% filter(pval < 0.05 & max_val > 0.25) %>% select(Taxa) %>% pull()
})


selected_species <- unique(unname(unlist(selected_species, recursive = FALSE)))

write(selected_species, file = paste(HOME_, file_path, "selected_species.txt", sep = "/"))

paste(HOME_, file_path, "selected_species.txt", sep = "/")
