library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(labdsv)
library(openxlsx)
library(readxl)
library(tibble)


HOME_ <- "/mnt/d"
file_path <- "PHD/ISPRA_20152017_Analysis"
file_name <- "site_taxa_matrix.csv"
sites_taxa = read.csv(paste(HOME_, file_path, file_name, sep="/"), check.names = FALSE)


#make column of date objects
sites_taxa["Date"] <- as.Date(sites_taxa$Date)

# add season column 
seasons <- c(rep("Winter", 3), rep("Spring", 3), rep("Summer", 3), rep("Autumn", 3))
names(seasons) <- c(1:12)
sites_taxa$Season <- sapply(sites_taxa$Date, FUN = function(x) seasons[month(x)])

sites_taxa <- sites_taxa %>% relocate(Season, .after = Date)


regional_matrix <- list()
for (region in unique(sites_taxa$Region)){
  regional_matrix[[region]] <- sites_taxa %>% filter(Region == region) %>% select(-c(Region, Date, id))
  regional_seasons <- regional_matrix[[region]]$Season
  regional_matrix[[region]] <- regional_matrix[[region]][, -1] %>% select_if(~sum(.) != 0)
  regional_matrix[[region]]$Season <- regional_seasons
}

# calculate indicator value
ind_val <- list()
for (region in unique(sites_taxa$Region)){
  ind_val[[region]] <- labdsv::indval(regional_matrix[[region]] %>% select(-c(Season)), regional_matrix[[region]]$Season)
}

list_dfs <- list()
for (region in names(ind_val)){
    df <- as.data.frame(ind_val[[region]]$indval)
    df["pval"] <- ind_val[[region]]$pval
    df["max_val"] <- ind_val[[region]]$indcls
    df <- rownames_to_column(df, var = "Taxa")
    list_dfs[[region]] <- df
}

write.xlsx(list_dfs, file = paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"))



# read excel file
sheet_names <- excel_sheets(paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"))
list_dfs <- lapply(sheet_names, function(sheet) {
  read_excel(paste(HOME_, file_path, "IndVal_per_region.xlsx", sep = "/"), sheet = sheet, chekc.names = FALSE)
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

