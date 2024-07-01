library(jsonlite)
library(dunn.test)
library(dplyr)
library(openxlsx)

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")

nutrients <- openxlsx::getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned.xlsx", sep = "/"))
nutrients <- nutrients[-which(nutrients == "E_cond")]
env_data <-  lapply(nutrients, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned.xlsx", sep = "/"), sheet = sheet))
names(env_data) <- nutrients

index_cluster <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))

index_cluster
method <- "ward_6"
env_data <- lapply(
    env_data,
    function(df_nut) merge(df_nut, index_cluster[, c("id", "Date", method)],  
    by = c("id", "Date"), 
    all.x = TRUE
    )
)
names(env_data) <- nutrients



