library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(openxlsx)
library(parallel)
library(indicspecies)


HOME_ <- paste(paste(path.expand("~"), "PHD", sep = "/"))
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/sites_taxa.csv", sep = "/"))
basins <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Description/Stations_Basins.csv", sep = "/"))
rownames(basins) <- basins$id


result <- multipatt(
    sites_taxa %>% select(-c(Date, id)), 
    basins[sites_taxa %>% pull(id), "Basin"], 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5,6, 14, 15, 16), #considering also combination of basis 14 = Lig+STyr, 15 = Lig+WMed, 16 = NA+SA
    control = how(nperm = 999) #duleg = TRUE
    )

taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_per_basin.xlsx", sep = "/"), rowNames = TRUE)



result <- multipatt(
    sites_taxa %>% select(-c(Date, id)) %>% mutate(across(everything(), ~ log1p(.))), 
    basins[sites_taxa %>% pull(id), "Basin"], 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5,6, 14, 15, 16), #considering also combination of basis 14 = Lig+STyr, 15 = Lig+WMed, 16 = NA+SA
    control = how(nperm = 999) #duleg = TRUE
    )

taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_per_basin_log_trasf.xlsx", sep = "/"), rowNames = TRUE)