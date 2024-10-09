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
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/sites_taxa.csv", sep = "/"))
basins <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Description/Station_basin.csv", sep = "/"))

names(basins)

#set as index of basin id
rownames(basins) <- basins$id

sites_taxa %>% head()

result <- multipatt(sites_taxa %>% select(-c(Date, id, Unknown)), basins[sites_taxa %>% pull(id), "Basin"], func = "IndVal.g", duleg = TRUE, restcomb = NULL)[c("str", "sign")]
mat <- multipatt(sites_taxa %>% select(-c(Date, id, Unknown)), basins[sites_taxa %>% pull(id), "Basin"], func = "IndVal.g", duleg = TRUE, restcomb = NULL)$str

taxa <- result$sign %>% filter(stat > 0.5) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()
result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))) %>% 
    select(c(NordAdr, SouthAdr, Ion, SouthTyr, Lig, WestMed)) %>%
    write.csv(paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_along_basins.csv", sep = "/"), row.names = TRUE)
df <- result %>% filter(stat > 0.5, p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index)



test <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_along_basins.csv", sep = "/"), row.names =1)
