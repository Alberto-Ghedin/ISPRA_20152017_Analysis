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

result <- multipatt(sites_taxa %>% select(-c(Date, id, Unknown)), basins[sites_taxa %>% pull(id), "Basin"], func = "IndVal.g", duleg = TRUE, restcomb = NULL)$sign

df <- result %>% filter(stat > 0.5, p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index)

result

index_to_basin <- c("Ion", "Lig","NordAdr", "SouthAdr", "SouthTyr", "WestMed")
index_to_basin
df %>% mutate(basin = index_to_basin[index]) %>% select(basin, stat, p.value) %>% write.csv(paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_per_basin.csv", sep = "/"), row.names = TRUE)

index_to_basin[1]
