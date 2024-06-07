library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(labdsv)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(openxlsx)
library(parallel)
library(indicspecies)

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")
site_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/site_taxa_matrix.csv", sep = "/"))

params <- fromJSON(txt = paste(HOME_, "/ISPRA_20152017_Analysis/params.json", sep = "/"))
site_taxa$Date <- as.Date(site_taxa$Date, format = "%Y-%m-%d")
site_taxa$Region <- factor(site_taxa$Region, levels = params$ordered_regions, ordered = TRUE)
site_taxa$id <- factor(site_taxa$id, levels = params$ordered_id, ordered = TRUE)
site_taxa$Season <- factor(site_taxa$Season, levels = names(params$seasons), ordered = TRUE)


apply_log <- TRUE
if (apply_log) {
  site_taxa[, -c(1:4)] <- log(site_taxa %>% select(-c(1:4)) + 1, 10)
}


index_clusters <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))
methods <- index_clusters %>% select(contains("ward"), contains("spectral")) %>% names()



compute_indval <- function(method, eco_matrix, duleg = TRUE, restcomb = NULL) {
  result <- multipatt(site_taxa[, -c(1:4)], index_clusters[[method]], func = "IndVal.g", duleg = duleg, restcomb = restcomb)$sign
  df <- result %>% filter(stat > 0.5, p.value < 0.05) %>% select(index, stat, p.value) %>% arrange(index)
  return(df)
}

## INDVAL log all_vs_sall (ava)
indval_list <- list()
indval_list <- mclapply(methods, compute_indval, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"), rowNames = TRUE)


## INDVAL log one_vs_sall (ova)
indval_list <- list()
f <- function(method, eco_matrix) {
  n_clusters <- as.numeric(tail(strsplit(method, "_")[[1]], n = 1))
  print(n_clusters)
  compute_indval(method, eco_matrix = eco_matrix, duleg = FALSE, restcomb = c(1:n_clusters))
  print(paste("finished", method))
}
indval_list <- mclapply(methods, f, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ova.xlsx", sep = "/"), rowNames = TRUE)

## INDVAL no_log all_vs_all (ava)
indval_list <- list()
indval_list <- mclapply(methods, compute_indval, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_no_log_ava.xlsx", sep = "/"), rowNames = TRUE)

## INDVAL no_log one_vs_all (ova)
indval_list <- list()
f <- function(method, eco_matrix) {
  n_clusters <- as.numeric(tail(strsplit(method, "_")[[1]], n = 1))
  print(n_clusters)
  df <- compute_indval(method, eco_matrix = eco_matrix, duleg = FALSE, restcomb = c(1:n_clusters))
  print(paste("finished", method))
  return(df)
}
indval_list <- mclapply(methods, f, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_no_log_ova.xlsx", sep = "/"), rowNames = TRUE)


sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"))

indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"), sheet = sheet))
names(indval_list) <- sheets

species_per_cluster <- list()
for (method in c("ward", "spectral")) {
  cols <- names(indval_list)[grep(method, names(indval_list))]
  n_clusters <- c()
  n_species <- c()
  for (col in cols) {
    n_clusters <- c(n_clusters, as.numeric(tail(strsplit(col, "_")[[1]], n = 1)))
    n_species <- c(n_species, nrow(indval_list[[col]]))
  }
  species_per_cluster[[method]] <- data.frame(n_clusters = n_clusters, n_species = n_species)
}


df <- do.call(rbind, species_per_cluster)
df$method <- rep(names(species_per_cluster), sapply(c("ward", "spectral"), function(x) {length(grep(x, names(indval_list)))}, simplify = "array"))


p<- ggplot(df, aes(x = n_clusters, y = n_species, color = method, group = method)) +
  geom_line() +
  labs(x = "N_clusters", y = "N_species", title = "Total number of indicator species") +
  scale_x_continuous(breaks = seq(min(df$n_clusters), max(df$n_clusters), 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
#save plot as png
ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/Species_per_cluster_no_log_ava.pdf", sep = "/"), plot = p, width = 10, height = 5)










## try species combinations

phyto_abund <- site_taxa %>% 
  pivot_longer(
    cols = -c(Date, Region, id, Season), 
    names_to = "Taxon", 
    values_to = "Ind.L"
  )

phyto_abund %>% filter(Ind.L != 0) %>% group_by(Taxon) %>% summarise(count = n()) %>% arrange(desc(count)) %>% filter(count >= quantile(count, 0.95)) %>% pull(Taxon) -> taxa

phyto_abund %>% filter(Ind.L != 0) %>% group_by(Taxon) %>% summarise(count = mean(Ind.L)) %>% arrange(desc(count)) %>% filter(count >= quantile(count, 0.95)) %>% pull(Taxon) -> taxa

length(taxa)
dim(site_taxa)
any(site_taxa %>% select(all_of(taxa)) %>% rowSums() == 0)

choose(39, 3)

taxa
