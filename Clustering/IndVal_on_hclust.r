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
library(gridExtra)

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


compute_indval <- function(method, eco_matrix, duleg = TRUE, restcomb = NULL, stat_threshold = 0) {
  result <- multipatt(site_taxa[, -c(1:4)], index_clusters[[method]], func = "IndVal.g", duleg = duleg, restcomb = restcomb)$sign
  df <- result %>% filter(stat > stat_threshold, p.value < 0.05) %>% select(index, stat, p.value) %>% arrange(index)
  return(df)
}

## INDVAL log all_vs_all (ava)
f <- function(method, eco_matrix) {
  n_clusters <- as.numeric(tail(strsplit(method, "_")[[1]], n = 1))
  print(paste("starting", n_clusters))
  result <- compute_indval(method, eco_matrix = eco_matrix, duleg = TRUE, stat_threshold = 0)
  print(paste("finished", method))
  return(result)
}
indval_list <- list()
indval_list <- mclapply(methods, f, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"), rowNames = TRUE)


## INDVAL log all combinations

f <- function(method, eco_matrix) {
  n_clusters <- as.numeric(tail(strsplit(method, "_")[[1]], n = 1))
  print(paste("starting", n_clusters))
  result <- compute_indval(method, eco_matrix = eco_matrix, duleg = FALSE)
  print(paste("finished", method))
  return(result)
}

subset <- c(1:6)
indval_list <- list()
indval_list <- mclapply(methods[grep("spectral", methods)][subset], f, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods[grep("spectral", methods)][subset]
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_1_6.xlsx", sep = "/"), rowNames = TRUE)

indval_list["spectral_1.5_4"]
subset <- c(7:12)
indval_list <- list()
indval_list <- mclapply(methods[grep("spectral", methods)][subset], f, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods[grep("spectral", methods)][subset]
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_7_12.xlsx", sep = "/"), rowNames = TRUE)



## INDVAL no_log all_vs_all (ava)
indval_list <- list()
indval_list <- mclapply(methods, compute_indval, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_no_log_ava.xlsx", sep = "/"), rowNames = TRUE)


sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_log_ava.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
names(indval_list) <- sheets


args_list <- list(
  list(indval_list, c("ward", "spectral"), length, "Total number of indicator species", "N_species", return_plot = TRUE),
  list(indval_list, c("ward", "spectral"), length, "Total number of indicator species above 0.5", "N_species", 0.5, return_plot = TRUE),
  list(indval_list, c("ward", "spectral"), sum, "Sum of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, c("ward", "spectral"), sum, "Sum of IndVal (species above 0.5)", "IndVal", 0.5, return_plot = TRUE), 
  list(indval_list, method = c("ward", "spectral"), mean, "Mean of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, methods = c("ward", "spectral"),mean, plot_title = "Mean of IndVal (species above 0.5)", y_label = "IndVal", threshold = 0.5, return_plot = TRUE)
)


plot_indval_statistic <- function(indval, methods, statistic, plot_title, y_label, threshold = 0, return_plot = FALSE, plot_name = NULL) {
  stat_per_cluster <- list()
  for (method in methods) {
    cols <- names(indval_list)[grep(method, names(indval_list))]
    n_clusters <- sapply(cols, function(col) {as.numeric(tail(strsplit(col, "_")[[1]], n = 1))}, USE.NAMES = FALSE)
    stat_per_cluster[[method]] <- do.call(rbind, mapply(function(df, n) {
    df$n_clusters <- n
    return(df)
  }, indval_list[cols], n_clusters, SIMPLIFY = FALSE))
  }
  df <- do.call(rbind, stat_per_cluster)
  df$method <- rep(names(stat_per_cluster), sapply(stat_per_cluster, nrow), simplify = "array")
  df <- df %>% filter(stat >= threshold) %>% group_by(method, n_clusters) %>%
    summarise(stat_column = {{statistic}}(stat), .groups = "drop")
  p<- ggplot(df, aes(x = n_clusters, y = stat_column, color = method, group = method)) +
  geom_line() +
  labs(x = "N_clusters", y = y_label, title = plot_title) +
  scale_x_continuous(breaks = seq(min(df$n_clusters), max(df$n_clusters), 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), plot.background = element_rect(fill = 'white'))
#save plot as png
if (!is.null(plot_name)) {
  ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results", plot_name, sep = "/"), plot = p, width = 10, height = 5)
}
if (return_plot) {
  return(p)
}
  }

plots <- sapply(args_list, function(args) {
  do.call(plot_indval_statistic, args)
}, simplify = FALSE, USE.NAMES = FALSE)
# Arrange the plots
condensed_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)

ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_statistic_plots.png", sep = "/"), plot = condensed_plot, width = 12, height = 10)


metric <- "n"
print(df %>% group_by(method, n_clusters) %>% length(.), n = 40)

df %>% filter(n_clusters <= 2)


df %>% group_by(method, n_clusters) %>% summarise(stat_column = sum(stat))
df %>% group_by(method, n_clusters)
values <- sapply(cols, function(col) {indval_list[[col]] %>% filter(stat >= threshold) %>% pull(stat) %>% {{statistic}} }, simplify = "array", USE.NAMES = FALSE)
method
cols





## ALL COMBINATIONS SPECTRAL
#merge all indval results
sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_1_6.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_1_6.xlsx", sep = "/"), sheet = sheet))
names(indval_list) <- sheets
sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_7_12.xlsx", sep = "/"))
temp <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_spectral_combinations_7_12.xlsx", sep = "/"), sheet = sheet)) 
names(temp) <- sheets
indval_list <- c(indval_list, temp)

cols <- names(indval_list)
df <- data.frame()
for (col in cols) {
  n_clusters <- as.numeric(tail(strsplit(col, "_")[[1]], n = 1))
  df <- rbind(df, indval_list[[col]] %>% select(index, stat) %>% mutate(n_clusters = n_clusters))

} 


species_per_cluster["spectral"]
df_merged <- merge(
  merge(df %>% group_by(n_clusters) %>% summarise(n = n()), 
  df %>% filter(index <= n_clusters) %>% group_by(n_clusters) %>% summarise(n_single = n()),
  by = "n_clusters"), 
  species_per_cluster[["spectral"]], 
  by = "n_clusters"
  )


p<- ggplot(df_merged) +
  geom_line(aes(x = n_clusters, y = n), color = "blue") +
  geom_line(aes(x = n_clusters, y = n_species), color = "red") +
  geom_line(aes(x = n_clusters, y = n_single), color = "green") +
  labs(x = "N_clusters", y = "N_species", title = "Total number of indicator species") +
  scale_x_continuous(breaks = seq(min(df$n_clusters), max(df$n_clusters), 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))
plot(p)



df %>% filter(n_clusters == 2)

ava <- indval_list["spectral_1.5_2"]
merge(indval_list["spectral_1.5_2"]


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
