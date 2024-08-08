library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(openxlsx)
library(parallel)
library(indicspecies)
library(gridExtra)


#getting abund values
HOME_ <- paste(paste(path.expand("~"), "PHD", sep = "/"))
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/sites_taxa.csv", sep = "/"))
params <- fromJSON(txt = paste(HOME_, "/ISPRA_20152017_Analysis/params.json", sep = "/"))


#getting cluster indices
sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/cluster_indices.xlsx", sep = "/"))
cluster_methods <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/cluster_indices.xlsx", sep = "/"), sheet = sheet))
names(cluster_methods) <- sheets



compute_indval <- function(method, eco_matrix, index_clusters, duleg = TRUE, restcomb = NULL, stat_threshold = 0) {
  result <- multipatt(eco_matrix, index_clusters[[method]], func = "IndVal.g", duleg = duleg, restcomb = restcomb)$sign
  df <- result %>% filter(stat > stat_threshold, p.value < 0.05) %>% select(index, stat, p.value) %>% arrange(index)
  return(df)
}

wrapper <- function(method, eco_matrix, index_clusters) {
  n_clusters <- as.numeric(tail(strsplit(method, "_")[[1]], n = 1))
  print(paste("starting", n_clusters))
  result <- compute_indval(method, eco_matrix = eco_matrix, index_clusters, duleg = TRUE, stat_threshold = 0)
  print(paste("finished", method))
  return(result)
}

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
  ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect", plot_name, sep = "/"), plot = p, width = 10, height = 5)
}
if (return_plot) {
  return(p)
}
  }

#### METHOD 1 ###
indval_list <- list()
indval_list <- mclapply(
  names(cluster_methods[["method_1"]])[-c(1,2)], 
  function (methods) wrapper(methods, eco_matrix = sites_taxa %>% select(-c(Date, id, Unknown)), index_cluster = cluster_methods[["method_1"]]),   
  mc.cores = 2
)
names(indval_list) <- names(cluster_methods[["method_1"]])[-c(1,2)]
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_1.xlsx", sep = "/"), rowNames = TRUE)


sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_1.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_1.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
names(indval_list) <- sheets


args_list <- list(
  list(indval_list, "ward", length, "Total number of indicator species", "N_species", return_plot = TRUE),
  list(indval_list, "ward", length, "Total number of indicator species above 0.5", "N_species", 0.5, return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (species above 0.5)", "IndVal", 0.5, return_plot = TRUE), 
  list(indval_list, method = "ward", mean, "Mean of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, methods = "ward",mean, plot_title = "Mean of IndVal (species above 0.5)", y_label = "IndVal", threshold = 0.5, return_plot = TRUE)
)


plots <- sapply(args_list, function(args) {
  do.call(plot_indval_statistic, args)
}, simplify = FALSE, USE.NAMES = FALSE)
condensed_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)

ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_statistic_method_1.pdf", sep = "/"), plot = condensed_plot, width = 12, height = 10)



#### METHOD 2 ###
indval_list <- list()
indval_list <- mclapply(
  names(cluster_methods[["method_2"]])[-c(1,2)], 
  function (methods) wrapper(methods, eco_matrix = sites_taxa %>% select(-c(Date, id, Unknown)), index_cluster = cluster_methods[["method_2"]]),   
  mc.cores = 2
)
names(indval_list) <- names(cluster_methods[["method_2"]])[-c(1,2)]
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"), rowNames = TRUE)


sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
names(indval_list) <- sheets


args_list <- list(
  list(indval_list, "ward", length, "Total number of indicator species", "N_species", return_plot = TRUE),
  list(indval_list, "ward", length, "Total number of indicator species above 0.5", "N_species", 0.5, return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (species above 0.5)", "IndVal", 0.5, return_plot = TRUE), 
  list(indval_list, method = "ward", mean, "Mean of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, methods = "ward",mean, plot_title = "Mean of IndVal (species above 0.5)", y_label = "IndVal", threshold = 0.5, return_plot = TRUE)
)


plots <- sapply(args_list, function(args) {
  do.call(plot_indval_statistic, args)
}, simplify = FALSE, USE.NAMES = FALSE)
condensed_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)

ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_statistic_method_2.pdf", sep = "/"), plot = condensed_plot, width = 12, height = 10)



#### METHOD 3 ###
indval_list <- list()
indval_list <- mclapply(
  names(cluster_methods[["method_3"]])[-c(1,2)], 
  function (methods) wrapper(methods, eco_matrix = sites_taxa %>% select(-c(Date, id)), index_cluster = cluster_methods[["method_3"]]),   
  mc.cores = 2
)
names(indval_list) <- names(cluster_methods[["method_3"]])[-c(1,2)]
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_3.xlsx", sep = "/"), rowNames = TRUE)


sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_3.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_3.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
names(indval_list) <- sheets


args_list <- list(
  list(indval_list, "ward", length, "Total number of indicator species", "N_species", return_plot = TRUE),
  list(indval_list, "ward", length, "Total number of indicator species above 0.5", "N_species", 0.5, return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, "ward", sum, "Sum of IndVal (species above 0.5)", "IndVal", 0.5, return_plot = TRUE), 
  list(indval_list, method = "ward", mean, "Mean of IndVal (all species)", "IndVal", return_plot = TRUE),
  list(indval_list, methods = "ward",mean, plot_title = "Mean of IndVal (species above 0.5)", y_label = "IndVal", threshold = 0.5, return_plot = TRUE)
)


plots <- sapply(args_list, function(args) {
  do.call(plot_indval_statistic, args)
}, simplify = FALSE, USE.NAMES = FALSE)
condensed_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)

ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_statistic_method_3.pdf", sep = "/"), plot = condensed_plot, width = 12, height = 10)
