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

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")
site_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/site_taxa_matrix.csv", sep = "/"))

params <- fromJSON(txt = paste(HOME_, "/ISPRA_20152017_Analysis/params.json", sep = "/"))
site_taxa$Date <- as.Date(site_taxa$Date, format = "%Y-%m-%d")
site_taxa$Region <- factor(site_taxa$Region, levels = params$ordered_regions, ordered = TRUE)
site_taxa$id <- factor(site_taxa$id, levels = params$ordered_id, ordered = TRUE)
site_taxa$Season <- factor(site_taxa$Season, levels = names(params$seasons), ordered = TRUE)



#site_taxa.dist <- vegdist(decostand(site_taxa[, -c(1:4)], method = "hellinger"), method = "euclidean")
#
## Hierarchical clustering
#site_taxa.hclust <- hclust(site_taxa.dist, method = "ward.D2")
#
#
#
#index_clusters <- site_taxa %>% select(id, Region, Season, Date) %>%
#  mutate(ward_9 = cutree(site_taxa.hclust, 9), 
#            ward_10 = cutree(site_taxa.hclust, 10), 
#            ward_11 = cutree(site_taxa.hclust, 11)
#  ) 

index_clusters <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))

indval_list <- list()
for (method in c("spectral_1.5_10", "ward_9", "ward_10", "ward_15", "ward_20")){
  result <- labdsv::indval(site_taxa %>% select(-c(1:4)), index_clusters[[method]])
  df <- as.data.frame(result$indval)
  df["pval"] <- result$pval
  df["max_val"] <- result$indcls
  indval_list[[method]] <- df
}

for (method in names(indval_list)) {
    indval_list[[method]] <- indval_list[[method]] %>% 
      filter(pval < 0.05 & max_val > 0.25)
    indval_list[[method]] <- indval_list[[method]] 
}

dev.off()
for (method in names(indval_list)) {
    png(paste(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_", sep = "/"), method, ".png", sep = ""), width = 1200, height = 1200)
    df <- indval_list[[method]] %>% filter(pval < 0.05 & max_val > 0.0) %>% select(-c(pval, max_val)) %>% mutate(Taxa = rownames(.)) %>% 
    pivot_longer(cols = -Taxa, names_to = "Cluster", values_to = "IndVal") 
    df$Cluster <- as.numeric(df$Cluster)
    p <- ggplot(df, aes(x = factor(Cluster), y = Taxa, fill = IndVal)) + 
    xlab("")+ 
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    labs(fill = "Z")
    plot(p)
    dev.off()
}


ggplot(df, aes(x = factor(Cluster), y = Taxa, fill = IndVal)) + 
    xlab("")+ 
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red")


list_results <- list()
for (method in names(indval_list)) {
  list_results[[method]] <- indval_list[[method]] %>% select(-c(pval, max_val)) %>% mutate(Taxa = rownames(.)) %>% 
pivot_longer(cols = -Taxa, names_to = "Cluster", values_to = "IndVal") %>% 
group_by(Taxa) %>% filter(IndVal == max(IndVal)) %>% ungroup() %>% arrange(Cluster)
}

write.xlsx(list_results, file = paste(HOME_, "ISPRA_20152017_Analysis/IndVal_per_cluster.xlsx", sep = "/"))
