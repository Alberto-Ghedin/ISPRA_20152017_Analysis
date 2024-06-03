library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
#library(labdsv)
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

site_taxa[, -c(1:4)] <- log(site_taxa %>% select(-c(1:4)) + 1, 10)


index_clusters <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))


compute_indval <- function(method, eco_matrix) {
  result <- multipatt(site_taxa[, -c(1:4)], index_clusters[[method]], func = "IndVal.g", duleg = TRUE)$sign
  df <- result %>% filter(stat > 0.5, p.value < 0.05) %>% select(index, stat, p.value) 
  return(df)
}

methods_ward <- index_clusters %>% select(contains("ward")) %>% names()


indval_list <- list()
indval_list <- mclapply(methods_ward[-1], compute_indval, eco_matrix = site_taxa[, -c(1:4)], mc.cores = 4)
names(indval_list) <- methods_ward[-1]

#save indval_list as excel file
write.xlsx(indval_list, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_on_hclust.xlsx", sep = "/"), rowNames = TRUE)

species_per_cluster <- list()
for (method in c("ward")) {
  cols <- names(indval_list)[grep(method, names(indval_list))]
  n_clusters <- c()
  n_species <- c()
  for (col in cols) {
    n_clusters <- c(n_clusters, as.numeric(strsplit(col, "_")[[1]][2]))
    n_species <- c(n_species, nrow(indval_list[[col]]))
  }
  species_per_cluster[[method]] <- data.frame(n_clusters = n_clusters, n_species = n_species)
}

df <- do.call(rbind, species_per_cluster)
#df$method <- rep(names(species_per_cluster), each = nrow(df) / 2)

p<- ggplot(df, aes(x = n_clusters, y = n_species)) +
  geom_line() +
  labs(x = "N_clusters", y = "N_species") +
  scale_x_continuous(breaks = seq(min(df$n_clusters), max(df$n_clusters), 1)) +
  theme_minimal()

#save plot as png
ggsave(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/Species_per_cluster.pdf", sep = "/"), plot = p, width = 10, height = 5)

h_clust_sklearn <- as.matrix(read.table(paste(HOME_ , "ISPRA_20152017_Analysis/Clustering/Results/link_matrix_to_hclust.txt", sep = "/"), sep = " ", header = FALSE))



linkage

a$merge <- h_clust_sklearn[, 1:2]

a$height <- h_clust_sklearn[, 3]


a$order <- as.list(read.table(paste(HOME_ , "ISPRA_20152017_Analysis/Clustering/Results/leaves.txt", sep = "/"), sep = " ", header = FALSE))$V1
length(a$order)
a$labels <- 1:(nrow(a$merge) +1)
class(a) <- "hclust" 

a$merge
clustering <- hclust(hell_dist, method = "ward.D2")

cut <- cutree(a, k = 12)
plot(a, ylim = c(10, max(a$height) + 0.1)) 
abline(h = max(), col = "red")

a$height[cut == 12]
plot(cutree(a, h = 0.86))
a$height
cut == 12
table(labels)

cutree(a, h = 0.86)
dev.off()
for (method in names(indval_list)) {
    png(paste(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_", sep = "/"), method, ".png", sep = ""), width = 1200, height = 1200, res = 200)
    df <- indval_list[[method]] %>% filter(max_val > 0.25) %>% select(-c(pval, max_val)) %>% mutate(Taxa = rownames(.)) %>% 
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






df <- do.call(rbind, species_per_cluster)


# Plot the data
ggplot(df, aes(x = n_clusters, y = n_species, color = method, group = method)) +
  geom_line() +
  labs(x = "N_clusters", y = "N_species") +
  scale_x_continuous(breaks = seq(min(df$n_clusters), max(df$n_clusters), 1)) +
  theme_minimal()




list_results <- list()
for (method in names(indval_list)) {
  list_results[[method]] <- indval_list[[method]] %>% select(-c(pval, max_val)) %>% mutate(Taxa = rownames(.)) %>% 
pivot_longer(cols = -Taxa, names_to = "Cluster", values_to = "IndVal") %>% 
group_by(Taxa) %>% filter(IndVal == max(IndVal)) %>% ungroup() %>% arrange(Cluster)
}

write.xlsx(list_results, file = paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/IndVal_per_cluster.xlsx", sep = "/"))
