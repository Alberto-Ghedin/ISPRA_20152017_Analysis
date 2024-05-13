library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(labdsv)
library(openxlsx)
library(readxl)
library(tibble)
library(jsonlite)
library(zoo)

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")

sheet_names <- excel_sheets(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"))

data <- list()  # Create an empty list to store the data from each sheet
for (sheet in sheet_names) {
    data[[sheet]] <- read_excel(paste(HOME_, "ISPRA_20152017_Analysis/eco_matrix_region.xlsx", sep = "/"), sheet = sheet)
    data[[sheet]][, "Region"] <- sheet 
}

# Combine all the data into a single data frame, columns do not match 
data <- do.call(bind_rows, data) %>% mutate(id = na.locf(id), Season = na.locf(Season))
data <- data %>% relocate(Region) 
data <- data %>% replace(is.na(.), 0)

#load cluster index
cluster_index <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))
cluster_index["Date"] <- as.Date(cluster_index$Date, format = "%Y-%m-%d")

data <- data %>% merge(cluster_index, by = c("Date" = "Date", "id" = "id")) %>% relocate(ward, spectral_1.5_10, .after = Season)


data.indval <- vector("list", 2)
names(data.indval) <- c("ward", "spectral_1.5_10")
for (method in names(data.indval)) {
    data.indval[[method]] <- labdsv::indval(data %>% select(-c(Season, id, Date, Region, ward, spectral_1.5_10)), data[, method])
}

list_dfs <- list()
for (method in names(data.indval)){
    df <- as.data.frame(data.indval[[method]]$indval)
    df["pval"] <- data.indval[[method]]$pval
    df["max_val"] <- data.indval[[method]]$indcls
    df <- tibble::rownames_to_column(df, var = "Taxa")
    list_dfs[[method]] <- df
}


selected_species <- lapply(list_dfs, function(df) {
  df %>% filter(pval < 0.05 & max_val > 0.25)
})



#compute PCA for each cluster using selected species
pca_results <- vector("list", 2)
names(pca_results) <- c("ward", "spectral_1.5_10")
for (method in names(pca_results)) {
    selected_species[[method]]$Taxa <- as.character(selected_species[[method]]$Taxa)
    
    #loop over clusters
    for (cluster in unique(data[, method])) {
        cluster_data <- data %>% filter(data[, method] == cluster)
        cluster_data <- cluster_data %>% select(c("id", "Date", selected_species[[method]]$Taxa)) %>%  select_if(~ !is.numeric(.) || sum(.) != 0)
        is_all_zero <- rowSums(cluster_data[-c(1,2)]) == 0
        cluster_data <- cluster_data[!is_all_zero,]
        pca_results[[method]][[as.character(cluster)]] <- rda(decostand(cluster_data[, -c(1,2)], "hellinger"), scale = FALSE)
        }
}

for (method in names(pca_results)) {
    for (cluster in names(pca_results[[method]])) {
        #computed explained variance by first two axes
        pca <- pca_results[[method]][[cluster]]
        explained_variance <- sum(pca$CA$eig[c(1,2)]) / sum(pca$CA$eig)
        print(paste(method, cluster, explained_variance, sep = " "))
    }
}

summary(pca_results[["ward"]][["0"]])$cont
pca_results[["ward"]][["0"]]$CA$eig[c(1,2)]

sum(pca_results[["ward"]][["0"]]$CA$eig[c(1,2)])
sum(pca_results[["ward"]][["0"]]$CA$eig)
plot_path <- paste(HOME_, "ISPRA_20152017_Analysis/Clustering/PCA", sep = "/")
for (method in names(pca_results)) {
    dir.create(paste(plot_path, method, sep = "/"), recursive = TRUE, showWarnings = FALSE)

}

for (method in names(pca_results)) {
    for (cluster in names(pca_results[[method]])) {
        pca <- pca_results[[method]][[cluster]]
        explained_variance <- lapply(pca$CA$eig[c(1,2)], function(x) {
            x / sum(pca$CA$eig)
        }
        )
        png(paste(plot_path, method, paste("biplot_", cluster, ".png", sep = ""), sep = "/"), width = 1500, height = 960, res = 100)
        par(mfrow = c(1, 2))
        biplot(pca_results[[method]][[as.character(cluster)]], scaling = 1, type = "text", display = "sites", cex = 0.7, col = "blue", main = "scaling 1", xlab = "", ylab = "")
        title(xlab=paste("PC1", round(explained_variance$PC1, 2), sep = " "), ylab=paste("PC2", round(explained_variance$PC2, 2), sep = " "), mgp=c(2.2, 2.2, 0))
        biplot(pca_results[[method]][[as.character(cluster)]], scaling = 2, type = "text", cex = 0.7, display = "species", col = "blue", main = "scaling 2", xlab = "", ylab = "")
        title(xlab=paste("PC1", round(explained_variance$PC1, 2), sep = " "), ylab=paste("PC2", round(explained_variance$PC2, 2), sep = " "), mgp=c(2.2, 2.2, 0))
        dev.off()
    }
}

pca <- pca_results[[method]][[cluster]]
explained_variance <- lapply(pca$CA$eig[c(1,2)], function(x) {
    x / sum(pca$CA$eig)
}
)

ggplot() +
    geom_point(data = scores(test), aes(x = PC1, y = PC2)) +
    geom_text(data = scores(test), aes(x = PC1, y = PC2, label = rownames(scores(test))), size = 3) +
    labs(x = "PC1", y = "PC2") +
    theme_minimal()

test <- rda(decostand(cluster_data[, -c(1,2)], "hellinger"), scale = FALSE, scaling = 1)

with(decostand(cluster_data[, -c(1,2)], "hellinger"), points(test, scaling = 1, display = "sites"))

ggplot(as.data.frame(test$CA$u), aes(x = PC1, y = PC2)) +
    geom_point()

rownames(as.data.frame(test$CA$v))
ggplot() +  geom_text(data = as.data.frame(test$CA$v), aes(x = PC1, y = PC2, label = rownames(as.data.frame(test$CA$v))), col = 'red') +
geom_segment(data = as.data.frame(test$CA$v), aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
            alpha = 0.75, color = 'darkred')
as.data.frame(test$CA$u)
autoplot(test, scaling = 1, display = "sites")
