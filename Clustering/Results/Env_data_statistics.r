library(jsonlite)
library(dplyr)
library(openxlsx)
library(vegan)
library(stats)
library(factoextra)
library(tidyverse)
library(rstatix)

# Load data
params <- fromJSON(txt = paste(path.expand("~"), "sys_specific.json", sep = "/"))
HOME_ <- paste(params$home, "PHD", sep = "/")

nutrients <- openxlsx::getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned.xlsx", sep = "/"))
nutrients <- nutrients[-which(nutrients == "E_cond")]
env_data <-  lapply(nutrients, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned.xlsx", sep = "/"), sheet = sheet, detectDate = TRUE))
names(env_data) <- nutrients


kepp_only_good_data <- TRUE
env_data_good <- vector("list", length = length(env_data))
if (kepp_only_good_data){
    for (nut in names(env_data)){
        env_data_good[[nut]] <- env_data[[nut]] %>% filter(QF == "1" | is.na(QF)) 
    }
}



for (nut in names(env_data)){
    env_data[[nut]]$Date <- env_data[[nut]]$Date %>% convertToDate()
    env_data[[nut]] <- env_data[[nut]] %>% group_by(Region, id, Date) %>% summarise(value = mean(Concentration, na.rm = TRUE), .groups = "drop") %>% rename(!!nut := value)
}
env_data <- Reduce(function(x, y) merge(x, y, by = c("id", "Date", "Region"), all = TRUE), env_data)


index_cluster <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Results/cluster_index.csv", sep = "/"))
index_cluster$Date <- as.Date(index_cluster$Date, format = "%Y-%m-%d")


env_data <- merge(env_data, index_cluster, by = c("Region", "id", "Date"))


diff_in_clusters <- function(method, nutrients) {
    kruskal <- setNames(seq(1:length(nutrients)), nutrients)
for (nut in nutrients){
    kruskal[nut] <- kruskal.test(as.formula(paste(nut, method, sep = " ~ " )) , data = env_data)$p.value
}


dunn_res <- vector("list", length(nutrients))
names(dunn_res) <- nutrients
for (nut in nutrients) {
   dunn_res[[nut]] <- rstatix::dunn_test(data = env_data, formula = as.formula(paste(nut, method, sep = " ~ ")), p.adjust.method = "bonferroni")  %>% 
select(all_of(c("group1", "group2", "p.adj"))) %>% mutate(pair = paste(group1, group2, sep = " - ")) %>% select(pair, p.adj) %>% rename(!!nut := p.adj) 

}

dunn_res <- Reduce(function(x, y) merge(x, y, by = "pair", all = TRUE), dunn_res) %>% column_to_rownames("pair")


cluster_statistic <- rbind(data.frame(t(kruskal), row.names = "kruskal"),dunn_res) %>% mutate_all(function (x) as.numeric(format(as.numeric(x), digits = 3)))
return(cluster_statistic)
}


selected_nutrinets <- c("NH4", "NO3", "TN", "PO4", "TP", "O_sat", "Salinity", "SiO4", "pH", "Chla", "T")
diff_in_clusters("spectral_1.5_9", selected_nutrinets)
cluster_statistic <- diff_in_clusters("ward_8", selected_nutrinets)
apply(cluster_statistic, 1, function(row) sum(row > 0.025, na.rm = TRUE))

data.frame(n_accepted = apply(cluster_statistic, 1, function(row) sum(row > 0.025, na.rm = TRUE)))

write.csv(cluster_statistic, paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Results/cluster_statistic.csv", sep = "/"), row.names = TRUE)




names(env_data)
pca_res <- prcomp(env_data %>% select(all_of(selected_nutrinets)) %>% na.omit() %>% na.omit(), scale = TRUE, center = TRUE)
eigs <- pca_res$sdev[c(1,2)]**2 / sum(pca_res$sdev**2) 



#groups <- env_data %>% select(all_of(c(selected_nutrinets, method))) %>% na.omit() %>% select(!!method) %>% pull() %>% as.factor()
fviz_pca_biplot(pca_res,
                #col.ind = groups,
                #geom.ind = "point",
               pointshape = 21,
                #fill.ind = groups,
                col.var = "black",
                repel = TRUE,
                #addEllipses = TRUE,
                label = c("ind", "var"),
                legend.title = "Clusters") + 
                labs(x = paste("PC1(", format(eigs[1] * 100, digits = 2), "%)", sep = ""), y = paste("PC2(", format(eigs[2] * 100, digits = 2), "%)", sep = ""), title = "PCA ENV") + 
                theme(plot.title = element_text(hjust = 0.5)) 

fviz_pca_biplot(pca_res, label = c("ind"), pointshape = 21, invisible = "var")
env_data %>%  summary()
method <- "ward_8"
dev.off()
for (index in unique(env_data[[method]])) {
    pca_res <- prcomp(env_data %>% filter((!!sym(method)) == index) %>% select(all_of(selected_nutrinets)) %>% na.omit(), scale = TRUE, center = TRUE)
    eigs <- pca_res$sdev[c(1,2)]**2 / sum(pca_res$sdev**2) 
    p <- fviz_pca_biplot(pca_res,
                geom.ind = "point",
                pointshape = 21,
                col.var = "black",
                repel = TRUE,
                addEllipses = TRUE,
                label = "var",
                legend.title = "Clusters") + 
                labs(x = paste("PC1(", format(eigs[1] * 100, digits = 2), "%)", sep = ""), y = paste("PC2(", format(eigs[2] * 100, digits = 2), "%)", sep = ""), title = paste("PCA ENV", index, sep = " - ")) + 
                theme(plot.title = element_text(hjust = 0.5)) 
    plot(p)
}
index <- 1
pca_res <- vegan::rda(env_data %>% filter((!!sym(method)) == index) %>% select(all_of(selected_nutrinets)) %>% na.omit() %>% slice(-c(415, 186)), scale = TRUE, center = TRUE)
eigs <- pca_res$CA$eig[c(1,2)] / sum(pca_res$CA$eig)
biplot(pca_res, scaling = 2, display = c("sites", "species"), type = c("text", "text"), cex = 0.7,
xlab = paste("PC1(", format(eigs[1] * 100, digits = 2), "%)", sep = ""),
ylab = paste("PC2(", format(eigs[2] * 100, digits = 2), "%)", sep = ""))

method <- "ward_8"
selected_nutrinets <- c("NH4", "NO3", "TN", "PO4", "TP", "O_sat", "Salinity", "SiO4", "pH", "Chla", "T")

pca_res <- vegan::rda(env_data %>% filter((!!sym(method)) == index) %>% select(all_of(selected_nutrinets)) %>% na.omit(), scale = TRUE, center = TRUE)
    eigs <- pca_res$CA$eig[c(1,2)] / sum(pca_res$CA$eig) 
biplot(pca_res, scaling = 2, display = c("sites", "species"), type = c("text", "text"), cex = 0.7,
xlab = paste("PC1(", format(eigs[1] * 100, digits = 2), "%)", sep = ""),
ylab = paste("PC2(", format(eigs[2] * 100, digits = 2), "%)", sep = ""))
for (index in unique(env_data[[method]])) {
    pca_res <- vegan::rda(env_data %>% filter((!!sym(method)) == index) %>% select(all_of(selected_nutrinets)) %>% na.omit(), scale = TRUE, center = TRUE)
    eigs <- pca_res$CA$eig[c(1,2)] / sum(pca_res$CA$eig) 
    print(index)
    p <- biplot(pca_res, scaling = 2, display = c("sites", "species"), type = c("text", "text"), cex = 0.7,
    xlab = paste("PC1(", format(eigs[1] * 100, digits = 2), "%)", sep = ""),
    ylab = paste("PC2(", format(eigs[2] * 100, digits = 2), "%)", sep = ""), title = index)
}
