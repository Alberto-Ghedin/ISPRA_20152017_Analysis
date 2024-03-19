library(vegan)


# Load the data
HOME_ <- "/mnt/d"
phyto_file_path <- "/PHD/ISPRA_20152017_Analysis/site_taxa_matrix.csv"
phyto.matrix <- read.csv(paste(HOME_, phyto_file_path, sep=""))

# Split the dataframe into two
phyto.sites_info <- phyto.matrix[, c(1,2,3)]
phyto.matrix <- phyto.matrix[, -c(1,2,3)]

phyto.hell <- vegan::decostand(phyto.matrix, method="hellinger")

# 1. Calculate the distance matrix
phyto.dist <- vegdist(phyto.hell)

phyto.hell[1, c(1:20)]
rm(clusters)

clusters.ward <- hclust(phyto.dist, method="ward.D2")
clusters.singlelink <- hclust(phyto.dist, method="single")

plot(clusters.ward)
plot(clusters.singlelink)

plot(cutree(clusters, k = 20))

plot(as.dendrogram(clusters))
rect.hclust(clusters, k=10, border="red")

dev.off()
plot(
    clusters.ward$height, 
    nrow(phyto.hell):2
)
plot(
    clusters.singlelink$height, 
    nrow(phyto.hell):2
)
