library(dplyr)
library(tibble)
library(iNEXT)
# Set the seed for reproducibility
set.seed(123)

# Generate 300 probabilities from a negative exponential distribution
probabilities <- dexp(seq(0, 10, length.out = 300), rate =2.5)
probabilities <- probabilities/sum(probabilities)
plot(sort(probabilities))
probabilities
n_individuals <- 1000
length(probabilities)

#sample n_individuals from a pool of 300 elements with the given probabilities

phyto_abund <- read.csv("phyto_abund.csv", header = TRUE)
sample_abund_list <- phyto_abund %>% filter(Det_level == "Species") %>% 
    group_by(Region, Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l))) %>% 
    group_by(Region) %>% summarise(sample_abund = list(sample_abund)) %>% deframe()

n_samples_list <- sapply(sample_abund_list, length)
n_samples_df <- data.frame(Region = names(n_samples_list), n_samples = n_samples_list)
probabilities <- phyto_abund %>% filter(Det_level == "Species") %>% group_by(Taxon) %>% summarise(n_obs = n_distinct(Date, id)) %>% pull(n_obs)
probabilities <- probabilities/sum(probabilities)
n_pool <- length(probabilities)

unique_observations <- list()
for (index in seq_along(1:nrow(n_samples_df))) {
    region <- n_samples_df$Region[index]
    n_samples <- n_samples_df$n_samples[index]
    
    observed_taxa <- c()
    for (n_individuals in sample_abund_list[[region]]) {
        observed_taxa <- c(
            observed_taxa,
            sample(1:n_pool, n_individuals, replace = TRUE, prob = probabilities) %>% unique()
        )
    }
    unique_observations[[region]] <- c(n_samples, unname(table(observed_taxa)))
}

unique_observations

endpoint <- max(n_samples_df$n_samples)
res <- iNEXT(unique_observations, datatype = "incidence_freq", q = 0, conf = 0.95, endpoint = endpoint)

acc_curve <- res$iNextEst$size_based %>% rename(Region = Assemblage) 

acculumated_richness <- acc_curve %>%
  group_by(Region) %>%
  filter(abs(t - 50) == min(abs(t - 50))) %>% slice(1)
acculumated_richness$Region <- factor(acculumated_richness$Region, levels = c("FVG", "VEN", "EMR", "MAR", "ABR", "MOL", "PUG", "BAS", "CAL", "SIC", "CAM","LAZ", "TOS", "LIG", "SAR"))
custom_palette <- rep("dodgerblue3", length(unique(acculumated_richness$Region)))
names(custom_palette) <- unique(acculumated_richness$Region)
p <- ggplot(acculumated_richness) + 
    geom_bar(stat = "identity", aes(x = Region, y = qD, fill = Region), colour = "black", linewidth = 1) +
    geom_errorbar(aes(x = Region, ymin = qD.LCL, ymax = qD.UCL), linewidth = 0.5, width = 0.5)  +
    theme_bw(base_size = 18) +
    scale_fill_manual(values = custom_palette) +
    labs(x = "Region", y = paste("Estimated richness", sep = " ")) +
    theme(legend.position = "none") +
    theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) 
p
