library(vegan)
library(tidyverse)
library(ggvegan)

from_region_to_abreviation <- c(
    "Friuli-Venezia-Giulia" = "FVG",
    "Veneto" = "VEN", 
    "Emilia-Romagna" = "EMR",
    "Marche" = "MAR",
    "Abruzzo" = "ABR",
    "Molise" = "MOL",
    "Puglia" = "PUG",
    "Basilicata" = "BAS",
    "Calabria" = "CAL",
    "Sicilia" = "SIC",
    "Campania" = "CAM", 
    "Lazio" = "LAZ",
    "Toscana" = "TOS",
    "Liguria" = "LIG",
    "Sardegna" = "SAR"
)

log_trans <- function(x) {
    eps <- x[x != 0] %>% na.omit() %>% min()
    return(as.numeric(log10(x + eps)))
}

chem_phys <- read.csv("./df_chem_phys.csv")
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))

phyto_abund <- read.csv("./phyto_abund.csv")

phyto_abund %>% head()
chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% is.na() %>% colSums() / 2219

vars <- c("NH4", "NO3", "PO4", "SiO4", "TN", "TP")
env.data <- merge(
    chem_phys %>% dplyr::select(-c(E_cond, Secchi_depth, NO2, Chla)) %>% na.omit(), 
    phyto_abund %>% dplyr::select(c(id, Date, Season, Basin)) %>% distinct(),
    by = c("id", "Date"), 
    all.x = TRUE
    )
pca.env <- vegan::rda(env.data %>% dplyr::select(-c(Region, id, Date, Season, Basin)) %>% as.matrix(), scale = TRUE)

summary(pca.env)
screeplot(pca.env, main = "PCA of chemical and physical data", bstick = TRUE)
autoplot(pca.env)

biplot(pca.env, scaling = 2, cex = 0.8)

ggvegan.species <- fortify(pca.env) %>% dplyr::filter(score == "species")
ggvegan.sites <- fortify(pca.env) %>% dplyr::filter(score == "sites")

ggvegan.sites$Region <- env.data$Region
ggvegan.sites$Region <- factor(ggvegan.sites$Region, levels = unname(from_region_to_abreviation))
ggvegan.sites$Season <- env.data$Season
ggvegan.sites$Season <- factor(ggvegan.sites$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
ggvegan.sites$Basin <- env.data$Basin


ggplot() + 
geom_point(
    data = ggvegan.sites, 
    aes(x = PC1, y = PC2, color = Season), 
    size = 3) + 
geom_segment(
    data = ggvegan.species, 
    aes(x = 0, y = 0, xend = PC1, yend = PC2), 
    arrow = arrow(length = unit(0.1, "inches")), 
    color = "black") +
geom_text(
    data = ggvegan.species, 
    aes(x = PC1, y = PC2), label = ggvegan.species$label, 
    hjust = 0, vjust = 0) 



