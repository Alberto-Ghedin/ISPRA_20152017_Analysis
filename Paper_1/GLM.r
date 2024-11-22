library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(lme4)
library(glmmTMB)

HOME_ <- "."
phyto_abund <- read.csv(paste(HOME_, "phyto_abund.csv", sep = "/"))
chem_phys <- read.csv("./df_chem_phys.csv")
sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l)), Region = first(Region), Season = first(Season), Basin = first(Basin)) 
phyto_abund %>% head()
model_nb <-  MASS::glm.nb(sample_abund ~ Region, data = sample_abund)
model_lmer <- lme4::lmer(sample_abund ~ Region + (1| Season), data = sample_abund)
model_nbmx <- lme4::glmer.nb(sample_abund ~ Season + (1| Region), data = sample_abund)

summary(model_nb)
summary(model_lmer)
summary(model_nbmx)
qqnorm(resid(model_nb))
qqline(resid(model_nb))
model_nb
res <- residuals(model_nb)
fitted <- fitted(model_nb)

plot(log(fitted), res, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, col = "red")

dev.off()
plot(model_lmer)

qqnorm(res)
abline(a = 0, b = 1, col = "red")
taxon_abund <- phyto_abund %>% group_by(Taxon) %>% summarise(mean = mean(Num_cell_l), sd = sd(Num_cell_l))


ggplot(taxon_abund) + geom_point(aes(x = log(mean), y = log(sd))) + theme_minimal() + labs(x = "Mean", y = "SD")

sample_abund

merge(chem_phys %>% dplyr::select(Date, id, Salinity) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id")) %>% 
ggplot() + geom_point(aes(x = Salinity, y = log10(sample_abund))) + theme_minimal() + labs(x = "Salinity", y = "Sample abundance")

model_sal <- MASS::glm.nb(sample_abund ~ Salinity + Season + Region, data = merge(chem_phys %>% dplyr::select(Date, id, Salinity) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id")))

summary(model_sal)
res <- residuals(model_sal)
fitted <- fitted(model_sal)
qqnorm(res)
qqline(res)


plot(log10(fitted), res, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, col = "red")

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

ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")
colors <- scales::hue_pal()(length(unique(phyto_abund$Region)))
palette <- setNames(colors, unique(phyto_abund$Region))
data <- merge(chem_phys %>% dplyr::select(Date, id, Salinity) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
data$Region <- factor(data$Region, levels = unname(from_region_to_abreviation))
data$Season <- factor(data$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
data$Basin <- factor(data$Basin, levels = ordered_basins)
p <- data %>% ggplot(aes(x = Region, y = Salinity, fill = Region)) +
    geom_boxplot(width = 0.5, position = position_dodge("preserve")) +
    #scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = palette) +
    facet_grid(Season ~ Basin, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = "none", 
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        panel.spacing = unit(1, "lines")
    ) +
    labs(
        y = "Abundance [cells/L]",
        title = "Sample abundance per basin and season"
    ) 

p
