library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(lme4)
library(glmmTMB)

model_val <- function(model) {
    op <- par(mfrow = c(2, 2))
    plot(model)
    return(par(op))   
}

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
HOME_ <- "."
phyto_abund <- read.csv(paste(HOME_, "phyto_abund.csv", sep = "/"))
chem_phys <- read.csv("./df_chem_phys.csv")
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))
sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l)), Region = first(Region), Season = first(Season), Basin = first(Basin)) 

sample_abund %>% group_by(Region, Season) %>% summarise(mean = mean(sample_abund), sd = sd(sample_abund)) %>% 
ggplot() + geom_point(aes(x = log10(mean), y = log10(sd), color = Region, shape = Season)) + theme_minimal() + labs(x = "Mean", y = "SD")

chem_phys %>% dplyr::select(-c(Secchi_depth, id, Date, E_cond)) %>% pivot_longer(cols = -Region, names_to = "var", values_to = "value") %>% 
ggplot() +
geom_boxplot(aes(x = Region, y = value))  + facet_wrap(~ var, scale = "free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Variable", y = "Value")

data_fit <-merge(chem_phys %>% dplyr::select(-c(Secchi_depth, E_cond, Region)) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id")) %>% dplyr::filter(TN / TP < 120)

#NB model ##
model_nb <-  MASS::glm.nb(sample_abund ~ Region + Season, data = sample_abund[-c(1554, 966, 1348),])
model_mixed_nb <- MASS::glmmPQL(sample_abund ~ Region, random = ~ 1 | Season, data = sample_abund, family = negative.binomial(theta = 1))
summary(model_mixed_nb)
model_val(model_mixed_nb)
predicted <- data.frame(
    Region = sample_abund[-c(1554, 966, 1348),]$Region,
    Season = sample_abund[-c(1554, 966, 1348),]$Season,
    predicted = predict(model_nb, newdata = sample_abund[-c(1554, 966, 1348),], type = "response"), 
    residuals = residuals(model_nb, type = "pearson")
)

predicted <- data.frame(
    Region = sample_abund$Region,
    Season = sample_abund$Season,
    predicted = predict(model_mixed_nb, newdata = sample_abund, type = "response"), 
    residuals = residuals(model_mixed_nb, type = "pearson")
)
predicted %>% ggplot() + 
geom_boxplot(aes(x = Region, y = log10(residuals+1), group = Region))
predicted %>% ggplot() + 
geom_point(aes(x = log10(predicted+1), y = log10(residuals+1), color = Region))

mean(predicted$residuals)
sd(predicted$residuals)
#### Chla only ####
data_fit <-merge(chem_phys %>% dplyr::select(c(Date, id, Chla)) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
model_nb_chla <- MASS::glm.nb(sample_abund ~ Chla, data = data_fit)
model_nb_chla <- MASS::glmmPQL(sample_abund ~ Chla, random = ~ 1 | Season, data = data_fit, family = negative.binomial(theta = 1))
summary(model_nb_chla)
model_val(model_nb_chla)
# DO, SiO4, pH 


data_fit %>% pull(DO) %>% is.finite() %>% table()
data_fit <-merge(chem_phys %>% dplyr::select(c(Date, id, DO, NO3, PO4, Salinity, -Region)) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
model_nb_w_nps <- MASS::glm.nb(
    sample_abund ~ DO,
    data = data_fit %>% na.omit()
)

summary(model_nb_w_nps)
drop1(model_nb_w_nps)
model_val(model_nb_w_nps)



model_lmer <- lme4::lmer(sample_abund ~ Region + (1| Season), data = sample_abund)
model_nbmx <- lme4::glmer.nb(sample_abund ~ Season + (1| Region), data = sample_abund)


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



ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")
colors <- scales::hue_pal()(length(unique(phyto_abund$Region)))
palette <- setNames(colors, unique(phyto_abund$Region))
data <- merge(chem_phys %>% dplyr::select(Date, id, Salinity) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
data <- merge(chem_phys %>% dplyr::select(Date, id, TN) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
data <- merge(chem_phys %>% dplyr::select(Date, id, TP) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
data <- merge(chem_phys %>% dplyr::select(Date, id, TN, TP) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))
data$Region <- factor(data$Region, levels = unname(from_region_to_abreviation))
data$Season <- factor(data$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
data$Basin <- factor(data$Basin, levels = ordered_basins)
p <- data %>% dplyr::filter(TN / TP < 120) %>% ggplot(aes(x = Region, y = TN / TP, fill = Region)) +
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

vars <- c("Salinity", "TN", "TP")
fit.data <-  merge(chem_phys %>% dplyr::select(Date, id, all_of(vars)) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id"))


