library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(lme4)
library(glmmTMB)
library(corrplot)
library(car)
library(tibble)

model_val <- function(model) {
    op <- par(mfrow = c(2, 2))
    plot(model)
    return(par(op))   
}


log_trans <- function(x) {
    eps <- x[x != 0] %>% na.omit() %>% min()
    return(as.numeric(log10(x + eps)))
}

boxcox_transform <- function(values) {
  
  if (any(is.na(values))) {
    clean_values <- as.numeric(na.omit(values))
  } else {
     clean_values <- as.numeric(values)
  }

  if (any(is.infinite(clean_values))) {
    clean_values <- as.numeric(clean_values[is.finite(clean_values)])
  }


  if (min(clean_values) == 0.0) {
    min_val <- clean_values[clean_values != 0] %>% min()
    clean_values <- clean_values + min_val * 0.1
  }
  lambda <- MASS::boxcox(clean_values ~ 1, plotit = FALSE)
  max_lambda <- lambda$x[which(lambda$y == max(lambda$y))]
  if (max_lambda == 0) {
    clean_values <- log(clean_values)
  } else {
    clean_values <- (clean_values^max_lambda - 1) / max_lambda
  }
  values[which(!(is.na(values) | is.infinite(values)))] <- clean_values
  return(values)
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
phyto_abund <- read.csv("./phyto_abund.csv") %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% mutate(
    New_basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Region == "PUG" & Basin == "SouthAdr" ~ "SA",
        Region == "PUG" & Basin == "Ion" ~ "SM",
        Region == "BAS" ~ "SM",
        Region == "CAL" & Basin == "Ion" ~ "SM",
        Region == "SIC" ~ "SIC", 
        Region == "CAL" & Basin == "SouthTyr" ~ "ST",
        Region == "CAM" ~ "ST",
        Region == "LAZ" & Basin == "SouthTyr" ~ "ST",
        Region == "LAZ" & Basin == "NorthTyr" ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
    )
chem_phys <- read.csv("./df_chem_phys.csv")
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))


sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l)), Region = first(Region), Season = first(Season), Basin = first(Basin), Closest_coast = first(Closest_coast)) 
sample_abund <- sample_abund %>% mutate(
    New_basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Region == "PUG" & Basin == "SouthAdr" ~ "SA",
        Region == "PUG" & Basin == "Ion" ~ "SM",
        Region == "BAS" ~ "SM",
        Region == "CAL" & Basin == "Ion" ~ "SM",
        Region == "SIC" ~ "SIC", 
        Region == "CAL" & Basin == "SouthTyr" ~ "ST",
        Region == "CAM" ~ "ST",
        Region == "LAZ" & Basin == "SouthTyr" ~ "ST",
        Region == "LAZ" & Basin == "NorthTyr" ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
)


chem_phys <- chem_phys %>% mutate(NO_rat = NO2 / NO3, 
                    DIN_TN = (NH4 + NO3) / TN,
                    P_rat = PO4 / TP, 
                    N_star = NO3 - 16 * PO4
)


png("./chem_phys.png", width = 800, height = 800)
chem_phys %>% dplyr::select(DO, O_sat, NH4, NO3, TN, PO4, TP, SiO4, Salinity, pH, NO_rat, DIN_TN, P_rat, N_star) %>% 
pairs()
dev.off()



summary(lm(chem_phys$O_sat ~ chem_phys$DO))
vars_to_transform <- c("Chla", "NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "DIN_TN", "pH", "O_sat", "P_rat")
data_fit <-merge(
    chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, NO_rat))  %>% 
    dplyr::filter(pH > 7) %>% 
    mutate(across(all_of(vars_to_transform), boxcox_transform)), 
    sample_abund, how = "inner", by = c("Date", "id")
    )  %>% 
    dplyr::filter(Region != "BAS")


vars <- c("T", "O_sat", "TP", "P_rat", "Salinity", "Closest_coast", "N_star") #"NO3", "PO4", "DO"
cleaned_data <- data_fit %>% 
    na.omit() %>% 
    dplyr::select(all_of(c("Region", "Season", "New_basin", vars, "sample_abund"))) %>% 
    dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) %>% 
    dplyr::mutate_at(all_of(vars), ~ (scale(.) %>% as.vector())) 
  
cleaned_data$Region <- factor(cleaned_data$Region, levels = unname(from_region_to_abreviation))
cleaned_data$Season <- factor(cleaned_data$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
cleaned_data$New_basin <- factor(cleaned_data$New_basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"))


## Testing lmer4 ##
### Probably Regional effect, can be estimated with a mixed term 
### Season need not be included in the model
### It is better to use full fixed model (AIC)

model_fixed <- lm(
    paste("log10(sample_abund) ~ ",  paste(vars, collapse = " + "), "+ Region"), 
    data = cleaned_data
    )    
summary(model_fixed)
car::vif(model_fixed)
step(model_fixed)

model_fixed <- lm(
    paste("log10(sample_abund) ~  Region + T + O_sat + TP + P_rat + Salinity + Closest_coast + N_star"), 
    data = cleaned_data
    )
summary(model_fixed)

cleaned_data %>% dplyr::select(P_rat, N_star) %>% cor()
cleaned_data %>%
ggplot() + 
geom_point(aes(y = log10(sample_abund), x = DO, col = Region)) + 
geom_smooth(aes(y = log10(sample_abund), x = DO), method = "lm") +
facet_wrap(~New_basin)


evaluation <- data.frame(
    fitted = fitted(model_fixed), 
    residuals = resid(model_fixed)
) %>% cbind(cleaned_data)

evaluation %>% 
ggplot() + 
geom_point(aes(Season, residuals))

evaluation %>% pivot_longer(
    cols = -c(Region, Season, New_basin, sample_abund, fitted, residuals), 
    names_to = "var",
    values_to = "value"
) %>% ggplot() + 
    geom_point(aes(x = value, y = log10(sample_abund), color = Region)) + 
    geom_smooth(aes(x = value, y = log10(sample_abund)), method = "lm") +
    facet_wrap(~var, scales = "free")

summary(model_fixed)$coefficients




##GLM ##
m1 <- glmmTMB(
    sample_abund ~ (1 | Region) + Closest_coast + T , 
    data = cleaned_data, 
    family = nbinom2
)
summary(m1)
sigma(m1)
qqnorm(resid(m1, type = "pearson"))
qqline(resid(m1, type = "pearson"))

resid(m1) %>% min()
data.frame(
    res = resid(m1, type = "pearson"), 
    fitted = fitted(m1)
) %>% cbind(cleaned_data)  %>% 
ggplot() +
geom_boxplot(aes(x = Region, y = res, group = Region)) +
#geom_point(aes(x = log10(fitted), y = res, col = Region)) +
scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))
stat_ellipse(aes(x = DO, y = fitted, col = Basin))
geom_smooth(aes(x = DO, y = res, col = Basin), method = "lm") 


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



library(vegan)


pca_NA <- rda(
    chem_phys %>% 
    dplyr::filter(Region %in% c("FVG", "VEN", "EMR")) %>%
    dplyr::select(DO, NH4, NO3, O_sat, PO4, Salinity, SiO4, T, TN, TP, pH) %>% mutate(across(where(is.numeric), ~ decostand(., "standardize"))) %>% na.omit(), 
    tidy = TRUE
    )
summary(pca_NA)
sites <- cbind(
    scores(pca_NA, display = "sites"), 
    chem_phys %>% 
    dplyr::filter(Region %in% c("FVG", "VEN", "EMR")) %>%
    dplyr::select(Region, DO, NH4, NO3, O_sat, PO4, Salinity, SiO4, T, TN, TP, pH) %>% mutate(across(where(is.numeric), ~ decostand(., "standardize"))) %>% na.omit()
    )
env_arrows <- scores(pca_NA, display = "species")

sites %>% 
ggplot() + 
geom_point(data = sites, aes(x = PC1, y = PC2, col = Region)) +
geom_segment(data = env_arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), col = "black") + 
geom_text(data = env_arrows, aes(x = PC1, y = PC2, label = rownames(env_arrows)), size = 3, hjust = 0.5, vjust = 0.5)

pca_not_NA <- rda(
    chem_phys %>% 
    dplyr::filter(!Region %in% c("FVG", "VEN", "EMR", "CAM")) %>%
    dplyr::select(DO, NH4, NO3, O_sat, PO4, Salinity, SiO4, T, TN, TP, pH) %>% mutate(across(where(is.numeric), ~ decostand(., "standardize"))) %>% na.omit()
    )
summary(pca_not_NA)
sites <- cbind(
    scores(pca_not_NA, display = "sites"), 
    chem_phys %>% 
    dplyr::filter(!Region %in% c("FVG", "VEN", "EMR", "CAM")) %>%
    dplyr::select(Region, DO, NH4, NO3, O_sat, PO4, Salinity, SiO4, T, TN, TP, pH) %>% mutate(across(where(is.numeric), ~ decostand(., "standardize"))) %>% na.omit()
    )
env_arrows <- scores(pca_not_NA, display = "species")
sites %>% 
ggplot() + 
geom_point(data = sites, aes(x = PC1, y = PC2, col = Region)) +
geom_segment(data = env_arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), col = "black") + 
geom_text(data = env_arrows, aes(x = PC1, y = PC2, label = rownames(env_arrows)), size = 3, hjust = 0.5, vjust = 0.5)


chem_phys %>% dplyr::select(-c(E_cond, Secchi_depth, NO2, DIN_TN, NO_rat, P_rat)) %>% 
pivot_longer(cols = -c(Date, id, Region), names_to = "var", values_to = "value") %>% 
ggplot() + 
geom_boxplot(aes(x = Region, y = value)) + 
facet_wrap(~var, scales = "free") 
