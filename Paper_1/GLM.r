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
sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l)), Region = first(Region), Season = first(Season), Basin = first(Basin)) 
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
  values[which(!is.na(values))] <- clean_values
  return(values)
}

vars <- c("DO", "NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "T", "pH")
vars_to_transform <- c("NH4", "NO3", "pH", "PO4", "Salinity", "SiO4", "TN", "TP", "T")
chem_phys %>% dplyr::select(-c(Secchi_depth, id, Date, E_cond, Chla, NO2, O_sat)) %>% 
mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
pivot_longer(cols = -Region, names_to = "var", values_to = "value") %>% 
ggplot() +
geom_histogram(aes(x = value))  + facet_wrap(~ var, scale = "free") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Variable", y = "Value")

data_fit <-merge(
    chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla, O_sat))  %>% 
    dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120) %>% 
    mutate(across(all_of(vars_to_transform), boxcox_transform)), 
    sample_abund, how = "inner", by = c("Date", "id")
    ) %>% dplyr::select(-c(Date, id))
    






abund_only_genera <- phyto_abund %>% mutate(
    Det_level = case_when(
        Class == "nan" ~ Taxon,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
) %>% group_by(Date, id, Det_level) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    basin = first(New_basin),
    Season = first(Season),
    .groups = "drop"
) %>% pivot_wider(names_from = Det_level, values_from = Abund, values_fill = 0)




merge(
    sapply(
    names(estimates),
    function(name) {
    return(
estimates[[name]][1, ])
}
) %>% as.data.frame() %>% rownames_to_column(var = "Region") %>% pivot_longer(cols = -Region, names_to = "Var", values_to = "Estimate"), 
sapply(
    names(estimates),
    function(name) {
    return(
estimates[[name]][2, ])
}
) %>% as.data.frame() %>% rownames_to_column(var = "Region") %>% pivot_longer(cols = -Region, names_to = "Var", values_to = "p_value"), 
by = c("Region", "Var")
) %>% ggplot() + 
geom_col(aes(x = Region, y = Estimate, fill = ifelse(p_value < 0.05, "blue", "red")), position = "dodge") + 
facet_wrap(~Var, scales = "free") + 
scale_fill_manual(values = c("blue", "red"))






season <- "Winter"
correlation <- chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% 
mutate(across(all_of(vars), log_trans)) %>% na.omit() %>% cor()


corrplot.mixed(correlation, order = 'AOE')



lm_statistic <- function(data, var, spat_col = "Region") {
    region <- data[[spat_col]][1]
    season <- data$Season[1]
    data_cleaned <- data %>% dplyr::select(sample_abund, !!sym(var)) %>% na.omit()
    info <- c(Var = var, Region = region, Season = season)
    if (nrow(data_cleaned) < 7) {
        return(
            data.frame(
                Var = var, Region = region, Season = season,
                Estimate = NA, `Pr(>|t|)` = NA
            )
        )
    }
    model <- lm(formula(paste("log10(sample_abund) ~", var)), data = data_cleaned)
    return(
        data.frame(
            Var = var, Region = region, Season = season,
            summary(model)$coefficients[2, c(1,4)] %>% as.list()
        )
    )
}

variables <- c("DO", "NH4", "NO3","PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T")
model <- lm(formula(paste("log10(sample_abund) ~", paste(variables, collapse = "+"))), data = data_fit %>% na.omit())
summary(model)
car::vif(model)
summary(model)$coefficients[2, c(1,4)] %>% as.list()

paste("log10(sample_abund) ~", vars, collapse = "+")



statistics <- bind_rows(
    lapply(variables, function(x) {lapply(data_fit %>% group_split(Season, Region), function(y) lm_statistic(y, x)) %>% bind_rows()}) 
    )# %>% na.omit()
statistics$Region <- factor(statistics$Region, levels = unname(from_region_to_abreviation))
statistics$Season <- factor(statistics$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))



statistics %>% mutate(Estimate = case_when(
    Pr...t.. > 0.05 ~ NA,
    TRUE ~ Estimate
)) %>% 
ggplot(aes(x = Season, y = Region, fill = Estimate)) +
facet_wrap(~Var, scales = "free") +
geom_tile() + labs(x = "Season", y = "Region", fill = "Estimate") + 
scale_fill_gradient2(low = "blue", high = "red", limits = c(-1, 1))

statistics <- bind_rows(
    lapply(variables, function(x) {lapply(data_fit %>% group_split(Season, Basin), function(y) lm_statistic(y, x, spat_col = "Basin")) %>% bind_rows()}) 
    )
statistics$Region <- factor(statistics$Region, levels = c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed"))
statistics$Season <- factor(statistics$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))

statistics %>% mutate(Estimate = case_when(
    Pr...t.. > 0.05 ~ NA,
    TRUE ~ Estimate
)) %>% 
ggplot(aes(x = Season, y = Region, fill = Estimate)) +
facet_wrap(~Var, scales = "free") +
geom_tile() + labs(x = "Season", y = "Region", fill = "Estimate") + 
scale_fill_gradient2(low = "blue", high = "red", limits = c(-1, 1))

lm(formula(paste("log10(sample_abund) ~", ".")), data = data_fit %>% dplyr::filter(Basin == "NorthAdr", Season == "Summer") %>% dplyr::select(sample_abund, all_of(variables))) %>% summary()

data_fit %>% mutate(N_P = TN / TP) %>% dplyr::filter(N_P < 120) %>% dplyr::filter(!Region %in% c("EMR", "FVG", "VEN")) %>% ggplot() + 
geom_point(aes(y = log10(sample_abund), x = N_P, color = Region)) + theme_minimal() 
data_fit %>% mutate(N_P = TN / TP) %>% dplyr::filter(N_P < 120) %>% dplyr::filter(Region %in% c("EMR", "FVG", "VEN")) %>% ggplot() + 
geom_point(aes(y = log10(sample_abund), x = N_P, color = Region)) + theme_minimal()#+ labs(x = "TN", y = "TP")
data_fit %>% dim()
data_fit <-merge(chem_phys %>% dplyr::select(c(Date, id, -Region, Chla, Salinity, TP, TN)) %>% na.omit(), sample_abund, how = "inner", by = c("Date", "id")) #%>% dplyr::filter(TN / TP < 120)

## Testing lmer4 ##
model_lmer <- lme4::glmer(sample_abund ~Salinity + (1 | Region) + (1| Season), data = data_fit, family = negative.binomial(theta = 1))
summary(model_lmer)
make_prediction_df <- function(data, model) {
    predicted <- data.frame(
        Region = data$Region,
        Season = data$Season,
        predicted = predict(model, newdata = data, type = "response"), 
        residuals = residuals(model, type = "pearson")
    )
    return(predicted)
}
predicted <- make_prediction_df(data_fit, model_lmer)
predicted %>% ggplot() + 
geom_point(aes(x = log10(predicted+1), y = log10(residuals+1), fill = Region), colour = "black", pch = 21, size = 3, alpha = 0.8)
predicted %>% ggplot() + 
geom_boxplot(aes(x = Region, y = log10(residuals+1), group = Region))
model_val(model_lmer)
plot(model_lmer)

drop1(model_lmer)
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


