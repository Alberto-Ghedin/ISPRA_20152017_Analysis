library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(corrplot)
library(tibble)

IMAGE_FORMAT <- "svg"
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


regression_plot_region <- function(data, var, log_env = FALSE) {
    var_sym <- sym(var)
    data %>% dplyr::select(Region, New_basin, Season, sample_abund, !!var_sym) %>% na.omit() %>% 
    ggplot(aes(y = log10(sample_abund), x = !!var_sym, col = Region)) + 
    geom_point() + 
    facet_wrap(~New_basin, scales = "free") + 
    geom_smooth(aes(x = !!var_sym, y = log10(sample_abund), col = Region), method = "lm") + 
    labs(title = var)
}

chem_phys <- chem_phys %>% mutate(NO_rat = NO2 / NO3, 
                    DIN_TN = (NH4 + NO3 + NO2) / TN,
                    P_rat = PO4 / TP
)
chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% apply(2, function(x) shapiro.test(x[is.finite(x)])$statistic)

vars_to_transform <- c("Chla", "NH4", "NO2", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "NO_rat", "DIN_TN", "P_rat", "pH", "O_sat")
data_fit <-merge(
    chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth))  %>% 
    dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120) %>% 
    mutate(across(all_of(vars_to_transform), boxcox_transform)), 
    sample_abund, how = "inner", by = c("Date", "id")
    ) %>% dplyr::select(-c(Date, id))


chem_phys %>% mutate(No2_cont = NO2 / (NO2 + NO3 + NH4)) %>% pull(No2_cont) %>% na.omit() %>% hist(breaks = 20)

chem_phys %>% dplyr::select(NO2, NH4, NO3) %>% ggplot() + 
geom_point(aes(x = NO2 + NO3 + NH4, y = NO2 + NH4)) 

regression_plot_region(data_fit, "Chla")
regression_plot_region(data_fit, "DO")
regression_plot_region(data_fit, "O_sat")
regression_plot_region(data_fit, "NH4", log_env = TRUE)
regression_plot_region(data_fit, "NO3", log_env = FALSE)
regression_plot_region(data_fit, "NO2")
regression_plot_region(data_fit, "PO4", log_env = FALSE)
regression_plot_region(data_fit, "SiO4", log_env = TRUE)
regression_plot_region(data_fit, "Salinity")
regression_plot_region(data_fit, "TN", log_env = FALSE)
regression_plot_region(data_fit, "TP", log_env = TRUE)
regression_plot_region(data_fit, "pH", log_env = FALSE)
regression_plot_region(data_fit, "T")
regression_plot_region(data_fit, "NO_rat", log_env = FALSE)
regression_plot_region(data_fit, "DIN_TN", log_env = FALSE)
regression_plot_region(data_fit, "P_rat", log_env = FALSE)

fit_reg_basin <- function(group, group_col, var) {
    cleaned_data <- data_fit %>% dplyr::filter(!!sym(group_col) == group) %>% 
            dplyr::select(!!sym(group_col), !!sym(var), sample_abund) %>%
            na.omit() %>% dplyr::filter(is.finite(!!sym(var)))
    if (nrow(cleaned_data) < 3) {
        return(c(NA, NA))
    }
        model <- lm(
            paste("log10(sample_abund) ~", var), 
            data = cleaned_data
        )
        return(summary(model)$coefficients[2, c(1,4)])
    }


estimates <- mapply(
    function(var) {
        mapply(
            fit_reg_basin, 
            group = data_fit %>% pull(New_basin) %>% unique(),
            group_col = "New_basin",
            var = var
        )
    }, 
    c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "NO2", "NO_rat", "DIN_TN", "P_rat"), 
    SIMPLIFY = FALSE
)


p <- merge(
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
) %>%
mutate(
    Region = factor(Region, levels = c("NA", "CA", "SA", "SM", "ST", "NT", "LIG", "SAR", "SIC"))
) %>% 
ggplot() + 
geom_col(aes(x = Region, y = Estimate, fill = ifelse(p_value < 0.05, "<0.05", ">0.05")), position = "dodge") + 
facet_wrap(~Var, scales = "free") + 
scale_fill_manual(values = c("blue", "red")) + 
labs(fill = "Significance")

ggsave(
    filename = "lm_abund_chem_basins",
    device = IMAGE_FORMAT,
    plot = p,
    width = 18,
    height = 10,
    units = "in",
    dpi = 300
)


estimates <- mapply(
    function(var) {
        mapply(
            fit_reg_basin, 
            group = data_fit %>% pull(Region) %>% unique(),
            group_col = "Region",
            var = var
        )
    }, 
    c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "NO2", "NO_rat", "DIN_TN", "P_rat"), 
    SIMPLIFY = FALSE
)

p <- merge(
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
) %>% 
mutate(Region = factor(Region, levels = unname(from_region_to_abreviation)) ) %>%
ggplot() + 
geom_col(aes(x = Region, y = Estimate, fill = ifelse(p_value < 0.05, "<0.05", ">0.05")), position = "dodge") + 
facet_wrap(~Var, scales = "free") + 
scale_fill_manual(values = c("blue", "red")) + 
labs(fill = "Significance")

ggsave(
    filename = "lm_abund_chem_regions",
    device = IMAGE_FORMAT,
    plot = p,
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
)

vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "NO2")
chem_phys %>% dplyr::select(all_of(c("id", "Date", "Region", vars))) %>% na.omit() %>% 
merge(sample_abund %>% dplyr::select(id, Date, New_basin), 
by = c("id", "Date")) %>% 
dplyr::select(-c(id, Date)) %>% 
pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>% 
ggplot(aes(y = Value, x = New_basin)) +
geom_boxplot() +
facet_wrap(~Variable, scales = "free")

chem_phys$NO2 %>% log2 %>% hist()
