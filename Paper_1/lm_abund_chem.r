library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(corrplot)
library(tibble)
library(openxlsx)

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



sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(
    sample_abund = as.integer(sum(Num_cell_l)), 
    Region = first(Region), 
    Season = first(Season), 
    Basin = first(Basin),
    Closest_coast = first(Closest_coast),
    SeaDepth = first(SeaDepth)
    ) 
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

#Open file excel by reading all sheets
sheets <- getSheetNames("./MEMs_per_basin.xlsx")
mems <- sapply(
    sheets,
    function(sheet) {
    data <- read.xlsx("./MEMs_per_basin.xlsx", sheet = sheet)
    }, 
    simplify = FALSE
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

chem_phys %>% dplyr::select(NO2, NH4, NO3) %>% ggplot() + 
geom_point(aes(x = NO2 + NO3 + NH4, y = NO3 + NH4)) + 
xlim(c(0,20)) + 
ylim(c(0,20)) 

chem_phys <- chem_phys %>% mutate(NO_rat = NO2 / NO3, 
                    DIN_TN = (NH4 + NO3) / TN,
                    P_rat = PO4 / TP
)
chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% apply(2, function(x) shapiro.test(x[is.finite(x)])$statistic)

vars_to_transform <- c("Chla", "NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "DIN_TN", "pH", "O_sat")
data_fit <-merge(
    chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, NO_rat, P_rat))  %>% 
    dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120) %>% 
    mutate(across(all_of(vars_to_transform), boxcox_transform)), 
    sample_abund, how = "inner", by = c("Date", "id")
    )


cleaned_data <- data_fit %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))


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



fit_reg_basin <- function(group, group_col, vars) {
    if (nrow(cleaned_data) < 3) {
        return(c(NA, NA))
    }
        model <- lm(
            paste("log10(sample_abund) ~", paste(vars, collapse = " + ")), 
            data = cleaned_data %>% dplyr::filter(!!sym(group_col) == group)
        )
        return(
            summary(model)$coefficients[-1, c(1,4)] %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Var") %>% 
            rename(p_value = `Pr(>|t|)`)
        )
    }

data_fit %>% head()
res <- fit_reg_basin("NA", "New_basin", c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "DIN_TN", "P_rat"))
res

estimates <- sapply(
    data_fit %>% pull(New_basin) %>% unique(),
    function(basin) {
        fit_reg_basin(
            group = basin, 
            group_col = "New_basin", 
            vars = c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "DIN_TN", "P_rat")
        )
        },
    simplify = FALSE
    )



p <- do.call(rbind, estimates) %>% mutate(
    Region = rep(names(estimates), each = nrow(estimates[[1]])) 
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



#MEMs addition and Cloasest_coast and SeaDepth using MEM + VIF
library(car)
library(vegan)
vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "O_sat", "Closest_coast")
cleaned_data <- data_fit %>% dplyr::select(all_of(c("New_basin", vars, "sample_abund", "id"))) %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))


mem_models <- sapply(
    names(mems), 
    function(name) {
        expls <- c(
                    colnames(mems[[name]])[-1],
                    vars
                    )
        data_reduced <- merge(
            cleaned_data %>% dplyr::filter(New_basin == name), 
            mems[[name]],
            by = "id"
        ) %>% dplyr::mutate(across(all_of(expls), ~ decostand(., "standardize")))
    model <- lm(
        paste(
            "log10(sample_abund) ~", 
            paste(
                expls, 
                collapse = " + "
                )
                ), 
        data = data_reduced
    )   
    return(model) 
    }, 
    simplify = FALSE
)


selections <- 
sapply(
    names(mem_models), 
    function(name) {
        expls <- c(
                    colnames(mems[[name]])[-1],
                    vars
                    )
        data_reduced <- merge(
            cleaned_data %>% dplyr::filter(New_basin == name), 
            mems[[name]],
            by = "id"
        ) %>% dplyr::mutate(across(all_of(expls), ~ decostand(., "standardize")))
        newvars <- which(vif(mem_models[[name]]) < 10) %>% names()
    return(
        lm(
            paste(
            "log10(sample_abund) ~", 
            paste(
                newvars, 
                collapse = " + "
                )
                ), data = data_reduced)
    )
    }, 
    simplify = FALSE
)

coefs <- sapply(
    names(selections), 
    function(name) {
        return(
            summary(selections[[name]])$coefficients[-1, c(1,4)] %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Var") %>% 
            rename(p_value = `Pr(>|t|)`) %>% 
            mutate(Region = name)
        )
    }, 
    simplify = FALSE
)
do.call(rbind, coefs) %>% dplyr::filter(!grepl("MEM", Var)) %>% 
ggplot() + 
geom_col(aes(x = Region, y = Estimate, fill = ifelse(p_value < 0.05, "<0.05", ">0.05")), position = "dodge") +
facet_wrap(~Var, scales = "free_y")



#USING LASSO 
library(glmnet)
vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "O_sat", "Closest_coast", "DIN_TN")
cleaned_data <- data_fit %>% dplyr::select(all_of(c("New_basin", vars, "sample_abund", "id"))) %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))
lasso_models <- 
sapply(
    names(mem_models), 
    function(name) {
        expls <- c(
                    colnames(mems[[name]])[-1],
                    vars
                    )
        data_reduced <- merge(
            cleaned_data %>% dplyr::filter(New_basin == name), 
            mems[[name]],
            by = "id"
        ) %>% dplyr::mutate(across(all_of(expls), ~ decostand(., "standardize")))
        model <- cv.glmnet(
        x = as.matrix(data_reduced %>% 
        dplyr::select(
            all_of(expls)
            )
            ),
        y = log10(data_reduced$sample_abund),
        alpha = 1,
        nfolds = 40
    )
    model <- glmnet(
        x = as.matrix(data_reduced %>% 
        dplyr::select(
            all_of(expls)
            )
            ),
        y = log10(data_reduced$sample_abund),
        alpha = 1,
        lambda = model$lambda.min
    )
    return(
        model
    )
    }, 
    simplify = FALSE
)

lasso_coefs <- sapply(
    names(lasso_models), 
    function(name) {
        return(
            data.frame(
                as.matrix(coef(lasso_models[[name]]))[-1, ]
            ) %>% dplyr::mutate(Region = name) %>% 
            dplyr::rename(coef = 1) %>% 
            rownames_to_column(var = "Var")
        )
    }, 
    simplify = FALSE
)
selections <- 
sapply(
    names(lasso_models), 
    function(name) {
        expls <- lasso_coefs[[name]] %>% dplyr::filter(coef != 0) %>% pull(Var)
        data_reduced <- merge(
            cleaned_data %>% dplyr::filter(New_basin == name), 
            mems[[name]],
            by = "id"
        ) %>% dplyr::mutate(across(all_of(expls), ~ decostand(., "standardize")))
    return(
        lm(
            paste(
            "log10(sample_abund) ~", 
            paste(
                expls, 
                collapse = " + "
                )
                ), data = data_reduced)
    )
    }, 
    simplify = FALSE
)
lm_coefs <- sapply(
    names(selections), 
    function(name) {
        return(
            summary(selections[[name]])$coefficients[-1, c(1,4)] %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Var") %>% 
            rename(p_value = `Pr(>|t|)`) %>% 
            mutate(Region = name)
        )
    }, 
    simplify = FALSE
)

do.call(rbind, lm_coefs) %>% dplyr::filter(!grepl("MEM", Var)) %>%
ggplot() + 
geom_col(aes(x = Region, y = Estimate, fill = ifelse(p_value < 0.05, "<0.05", ">0.05")), position = "dodge") +
facet_wrap(~Var)


do.call(rbind, lm_coefs) %>% 
ggplot() + 
geom_tile(aes(x = Region, y = Var, fill = ifelse(p_value < 0.05, Estimate, NA))) + 
scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50") +
labs(fill = "Estimate")


chem_phys %>% 
ggplot() + 
geom_point(aes(x = PO4, y = NO3, col = Region)) + 
facet_wrap(~Region, scales = "free")


phyto_abund %>% dplyr::filter(Class == "Bacillariophyceae") %>%
group_by(Date, id) %>% summarise(
    sample_abund = sum(Num_cell_l), 
    Region = first(Region), 
    Season = first(Season), 
    Basin = first(Basin),
    Closest_coast = first(Closest_coast),
    SeaDepth = first(SeaDepth)
) %>% merge(
    chem_phys, 
    by = c("Date", "id")
) %>% ggplot() + 
geom_point(aes(x = log2(SiO4 + 1), y = log10(sample_abund), col = Region.x)) + 
geom_smooth(aes(x = log2(SiO4 + 1), y = log10(sample_abund), col = Region.x), method = "lm")
