library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(tibble)
library(openxlsx)
library(rjson)


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


compute_TRIX <- function(data) {
    max_chla <- data %>% pull(Chla) %>% max(na.rm = TRUE) 
    min_chla <- data %>% pull(Chla) %>% min(na.rm = TRUE)
    max_din <- data %>% mutate(DIN = NO3 +  NH4 ) %>% pull(DIN) %>% max(na.rm = TRUE)
    min_din <- data %>% mutate(DIN = NO3 +  NH4 ) %>% pull(DIN)
    min_din <- min_din[min_din != 0] %>% min(na.rm = TRUE)
    max_tp <- data %>% pull(TP[TP !=0]) %>% max(na.rm = TRUE) 
    min_tp <- data %>% pull(TP[TP !=0]) 
    min_tp <- min_tp[min_tp != 0] %>% min(na.rm = TRUE) 
    data %>% pull(TP[TP !=0])
    max_o_dev <- abs(data %>% pull(O_sat) %>% max(na.rm = TRUE) - 100)
    min_o_dev <- abs(data %>% pull(O_sat)  - 100)
    min_o_dev <- min_o_dev[min_o_dev != 0] %>% min(na.rm = TRUE)
 output <- data %>% 
 mutate(
    Chla = ifelse(Chla == 0, min_chla, Chla),
    DIN = ifelse(NO3 + NH4 == 0, min_din, NO3 + NH4),
    TP = ifelse(TP == 0, min_tp, TP),
    O_dev = ifelse(abs(O_sat - min_o_dev) == 0, min_o_dev, abs(O_sat - min_o_dev))
    ) %>% 
mutate(
    TRIX = 10 / 4 * (
        (log10(Chla) - log10(min_chla)) / (log10(max_chla) - log10(min_chla)) +
        (log10(DIN) - log10(min_din)) / (log10(max_din) - log10(min_din)) +
        (log10(TP) - log10(min_tp)) / (log10(max_tp) - log10(min_tp)) + 
        (log10(O_dev) - log10(min_o_dev)) / (log10(max_o_dev) - log10(min_o_dev))
    )
    )
    return(output)
}

plot_variable_along_coast <- function(data, var, group, title, ylab, ordered_latitude = ordered_latitude, ordered_longitude = ordered_longitude) {
    colors <- scales::hue_pal()(length(unique(data[[group]])))
    palette <- setNames(colors, sort(as.character(data[[group]] %>% unique())))
    p <- data %>% 
    mutate(index_id = match(id, params$ordered_id)) %>% 
    ggplot(aes(x = index_id, y = !!as.symbol(var))) + 
    geom_boxplot(aes(group = id, fill = !!as.symbol(group))) +
    scale_fill_manual(values = palette) + 
    labs(title = title, y = ylab) +
    scale_x_continuous(
        name = "Latitude",
        breaks =  seq(1, length(params$ordered_id), 3), 
        labels = one_every_n_item(ordered_latitude, 3),
        limits = c(-0.01, 162.01), 
        sec.axis = sec_axis(
          ~., 
          name = "Longitude",
          breaks = seq(1, length(params$ordered_id), 3),
          labels = one_every_n_item(ordered_longitude, 3)
        )
    ) + 
    ggplot2::theme(
            axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 17),
            axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 17),
            axis.text.y = element_text(angle = 0, hjust = 0, size = 17),
            axis.title.y = element_text(size = 22),
            axis.title.x = element_text(size = 22),
            plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18, face = "bold"),
        ) 
    return(p)
} 

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
HOME_ <- "./Paper_1"
phyto_abund <- read.csv(paste(HOME_, "phyto_abund.csv", sep = "/")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% mutate(
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
chem_phys <- read.csv(paste(HOME_, "df_chem_phys.csv", sep = "/"))
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))


params <- fromJSON(file = file.path(HOME_, "params.json"))


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
sheets <- getSheetNames(paste(HOME_, "MEMs_per_basin.xlsx", sep = "/"))
mems <- sapply(
    sheets,
    function(sheet) {
    data <- read.xlsx(paste(HOME_, "MEMs_per_basin.xlsx", sep = "/"), sheet = sheet)
    }, 
    simplify = FALSE
    )


one_every_n_item <- function(list, n) {
    return(
        list[seq(1, length(list), n)]
    )
}
ordered_latitude <- phyto_abund %>% dplyr::select(id, Latitude) %>% distinct() %>% 
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Latitude) %>% round(2)
ordered_longitude <- phyto_abund %>% dplyr::select(id, Longitude) %>% distinct() %>%
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Longitude) %>% round(2)
colors <- scales::hue_pal()(length(unique(phyto_abund$New_basin)))
palette <- setNames(colors, sort(as.character(phyto_abund$New_basin %>% unique())))
colors <- scales::hue_pal()(length(unique(phyto_abund$Region)))
palette <- setNames(colors, sort(as.character(phyto_abund$Region %>% unique())))


station_sufficient_samples <- phyto_abund %>% group_by(id) %>% summarise(
    n_samples = n_distinct(Date)
) %>% arrange(desc(n_samples)) %>% dplyr::filter(n_samples >= 10) %>% pull(id)

p <- chem_phys %>% merge(phyto_abund %>% dplyr::select(Date, id, New_basin), by = c("Date", "id")) %>% 
mutate(index_id = match(id, params$ordered_id)) %>% 
ggplot(aes(x = index_id, y = log10(Chla))) + 
geom_boxplot(aes(group = id, fill = Region)) +
scale_fill_manual(values = palette) + 
labs(title = "Sample Chla across all stations", y = "Chla [mumol/L] (log scale)") +
scale_x_continuous(
    name = "Latitude",
    breaks =  seq(1, length(params$ordered_id), 3), 
    labels = one_every_n_item(ordered_latitude, 3),
    limits = c(-0.01, 162.01), 
    sec.axis = sec_axis(
      ~., 
      name = "Longitude",
      breaks = seq(1, length(params$ordered_id), 3),
      labels = one_every_n_item(ordered_longitude, 3)
    )
) + 
ggplot2::theme(
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 17),
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 17),
        axis.text.y = element_text(angle = 0, hjust = 0, size = 17),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18, face = "bold"),
    ) 
ggsave(
    file.path(HOME_, "chla_per_region.svg"), 
    p, 
    width = 18, 
    height = 8.5, 
    dpi = 300
)



phyto_abund %>% group_by(id) %>% summarise(
    n_samples = n_distinct(Date), 
    Region = first(Region)
) %>% arrange(desc(n_samples)) %>% dplyr::filter(n_samples < 10) 

chem_phys <- read.csv(paste(HOME_, "df_chem_phys.csv", sep = "/"))
chem_phys <- compute_TRIX(chem_phys)
chem_phys <- chem_phys %>% mutate(
                   #  NO_rat = NO2 / NO3, 
                   #DIN_TN = (NH4 + NO3) / TN,
                   # P_rat = PO4 / TP, 
                    NP = NO3 / PO4, 
                    NP_tot = TN / TP
)

vars_to_transform <- c("NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "NP_tot", "pH")

chem_phys <- chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla)) %>%
dplyr::filter(
  DO < 400, 
  NH4 < 10,
  NO3 < 57, 
  O_sat < 160, 
  pH > 7, 
  PO4 < 2, 
  Salinity > 20, 
  SiO4 < 40, 
  TN < 120, 
  NP_tot < 204
  ) #%>% 
merge(
  phyto_abund %>% dplyr::distinct(id, Date, Closest_coast, SeaDepth, Season, New_basin), 
  by = c("Date", "id")
) 

chem_phys <- chem_phys %>% mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
dplyr::select(-c(O_sat, NP)) %>%
na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) 

ordered_transect <- c(
    "SMTS", "SMLG", "VENEZIA", "ROSOLINA", "PORTO_GARIBALDI", "CESENATICO", "RIMINI", "Chienti", "Esino", "GU",
    "VA", "R14001_B2", "FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO", "PORTO_CESAREO", "PUNTA_RONDINELLA", "SINNI", "Villapiana",
    "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche", "Isole_Ciclopi", "Plemmirio", "Isola_Correnti", "San_Marco", "Isole_Egadi", "Capo_Gallo","Vibo_marina", 
    "Cetraro", "Cilento", "Salerno", "Napoli", "Domizio", "m1lt01", "m1lt02",  "m1rm03",  "m1vt04", "Collelungo","Carbonifera",
    "Donoratico", "Fiume_Morto", "Mesco" , "Portofino"  ,  "Voltri"  , "Quiliano" , 
    "Olbia", "Arbatax", "Villasimius", "Cagliari", "Oristano", "Alghero", "Porto_Torres"

)
sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))

sample_abund <- merge(
    sample_abund, 
    sea_depth %>% dplyr::select(id, SeaDepth, Transect)
)
sample_abund$Transect <- factor(sample_abund$Transect, levels = ordered_transect, ordered = TRUE)

ids <- phyto_abund %>% dplyr::filter(Region == "CAL") %>% pull(id) %>% unique()
chem_phys %>% dplyr::filter(id %in% ids) %>% 
dplyr::select(-c(O_dev)) %>%
pivot_longer(
  cols = c("NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "TRIX", "NP_tot", "DO", "DIN"),
  names_to = "Variable",
  values_to = "Value"
) %>% 
merge(
    sample_abund %>% dplyr::select(Date, id, Transect)
) %>% 
mutate(id = factor(id, levels = params$ordered_id[params$ordered_id %in% ids], ordered = TRUE)) %>%
ggplot() + 
geom_boxplot(aes(x = id, y = Value, fill = Transect)) + 
facet_wrap(~ Variable, scales = "free")


chem_phys %>% dplyr::filter(id %in% cam_ids) %>% 
dplyr::select(all_of(c("id", "Date", "NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "TRIX", "NP_tot", "DO", "DIN", "SeaDepth"))) %>% 
merge(
    sample_abund %>% dplyr::select(Date, id, sample_abund, Transect), 
    by = c("Date", "id")
) %>% 
pivot_longer(
    cols = -all_of(c("Date", "id", "sample_abund", "Transect")),
    names_to = "Variable",
    values_to = "Value"
) %>% ggplot() + 
geom_point(aes(x = Value, y = log10(sample_abund), col = Transect)) +
facet_wrap(~ Variable, scales = "free")



cleaned_data <- data_fit %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))



chem_phys %>% dplyr::filter(Region %in% c("FVG", "VEN", "EMR", "MOL", "PUG", "ABR", "MAR")) %>% 
dplyr::select(Date, id, Region, Salinity) %>% 
merge(
    sample_abund %>% dplyr::select(Date, id, Season),
    by = c("Date", "id")
) %>% 
mutate(year = format(as.Date(Date), "%Y")) %>% 
ggplot() + 
geom_boxplot(aes(x = year, y = Salinity, fill = year)) + 
facet_wrap(~Season, scales = "free")




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


cleaned_data %>% colnames()
cleaned_data <- cleaned_data %>% mutate(NA_div = 
case_when(
    New_basin == "NA" ~ "NA",
    New_basin != "NA" ~ "ELSE",
    ))


estimates <- sapply(
    cleaned_data %>% pull(NA_div) %>% unique(),
    function(basin) {
        fit_reg_basin(
            group = basin, 
            group_col = "NA_div", 
            vars = c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "DIN_TN", "O_sat", "Closest_coast")
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
p

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

#Using region, no MEM 
library(car)
vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "Closest_coast")
cleaned_data <- data_fit %>% dplyr::select(all_of(c("New_basin", vars, "sample_abund", "id"))) %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))


#MEMs addition and Cloasest_coast and SeaDepth using MEM + VIF
library(car)
library(vegan)
vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "O_sat", "Closest_coast")
cleaned_data <- data_fit %>% dplyr::select(all_of(c("Region", "New_basin", vars, "sample_abund", "id"))) %>% 
    na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))


full_models <- sapply(
    cleaned_data %>% pull(New_basin) %>% unique(), 
    function(name) {
        data_reduced <- cleaned_data %>% 
        dplyr::filter(New_basin == name) %>% 
        dplyr::mutate(across(all_of(vars), ~ decostand(., "standardize")))
        model <- lm(
            paste(
                "log10(sample_abund) ~", 
                paste(
                    c("id", 
                    vars
                    ),
                    collapse = " + "
                )
            ), 
            data = data_reduced
        )
        }, 
    simplify = FALSE
)
data_reduced
selections <- sapply(
    names(full_models), 
    function(name) {
        data_reduced <- cleaned_data %>% 
        dplyr::filter(New_basin == name) %>% 
        dplyr::mutate(across(all_of(vars), ~ decostand(., "standardize")))
        new_model <- step(full_models[[name]])
    }
)

selections[[1]]

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
    cleaned_data %>% dplyr::select(Date, id, SiO4),
    by = c("Date", "id")
) %>% ggplot() + 
geom_point(aes(x =SiO4, y = log10(sample_abund), col = Region))# + 
geom_smooth(aes(x = log2(SiO4 + 1), y = log10(sample_abund), col = Region), method = "lm")

chem_phys %>% head()
chem_phys %>% 
ggplot() + 
geom_histogram(aes(x = PO4 / TP))
quantile(chem_phys %>% mutate(P_rat = PO4 / TP) %>% pull(P_rat),probs = seq(0, 1, 0.25), na.rm = TRUE)
#TRIX index 
chem_phys %>% mutate(
    chla_gm3 = Chla * 0.8935, 
    DIN_gm3 = NO3 * 0.062 +  NH4 * 0.018,
    TP_gm3 = TP * 0.031,
    O_dev = abs(O_sat - 100),
) %>% 
mutate(
    TRIX = 10 / 4 * (
        (log10(chla_gm3) - log10(min(chla_gm3, na.rm = TRUE))) / (log10(max(chla_gm3, na.rm = TRUE)) - log10(min(chla_gm3, na.rm = TRUE))) +
        (log10(DIN_gm3) - log10(min(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE))) / (log10(max(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE)) - log10(min(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE))) +
        (log10(TP_gm3) - log10(min(TP_gm3[TP_gm3 != 0], na.rm = TRUE))) / (log10(max(TP_gm3[TP_gm3 != 0], na.rm = TRUE)) - log10(min(TP_gm3[TP_gm3 != 0], na.rm = TRUE))) + 
        (log10(O_dev) - log10(min(O_dev[O_dev != 0], na.rm = TRUE))) / (log10(max(O_dev, na.rm = TRUE)) - log10(min(O_dev[O_dev != 0], na.rm = TRUE)))
    )
) %>% merge(
    sample_abund, 
    by = c("Date", "id")
) %>% ggplot() + 
geom_point(aes(x = TRIX, y = log10(sample_abund), col = New_basin)) + 
geom_smooth(aes(x = TRIX, y = log10(sample_abund), col = New_basin), method = "lm") + 
facet_wrap(~New_basin)


chem_phys %>% mutate(
    chla_gm3 = Chla * 0.8935, 
    DIN_gm3 = NO3 * 0.062 +  NH4 * 0.018,
    TP_gm3 = TP * 0.031,
    O_dev = abs(O_sat - 100),
) %>% 
mutate(
    TRIX = 10 / 4 * (
        (log10(chla_gm3) - log10(min(chla_gm3, na.rm = TRUE))) / (log10(max(chla_gm3, na.rm = TRUE)) - log10(min(chla_gm3, na.rm = TRUE))) +
        (log10(DIN_gm3) - log10(min(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE))) / (log10(max(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE)) - log10(min(DIN_gm3[DIN_gm3 != 0], na.rm = TRUE))) +
        (log10(TP_gm3) - log10(min(TP_gm3[TP_gm3 != 0], na.rm = TRUE))) / (log10(max(TP_gm3[TP_gm3 != 0], na.rm = TRUE)) - log10(min(TP_gm3[TP_gm3 != 0], na.rm = TRUE))) + 
        (log10(O_dev) - log10(min(O_dev[O_dev != 0], na.rm = TRUE))) / (log10(max(O_dev, na.rm = TRUE)) - log10(min(O_dev[O_dev != 0], na.rm = TRUE)))
    )
) %>% merge(
    sample_abund, 
    by = c("Date", "id")
) %>% merge(
    abund_groups %>% pivot_wider(
    names_from = Group, 
    values_from = Abund, 
    values_fill = 0
) %>% mutate(
    DIA_DIN = DIA / DIN, 
    COC_DIN = COC / DIN,
    CRY_DIN = CRY / DIN
) %>% dplyr::select(Date, id, DIA_DIN, COC_DIN, CRY_DIN), 

by = c("Date", "id")
) %>% ggplot() + 
geom_point(aes(x = TRIX, y = log10(DIA_DIN), col = Season)) + 
facet_wrap(~New_basin) 

merge(
    abund_groups %>% pivot_wider(
    names_from = Group, 
    values_from = Abund, 
    values_fill = 0
) %>% mutate(
    DIA_DIN = DIA / DIN
) %>% dplyr::select(Date, id, DIA_DIN, Season), 
    phyto_abund %>% dplyr::select(Latitude, Longitude, id) %>% distinct(), 
    by = "id"
) %>% 
ggplot() + 
geom_point(aes(x = Longitude, y = Latitude, col = log10(DIA_DIN), size = log10(DIA_DIN) + 3)) + 
facet_wrap(~Season, scales = "free") 



chem_phys <- read.csv(paste(HOME_, "df_chem_phys.csv", sep = "/"))
chem_phys <- compute_TRIX(chem_phys)
chem_phys <- chem_phys %>% mutate(
                   #  NO_rat = NO2 / NO3, 
                   #DIN_TN = (NH4 + NO3) / TN,
                   # P_rat = PO4 / TP, 
                    NP = NO3 / PO4, 
                    NP_tot = TN / TP
)


chem_phys <- chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla)) %>%
dplyr::filter(
  DO < 400, 
  NH4 < 10,
  NO3 < 57, 
  O_sat < 160, 
  pH > 7, 
  PO4 < 2, 
  Salinity > 20, 
  SiO4 < 40, 
  TN < 120, 
  NP_tot < 204
  )

data_test <- chem_phys %>%
merge(
    sample_abund %>% dplyr::select(Date, id, sample_abund, Season, New_basin, Region),
    by = c("Date", "id")
) %>% dplyr::filter(
    Region %in% c("FVG", "VEN", "EMR")
)

data_test %>% 
ggplot() + 
geom_point(aes(x = TRIX, y = log10(sample_abund))) + 
facet_wrap(~Region, scales = "free") 
dev.off()
