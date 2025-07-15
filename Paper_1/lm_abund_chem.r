library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(tibble)
library(openxlsx)
library(rjson)



IMAGE_FORMAT <- "svg"
HOME_ <- "./Paper_1"

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))
#Open file excel by reading all sheets
sheets <- getSheetNames(paste(HOME_, "MEMs_per_basin.xlsx", sep = "/"))
mems <- sapply(
    sheets,
    function(sheet) {
    data <- read.xlsx(paste(HOME_, "MEMs_per_basin.xlsx", sep = "/"), sheet = sheet)
    }, 
    simplify = FALSE
    )


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


sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(
    sample_abund = as.integer(sum(Num_cell_l)), 
    Region = first(Region), 
    Season = first(Season), 
    Basin = first(Basin),
    Closest_coast = first(Closest_coast),
    SeaDepth = first(SeaDepth), 
    Longitude = first(Longitude),
    Latitude = first(Latitude), 
    Region = first(Region),
    ) %>% mutate(
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
) %>% merge(
    sample_abund, 
    sea_depth %>% dplyr::select(id, SeaDepth, Transect)
) %>% 
dplyr::mutate(
    Transect = factor(Transect, levels = ordered_transect, ordered = TRUE)
) 

abund_groups <- process_abund_groups(phyto_abund)



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


chem_phys <- read.csv(paste(HOME_, "df_chem_phys.csv", sep = "/"))
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))

chem_phys <- chem_phys %>% mutate(
                    NP = NO3 / PO4, 
                    DIN = NO3 + NH4,
                    NP_tot = TN / TP
)


chem_phys <- chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla)) %>%
dplyr::filter(
  DO < 400 | is.na(DO), 
  NH4 < 10 | is.na(NH4),
  NO3 < 57 | is.na(NO3), 
  O_sat < 160 | is.na(O_sat), 
  pH > 7 | is.na(pH), 
  PO4 < 2 | is.na(PO4), 
  Salinity > 20 | is.na(Salinity), 
  SiO4 < 40 | is.na(SiO4), 
  TN < 120 | is.na(TN), 
  NP_tot < 204 | is.na(NP_tot)
  ) %>% 
chem_phys$TRIX <- compute_TRIX(chem_phys)
merge(
  phyto_abund %>% dplyr::distinct(id, Date, Closest_coast, SeaDepth, Season, New_basin), 
  by = c("Date", "id")
) 

vars_to_transform <- c("NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "NP_tot", "pH")
chem_phys <- chem_phys %>% mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
dplyr::select(-c(O_sat, NP)) %>%
na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) 


ids <- phyto_abund %>% dplyr::filter(Region %in% c("CAL", "SIC", "BAS")) %>% pull(id) %>% unique()
chem_phys %>% dplyr::filter(id %in% ids) %>% 
dplyr::filter(NP_tot < 1000) %>%
pivot_longer(
  cols = c("NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "TRIX", "NP_tot", "DO", "DIN"),
  names_to = "Variable",
  values_to = "Value"
) %>% 
merge(
    sample_abund %>% dplyr::select(id, Transect), 
    by = c("id")
) %>% 
mutate(id = factor(id, levels = params$ordered_id[params$ordered_id %in% ids], ordered = TRUE)) %>%
ggplot() + 
geom_boxplot(aes(x = id, y = Value, fill = Transect)) + 
facet_wrap(~ Variable, scales = "free")


plot_variable_along_coast(
    data = chem_phys, 
    var = "TRIX", 
    group = "Region", 
    title = "TRIX across all stations", 
    ylab = "TRIX index", 
    ordered_latitude = ordered_latitude, 
    ordered_longitude = ordered_longitude
)

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
selections <- sapply(
    names(full_models), 
    function(name) {
        data_reduced <- cleaned_data %>% 
        dplyr::filter(New_basin == name) %>% 
        dplyr::mutate(across(all_of(vars), ~ decostand(., "standardize")))
        new_model <- step(full_models[[name]])
    }
)


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
