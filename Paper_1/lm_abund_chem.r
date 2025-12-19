library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(tibble)
library(openxlsx)
library(rjson)



IMAGE_FORMAT <- "pdf"
HOME_ <- "."

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

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% dplyr::select(id, Transect,SeaDepth)
)
phyto_abund$Region <- from_region_to_abreviation[as.character(phyto_abund$Region)]
phyto_abund$Transect <- factor(phyto_abund$Transect, levels = ordered_transect, ordered = TRUE)
phyto_abund <- phyto_abund %>% mutate(
    Basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Transect %in% c("FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO") ~ "SA",
        Transect %in% c("PORTO_CESAREO", "PUNTA_RONDINELLA") ~ "SM",
        Region == "BAS" ~ "SM",
        Transect %in% c("Villapiana", "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche") ~ "SM",
        Region == "SIC" ~ "SIC", 
        Transect %in% c("Vibo_marina", "Cetraro") ~ "ST",
        Region == "CAM" ~ "ST",
        Transect %in% c("m1lt01", "m1lt02") ~ "ST",
        Transect %in% c("m1rm03", "m1vt04") ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
)
phyto_abund$Basin <- factor(phyto_abund$Basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)



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


chem_phys <- chem_phys %>% dplyr::select(-c(E_cond, Secchi_depth, NO2)) %>%
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
  ) %>% merge(
  phyto_abund %>% dplyr::distinct(id, Date, Closest_coast, SeaDepth, Season, Basin), 
  by = c("Date", "id")
) 
chem_phys$TRIX <- compute_TRIX(chem_phys)


vars_to_transform <- c("NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "NP_tot", "pH")
chem_phys <- chem_phys %>% mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
dplyr::select(-c(O_sat, NP)) %>%
na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) 

chem_phys %>% colnames()
vars <- c("NH4", "NO3", "PO4", "SiO4", "Salinity", "T", "TRIX", "NP_tot")
chem_phys %>%
dplyr::filter(Region == "SAR") %>%
merge(phyto_abund %>% dplyr::select(id, Date, Transect), by = c("id", "Date")) %>%
mutate(id = factor(id, levels = params$ordered_id[which(id %in% params$ordered_id)], ordered = TRUE)) %>%
dplyr::filter(SiO4 < 5, PO4 < 0.15) %>%
pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>%
ggplot() + 
geom_boxplot(aes(x = id, y = Value, fill = Transect)) + 
facet_wrap(~Variable, scales = "free")

dev.off()
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


# Compute 5 quantiles (min, 25%, median, 75%, max) for each variable in vars
quantiles_df <- sapply(vars, function(v) {
    quantile(chem_phys[[v]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
}) %>% 
    t() %>% 
    as.data.frame() %>% 
    setNames(c("Min", "Q1", "Median", "Q3", "Max")) %>% 
    tibble::rownames_to_column("Variable")

    # Assign quantile class for each variable's median per season and basin
    vars <- c("T", "NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "pH", "NP")
    summary_by_basin_season <- chem_phys %>%
        mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
        group_by(Basin, Season) %>%
        summarise(across(all_of(vars), ~ median(., na.rm = TRUE), .names = "median_{.col}")) %>%
        ungroup()

    # For each variable, assign quantile class based on quantiles_df
    for (v in vars) {
        summary_by_basin_season[[paste0(v, "_class")]] <- cut(
            summary_by_basin_season[[paste0("median_", v)]],
            breaks = quantiles_df[quantiles_df$Variable == v, c("Min", "Q1", "Median", "Q3", "Max")] %>% as.numeric(),
            labels = c("1", "2", "3", "4"),
            include.lowest = TRUE, right = TRUE
        )
    }

    data.frame(
        summary_by_basin_season %>% dplyr::select(Basin, Season, ends_with("_class")),
        stringsAsFactors = FALSE
    ) %>% pivot_longer(
        cols = -c(Basin, Season), 
        names_to = c("Variable", "Class"), 
        names_sep = "_", 
        values_to = "Quantile_Class"
    ) %>% ggplot() + 
    geom_tile(aes(x = Variable, y = Basin, fill = Quantile_Class)) +
    scale_y_discrete(limits = rev) +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~Season, scales = "free_y")

    # Save the result
    write.csv(summary_by_basin_season, file.path(HOME_, "median_class_by_basin_season.csv"), row.names = FALSE)

chem_phys %>% colnames()
vars <- c("T", "NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "pH", "NP")
chem_phys %>%
mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
    group_by(Basin, Season) %>%
    summarise(across(all_of(vars), 
    list(
        qt1 = ~ quantile(., 0.25, na.rm = TRUE),
        median = ~ median(., na.rm = TRUE),
        qt3 = ~ quantile(., 0.75, na.rm = TRUE),
        min = ~ min(., na.rm = TRUE),
        max = ~ max(., na.rm = TRUE)
        )
        )
        ) %>%
    ungroup() %>% 
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    write.csv(
        file.path(HOME_, "chem_phys_summary.csv"), 
        row.names = FALSE
    )

p<- chem_phys %>%
mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)) %>%
ggplot() + 
geom_boxplot(aes(x = TRIX, y = Season)) +
scale_y_discrete(limits = rev) + 
    facet_wrap(~Basin, scales = "free_y", ncol = 1)
ggsave(
    file.path(HOME_, paste("TRIX_boxplot", IMAGE_FORMAT, sep = ".")), 
    p, 
    width = 10, 
    height = 18, 
    dpi = 300
)


print(paste("Correlation between log10(NO3) and log10(NH4):", round(correlation, 3)))
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


cleaned_data <- sample_abund %>% dplyr::select(Basin, Region, Date, id, sample_abund) %>% 
merge(
    chem_phys %>% dplyr::mutate(
    `log_NH4` = log(NH4),
    `log_NO3` = log(NO3), 
    `log_PO4` = log(PO4),
    `log_SiO4` = log(SiO4),
) %>% dplyr::select(
    id, Date, `log_NH4`, `log_NO3`, `log_PO4`, `log_SiO4`, `Salinity`, `T`, `DO`), 
    by = c("id", "Date")
) %>% na.omit() %>% dplyr::filter(if_all(where(is.numeric), ~ is.finite(.)))

vars <- c("log_NH4", "log_NO3", "log_PO4", "log_SiO4", "Salinity",  "T", "DO")
full_model_all_italy <- lm(
    paste(
        "log10(sample_abund) ~", 
        paste( 
            vars
            ,
            collapse = " + "
        )
    ), 
    data = cleaned_data %>% dplyr::mutate(across(all_of(vars), ~ decostand(., "standardize")))
)

coeffs_all_italy <- summary(step(full_model_all_italy))$coefficients[-1, c(1,4)] %>% as.data.frame() %>% 
rownames_to_column(var = "Var") %>% 
rename(p_value = `Pr(>|t|)`)


selected_model_per_basin <- sapply(
    cleaned_data %>% pull(Basin) %>% unique(), 
    function(name) {
        data_reduced <- 
        model <- lm(
            paste(
                "log10(sample_abund) ~", 
                paste(
                    vars,
                    collapse = " + "
                )
            ), 
            data = cleaned_data %>% 
        dplyr::filter(Basin == name) %>% 
        dplyr::mutate(across(all_of(vars), ~ decostand(., "standardize")))
        )
        return(step(model))
        }, 
    simplify = FALSE
)
names(selected_model_per_basin) <- cleaned_data %>% pull(Basin) %>% unique()

coeffs_per_basin <- sapply(
    names(selected_model_per_basin), 
    function(name) {
        return(
            summary(selected_model_per_basin[[name]])$coefficients[-1, c(1,4)] %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Var") %>% 
            rename(p_value = `Pr(>|t|)`) %>% 
            mutate(Basin = name)
        )
    }, 
    simplify = FALSE
)   
coeffs_per_basin


p <- bind_rows(
    bind_rows(coeffs_per_basin), 
    coeffs_all_italy %>% mutate(
        Basin = "Italy"
    )
) %>% 
mutate(
    Basin = factor(Basin, levels = c("Italy", "NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE), 
    Var = case_when(
        Var == "log_NH4" ~ "log(NH4)", 
        Var == "log_NO3" ~ "log(NO3)", 
        Var == "log_PO4" ~ "log(PO4)", 
        Var == "log_SiO4" ~ "log(SiO4)", 
        TRUE ~ Var
    )
) %>% ggplot() + 
geom_col(aes(x = Basin, y = Estimate, fill = as.factor(sign(Estimate))), position = "dodge", color = "black") +
facet_wrap(~Var, ncol = 3) +
labs(x = "Basin", title = "Standardized coefficients of linear models \n predicting log10(Abundance)") +
theme(
    axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 18), 
    axis.title =  element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5), 
    legend.position = "none"
)
ggsave(
    file.path(HOME_, paste("LM_abundance_per_basin_coefficients", IMAGE_FORMAT, sep = ".")), 
    p, 
    width = 10, 
    height = 12, 
    dpi = 300
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



