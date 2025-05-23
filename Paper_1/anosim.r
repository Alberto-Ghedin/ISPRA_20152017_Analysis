library(vegan)
library(dplyr)
library(tidyr)
library(ggvegan)
library(ggpubr)

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

chem_phys <- read.csv("./df_chem_phys.csv") %>% dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120)
chem_phys$Region <- from_region_to_abreviation[chem_phys$Region]
chem_phys$Region <- factor(chem_phys$Region, levels = unname(from_region_to_abreviation))

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

chem_phys <- chem_phys %>% mutate(NO_rat = NO2 / NO3, 
                    DIN_TN = (NH4 + NO3) / TN,
                    P_rat = PO4 / TP
)
chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% apply(2, function(x) shapiro.test(x[is.finite(x)])$statistic)
vars_to_transform <- c("Chla", "NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "DIN_TN", "pH", "O_sat")

data_anosim <- chem_phys %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, NO_rat, P_rat))  %>% 
    dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120) %>% 
    mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
merge(sample_abund %>% dplyr::select(id, Date, New_basin), 
by = c("id", "Date")) %>% 
dplyr::select(-c(id, Date, Chla, DIN_TN)) %>% na.omit()



test <- vegan::anosim(
    data_anosim  %>% dplyr::select(-c(New_basin)) %>% vegan::decostand(method = "standardize"), 
    data_anosim  %>%  pull(New_basin), 
    distance = "euclidean",
    permutations = 400, 
    parallel = 4
    )
summary(test)
plot(test)
basin_combinations <- data_anosim %>% pull(New_basin) %>% unique() %>% combn(2)
first_combination <- basin_combinations[1,]
second_combination <- basin_combinations[2, ]


post_hoc_anosim <- mapply(
    function(x, y) {
        vegan::anosim(
            data_anosim %>% dplyr::filter(New_basin %in% c(x, y)) %>% dplyr::select(-c(New_basin)) %>% vegan::decostand(method = "standardize"), 
            data_anosim %>% dplyr::filter(New_basin %in% c(x, y)) %>% pull(New_basin), 
            distance = "mahalanobis",
            permutations = 400, 
            parallel = 6
        )
    },
    first_combination, 
    second_combination, 
    SIMPLIFY = FALSE
)

data.frame(
    Basin1 = first_combination, 
    Basin2 = second_combination, 
    R = sapply(post_hoc_anosim, function(x) x$statistic), 
    PValue = sapply(post_hoc_anosim, function(x) x$signif)
) %>% ggplot() + 
geom_tile(aes(x = Basin1, y = Basin2, fill = ifelse(PValue < 0.05, "<0.05", ">0.05")))

data.frame(
    Basin1 = first_combination, 
    Basin2 = second_combination, 
    R = sapply(post_hoc_anosim, function(x) x$statistic), 
    PValue = sapply(post_hoc_anosim, function(x) x$signif)
) %>% ggplot() + 
geom_tile(aes(x = Basin1, y = Basin2, fill = R)) + 
geom_text(aes(x = Basin1, y = Basin2, label = round(R, 2)), color = "white")

data.frame(
    Basin1 = first_combination, 
    Basin2 = second_combination, 
    R = sapply(post_hoc_anosim, function(x) x$statistic), 
    PValue = sapply(post_hoc_anosim, function(x) x$signif)
) %>% dplyr::filter(PValue < 0.05)
