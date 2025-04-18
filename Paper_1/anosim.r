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

vars <- c("DO", "NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "T", "pH")
chem_phys %>% dplyr::select(-c(Region, id, Date, E_cond)) %>% apply(2, function(x) shapiro.test(x)$statistic)
vars_to_transform <- c("Chla", "NH4", "NO2", "NO3","PO4", "Salinity", "SiO4", "TN", "TP")

chem_phys.transf <- chem_phys %>% dplyr::select(-c(Secchi_depth)) %>% 
mutate(across(all_of(vars_to_transform), boxcox_transform)) 

ggqqplot(log2(chem_phys$Salinity))
chem_phys.transf$Salinity  %>% hist()
data_anosim <- chem_phys.transf %>% dplyr::select(all_of(c("id", "Date", vars))) %>% na.omit() %>% 
merge(sample_abund %>% dplyr::select(id, Date, New_basin), 
by = c("id", "Date")) %>% 
dplyr::select(-c(id, Date))
test <- vegan::anosim(
    data_anosim  %>% dplyr::select(-c(New_basin)), 
    data_anosim  %>%  pull(New_basin), 
    distance = "mahalanobis",
    permutations = 400, 
    parallel = 4
    )

plot(test)
basin_combinations <- data_anosim %>% pull(New_basin) %>% unique() %>% combn(2)
first_combination <- basin_combinations[1,]
second_combination <- basin_combinations[2, ]


post_hoc_anosim <- mapply(
    function(x, y) {
        vegan::anosim(
            data_anosim %>% dplyr::filter(New_basin %in% c(x, y)) %>% dplyr::select(-c(New_basin)), 
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
) %>% dplyr::filter(PValue > 0.05)

data.frame(
    Basin1 = first_combination, 
    Basin2 = second_combination, 
    R = sapply(post_hoc_anosim, function(x) x$statistic), 
    PValue = sapply(post_hoc_anosim, function(x) x$signif)
) %>% dplyr::filter(PValue < 0.05)
