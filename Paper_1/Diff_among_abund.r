library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
library(vegan)
library(rjson)
library(openxlsx)
library(colorBlindness)

HOME_ <- "."
IMAGE_FORMAT <- "pdf"
source(file.path(HOME_, "utils.r"))

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% select(id, Transect,SeaDepth)
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

#are there differences aming regions?
anova(
    lm(log10(Num_cell_l + 1) ~ Region, data = phyto_abund)
)

phyto_abund %>% 
dplyr::filter(
    Basin == "SAR", 
    Class == "Cryptophyceae",
    #Genus %in% c("Plagioselmis"),
    #Season != "Summer"
     ) %>% 
group_by(Date, id) %>% 
summarise(
    abund = sum(Num_cell_l, na.rm = TRUE),
    Season = first(Season), 
    .groups = "drop"
) %>% 
mutate(
    Season = case_when(
        Season == "Winter" ~ "Winter", 
        Season == "Spring" ~ "Winter",
        Season == "Summer" ~ "Else",
        Season == "Autumn" ~ "Else" 
    )
) %>% 
    lm(log10(abund + 1) ~ Season, data = .) %>% 
    anova()

phyto_abund %>% 
dplyr::filter(Basin == "SAR", Class %in% c("Coccolithophyceae")) %>% 
group_by(Date, id) %>% 
summarise(
    abund = sum(Num_cell_l, na.rm = TRUE),
    Season = first(Season), 
    .groups = "drop"
) %>%
#mutate(
#    Season = case_when(
#        Season == "Winter" ~ "Cold", 
#        Season == "Spring" ~ "Warm",
#        Season == "Summer" ~ "Warm",
#        Season == "Autumn" ~ "Cold" 
#    )
#) %>% 
ggplot() + 
geom_boxplot(aes(x = Season, y = log10(abund + 1))) 


phyto_abund %>% 
dplyr::filter(
    Basin == "SAR", 
    Class == "Coccolithophyceae",
    #Genus %in% c("Plagioselmis"),
    #Season != "Summer"
     ) %>% 
group_by(Date, id) %>% 
summarise(
    abund = sum(Num_cell_l, na.rm = TRUE),
    Season = first(Season), 
    .groups = "drop"
) %>% 
group_by(Season) %>% 
summarise(
    abund = mean(abund, na.rm = TRUE)
) %>% arrange(desc(abund))

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
  phyto_abund %>% dplyr::distinct(id, Date, Closest_coast, SeaDepth, Season, Basin, Transect), 
  by = c("Date", "id")
) 
chem_phys$TRIX <- compute_TRIX(chem_phys)



ids <- phyto_abund %>% dplyr::filter(Region %in% c("SAR")) %>% pull(id) %>% unique()
chem_phys %>% dplyr::filter(id %in% ids) %>% 
dplyr::filter(NP_tot < 1000) %>%
pivot_longer(
  cols = c("NH4", "NO3", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "TRIX", "NP_tot", "DO", "DIN"),
  names_to = "Variable",
  values_to = "Value"
) %>% 
mutate(id = factor(id, levels = params$ordered_id[params$ordered_id %in% ids], ordered = TRUE)) %>%
ggplot() + 
geom_boxplot(aes(x = id, y = Value, fill = Transect)) + 
facet_wrap(~ Variable, scales = "free")


sample_abund <- phyto_abund %>%
    group_by(Date, id) %>%
    summarise(
        Abund = sum(Num_cell_l),
        Basin = first(Basin),
        Region = first(Region), 
        Season = first(Season), 
        Closest_coast = first(Closest_coast),
        Longitude = first(Longitude), 
        Latitude = first(Latitude), 
        .groups = "drop"
    )
sample_abund$Region <- factor(sample_abund$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
sample_abund$Basin <- factor(sample_abund$Basin, levels = ordered_basins, ordered = TRUE)
sample_abund$id <- factor(sample_abund$id, levels = params$ordered_id, ordered = TRUE)
sample_abund <- merge(
    sample_abund, 
    sea_depth %>% select(id, SeaDepth, Transect)
)
sample_abund$Transect <- factor(sample_abund$Transect, levels = ordered_transect, ordered = TRUE)
sample_abund <- sample_abund %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) 

# Compute distance matrix between MAR and MOL regions based on log10 sample abundances


abund_data <- sample_abund %>%
    filter(Region %in% c("MAR", "CAL")) %>%
    select(Date, id, Region, Abund) %>%
    mutate(Abund = log10(Abund + 1)) %>% 
    mutate(sample = paste(Date, id, sep = "_"))  %>% 
    select(-Date, -id)
vars <- c("DO", "DIN", "PO4", "TN", "TP", "T", "pH")
samples_chem <- chem_phys %>% 
dplyr::filter(Region %in% c("MAR", "CAL")) %>%
dplyr::select(all_of(
    c("id", "Date", "Region", vars)
)) %>% na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) %>%
mutate(sample = paste(Date, id, sep = "_"))

common_samples <- intersect(abund_data$sample, samples_chem$sample)

# Compute distance matrix (Euclidean by default)
dist_matrix <- dist(abund_data %>%
                    dplyr::filter(sample %in% common_samples) %>% 
                    dplyr::select(-sample, -Region) %>% 
                    vegan::decostand(method = "standardize"),
                    method = "euclidean")

chem_dist <- samples_chem %>% 
dplyr::filter(sample %in% common_samples) %>%
dplyr::select(-c(id, Date, Region, sample)) %>% 
mutate(across(all_of(vars), \(x) vegan::decostand(x, method = "standardize"))) %>%
dist(method = "euclidean")

mantel(chem_dist, dist_matrix, method="spearman", permutations=999, strata = NULL,
    na.rm = FALSE, parallel = getOption("mc.cores"))



sample_abund %>% dplyr::filter(Region == "SAR") %>%
merge(
    chem_phys %>% dplyr::select(-Region),
    by = c("Date", "id")
) %>% ggplot() + 
geom_point(aes(x = T, y = log10(Abund + 1), colour = Region)) + 
geom_smooth(aes(x = T, y = log10(Abund + 1)), method = "lm", se = FALSE)

sample_abund %>% 
dplyr::filter(Region %in% c("CAL", "CAM")) %>% 
mutate(year = as.numeric(format(as.Date(Date), "%Y"))) %>%
group_by(Region, year) %>% 
summarise(
    sample_per_year = n()
)
