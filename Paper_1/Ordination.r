library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
library(vegan)
library(ggvegan)
library(readxl)
library(indicspecies)
library(MASS)
library(dplyr)

IMAGE_FORNMAT <- "svg"

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

HOME_ <- "."

log_trans <- function(x) {
    eps <- x[x != 0] %>% na.omit() %>% min()
    return(as.numeric(log10(x + eps)))
}

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv"))
phyto_abund <- phyto_abund %>% 
mutate(
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
vars_to_transform <- c("Chla", "NH4", "NO2", "NO3","PO4", "Salinity", "SiO4", "TN", "TP")

chem_phys <- read.csv("./df_chem_phys.csv") %>% dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla, O_sat)) %>% 
    dplyr::filter(pH > 7, PO4 < 2, TN / TP < 120) %>% na.omit() %>% 
    mutate(across(all_of(c("NH4", "NO3", "PO4", "SiO4", "TN", "TP", "Salinity")), boxcox_transform))

chem_phys %>% dplyr::select(-c(Date, id)) %>% pivot_longer(cols = vars, names_to = "Variable", values_to = "Value") %>%
ggplot(aes(x = Value)) + geom_histogram() + facet_wrap(~Variable, scales = "free")

sheets <- excel_sheets("./indval_only_genera_per_basin.xlsx")
all_data <- sapply(sheets, function(sheet) {
    data <- read_excel("./indval_only_genera_per_basin.xlsx", sheet = sheet, col_names = TRUE)
    colnames(data)[1] <- "Taxon"
    data <- data %>% dplyr::filter(Taxon != "Other phytoplankton")
    return(data)
}, simplify = FALSE)
names(all_data) <- sheets

abund_only_genera <- phyto_abund %>% dplyr::filter(Taxon != "Other phytoplankton") %>% mutate(
    Genus = case_when(
        Class == "nan" ~ Taxon,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
) %>% group_by(Date, id, Genus) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    Basin = first(New_basin),
    Season = first(Season),
    .groups = "drop"
) %>% pivot_wider(names_from = Genus, values_from = Abund, values_fill = 0)



create_resp_expl_data <- function(
  env_data,
  sites_taxa,
  covariates, 
  species,
  log_transform = FALSE, 
  db = NULL, 
  hellinger = TRUE
 ) {
    RDA.data <- sites_taxa %>% dplyr::select(all_of(c("id", "Date", species)))  %>% 
    merge(
        env_data,
        by = c("id", "Date"),
        all = FALSE
    ) %>% dplyr::filter(rowSums(across(all_of(species))) > 0)

    resp <- RDA.data %>% dplyr::select(all_of(species)) %>% mutate(across(everything(), ~ if (log_transform) log1p(.) else .))
    if (!is.null(db)) {
        resp <- vegan::vegdist(resp, method = db, binary = FALSE)
    } else if (hellinger) {
        resp <- vegan::decostand(resp, method = "hellinger") #%>% mutate(across(everything(), ~ . - mean(.)))
    }
    expl <- RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(c(covariates)))
    return(list(resp = resp, expl = expl, label = RDA.data %>% dplyr::select(all_of(c("id", "Date")))))
 }




n_row <- abund_only_genera %>% dplyr::filter(Basin == basin) %>% nrow()
selected_genera <- abund_only_genera %>% dplyr::filter(Basin == basin) %>% dplyr::select(where(~is.numeric(.))) %>% 
dplyr::select(where(~ sum(. != 0) / n_row > 0.05)) %>% colnames()

env_datas <- sapply(
  names(all_data),
  function(basin) {
    selected_genera <- all_data[[basin]] %>% rowwise() %>%
  dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
  ungroup() %>% pull(Taxon)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  covariates = vars,
  species = setdiff(selected_genera, c("Bacillariophyceae", "Dinoflagellata", "Dinophyceae")),
  log_transform = FALSE, 
  db =  NULL, 
  hellinger = TRUE
)
  return(RDA.data$expl)
  }, 
simplify = FALSE
)
names(env_datas) <- names(all_data)

RDA.data$expl %>% head()
pca.env <- rda(RDA.data$expl)
pca.env$species
biplot(pca.env, scaling = 2)
RDA.data$resp
pca.spec <- rda(RDA.data$resp)
summary(pca.spec)
model <- dbrda(RDA.data$resp ~ ., data = RDA.data$expl)
summary(model)
vegan::RsquareAdj(model)

model <- cca(RDA.data$resp ~ ., data = RDA.data$expl)
summary(model)
ordiplot(model, scaling = 2, type = "text")



cca_dir <- file.path(HOME_, "cca")
dir.create(cca_dir, showWarnings = FALSE)
sapply(
  names(all_data),
  function(basin) {
    selected_genera <- all_data[[basin]] %>% rowwise() %>%
  dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
  ungroup() %>% pull(Taxon)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  covariates = vars,
  species = setdiff(selected_genera, c("Bacillariophyceae", "Dinoflagellata", "Dinophyceae", "Coccolithophyceae", "Cryptophyceae", "Other phytoplankton")),
  log_transform = TRUE, 
  db =  NULL, 
  hellinger = FALSE
)
  model <- rda(RDA.data$resp ~ ., data = RDA.data$expl, tidy = TRUE)
  sites <- cbind(
  scores(model, display = "sites"), 
  RDA.data$label) %>% merge(
    abund_only_genera %>% dplyr::select(id, Date, Season),
    by = c("id", "Date")
  )
  env_arrows <- scores(model, display = "bp")
  species <- scores(model, display = "species")
  perc <- round(100*(summary(model)$cont$importance[2, 1:2]), 2)

  p <- sites %>% ggplot() + 
  geom_point(aes(x = RDA1, y = RDA2, col = Season)) + 
  geom_segment(data = species, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), col = "grey") +
  geom_text(data = species, aes(x = RDA1, y = RDA2, label = rownames(species)),
            size = 3, hjust = 0.5, vjust = 0.5) +
  geom_segment(data = env_arrows, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), col = "blue") +
  geom_text(data = env_arrows, aes(x = RDA1, y = RDA2, label = rownames(env_arrows)),
            size = 3, hjust = 0.5, vjust = 0.5) +
  labs(title = paste("RDA for", basin, "basin"), x = paste("RDA1 (",perc[1], "%)", sep = "") , y = paste("RDA2 (",perc[2], "%)", sep = "")) + 
  theme_bw() + 
  theme(legend.position = "bottom")
  ggsave(
      file.path(cca_dir, basin),
      device = IMAGE_FORNMAT,
      plot = p,
      width = 10,
      height = 10
    )
  }
)

abund_only_genera %>% head()

basin <- "NA"
selected_genera <- all_data[[basin]] %>% rowwise() %>%
  dplyr::filter(max(c_across(where(is.numeric))) > 0.5) %>%
  ungroup() %>% pull(Taxon)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  covariates = vars,
  species = setdiff(selected_genera, c("Bacillariophyceae", "Dinoflagellata", "Dinophyceae", "Other phytoplankton")),
  log_transform = TRUE, 
  db =  NULL, 
  hellinger = FALSE
)
apply(RDA.data$resp, 1, sum)
RDA.data$resp
RDA.data$expl
model <- rda(RDA.data$resp ~ ., data = RDA.data$expl, tidy = TRUE)






fitted_species <- fitted(model)


