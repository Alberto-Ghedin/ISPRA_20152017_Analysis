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
library(openxlsx)


IMAGE_FORNMAT <- "svg"
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

vars <- c("NH4", "NO3", "DO", "PO4", "SiO4", "Salinity", "TN", "TP", "pH", "T", "DIN_TN", "O_sat", "Closest_coast", "SeaDepth")
vars_to_transform <- c("NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "DIN_TN", "pH", "O_sat", "DIN_TN", "P_rat")


chem_phys <- read.csv("./df_chem_phys.csv") %>% 
dplyr::select(-c(Region, E_cond, Secchi_depth, NO2, Chla)) %>%
dplyr::filter(
  DO < 400, 
  NH4 < 10,
  NO3 < 57, 
  O_sat < 160, 
  pH > 7, 
  PO4 < 2, 
  Salinity > 20, 
  SiO4 < 40, 
  TN < 120
  ) %>% 
mutate(
  DIN_TN = (NH4 + NO3) / TN,
  P_rat = PO4 / TP, 
  N_star = NO3 - 16 * PO4
  ) %>% 
mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) %>% 
merge(
  phyto_abund %>% dplyr::distinct(id, Date, Closest_coast, SeaDepth, Season, New_basin), 
  by = c("Date", "id")
)

##Check env per basin, 
pca_basins <- sapply(
  chem_phys %>% pull(New_basin) %>% unique(),
  function(basin) {
    data <- chem_phys %>% 
    dplyr::select(-c(id, Date, P_rat, Closest_coast, SeaDepth, Season)) %>%
    dplyr::filter(New_basin == basin) %>% na.omit() %>% 
    dplyr::select(-New_basin)
    
    return(
      rda(
        data %>% mutate(across(everything(), ~ decostand(., "standardize"))),
        tidy = TRUE
      )
    )
  }, 
  simplify = FALSE
)
names(pca_basins) <- chem_phys %>% pull(New_basin) %>% unique()

dir.create("./pca_env", showWarnings = FALSE)
plot_dir <- file.path(HOME_, "pca_env")
sapply(
  names(pca_basins),
  function(basin) {
    pca <- pca_basins[[basin]]
    sites <- cbind(
      scores(pca, display = "sites"), 
      chem_phys %>% dplyr::filter(New_basin == basin) %>% dplyr::select(id, Date, Season)
    )
    env_arrows <- scores(pca, display = "species")
    perc <- round(100*(summary(pca)$cont$importance[2, 1:2]), 2)
    
    p <- sites %>% ggplot() + 
      geom_point(aes(x = PC1, y = PC2, col = Season)) + 
      geom_segment(data = env_arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                   arrow = arrow(length = unit(0.2, "cm")), col = "blue") +
      geom_text(data = env_arrows, aes(x = PC1, y = PC2, label = rownames(env_arrows)),
                size = 3, hjust = 0.5, vjust = 0.5) +
      labs(title = paste("RDA for", basin), x = paste("PC1 (",perc[1], "%)", sep="") , y=paste("PC2 (",perc[2], "%)", sep="")) + 
      theme_bw() + 
      theme(legend.position="bottom")
    
    ggsave(
      file.path(plot_dir, paste(basin, "pca_env", sep = "_")),
      device = IMAGE_FORNMAT,
      plot = p,
      width = 10,
      height = 10
    )
  }, 
  simplify=FALSE
)

chem_phys %>% head()




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


#MEMS 
sheets <- getSheetNames("./MEMs_per_basin.xlsx")
mems <- sapply(
    sheets,
    function(sheet) {
    data <- read.xlsx("./MEMs_per_basin.xlsx", sheet = sheet)
    }, 
    simplify = FALSE
    ) 
names(mems) <- sheets


create_resp_expl_data <- function(
  env_data,
  sites_taxa,
  mems, 
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
    ) %>% dplyr::filter(rowSums(across(all_of(species))) > 0) %>% 
    merge(
        mems,
        by = "id",
        all = FALSE
    )

    resp <- RDA.data %>% dplyr::select(all_of(species)) %>% mutate(across(everything(), ~ if (log_transform) log1p(.) else .))
    if (!is.null(db)) {
        resp <- vegan::vegdist(resp, method = db, binary = FALSE)
    } else if (hellinger) {
        resp <- vegan::decostand(resp, method = "hellinger") #%>% mutate(across(everything(), ~ . - mean(.)))
    }
    mems_cols <- mems %>% dplyr::select(where(~is.numeric(.))) %>% colnames()
    expl <- RDA.data %>% dplyr::mutate(across(all_of(c(covariates, mems_cols)), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(c(covariates, mems_cols)))
    return(list(resp = resp, expl = expl, label = RDA.data %>% dplyr::select(all_of(c("id", "Date")))))
 }

pull_taxa <- function(data, rate) {
  return(
    setdiff(
      data %>% rowwise() %>%
      dplyr::filter(max(c_across(where(is.numeric))) > rate) %>% 
      ungroup() %>% pull(Taxon), 
      c("Bacillariophyceae", "Dinoflagellata", "Dinophyceae", "Coccolithophyceae", "Cryptophyceae", "Other phytoplankton")
      )
      )
}


variation_partitionings <- sapply(
  names(all_data),
  function(basin) {
    selected_genera <- pull_taxa(all_data[[basin]], 0.5)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  mems = mems[[basin]],
  covariates = vars,
  species = selected_genera,
  log_transform = TRUE, 
  db =  NULL, 
  hellinger = FALSE
)
  mems_col <- RDA.data$expl %>% dplyr::select(dplyr::contains("MEM")) %>% colnames()
  var_par <- vegan::varpart(
    RDA.data$resp, 
    RDA.data$expl[, mems_col], 
    RDA.data$expl[, setdiff(colnames(RDA.data$expl), mems_col)], 
    transfo = "hel"
  )
  return(var_par)
  }, 
  simplify = FALSE
)

make_triplot <- function(model, labels, basin) {
  sites <- cbind(
  scores(model, display = "sites"), 
  labels
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
  
  return(p)
}

compute_rda_model <- function(resp, expl, covariates, partial = FALSE, partial_covariates = NULL) {
  if (partial) {
    mod <- rda(
      formula(
        paste(
          "resp ~", paste(covariates, collapse = " + "), 
          "+ Condition(", 
          paste(partial_covariates, collapse = " + "),
          ")",
          sep = " "
        )
      ), 
    expl
    )
  } else {
    mod <- rda(resp ~ ., expl)
  }
  return(mod)
}

rda_dir <- file.path(HOME_, "rda")
dir.create(rda_dir, showWarnings = FALSE)
models <- sapply(
  names(all_data),
  function(basin) {
    selected_genera <- pull_taxa(all_data[[basin]], 0.5)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  mems = mems[[basin]],
  covariates = vars,
  species = selected_genera,
  log_transform = TRUE, 
  db =  NULL, 
  hellinger = FALSE
)

  resp <- RDA.data$resp
  expl <- RDA.data$expl
  mems_col <- RDA.data$expl %>% dplyr::select(dplyr::contains("MEM")) %>% colnames()
  model <- compute_rda_model(
    resp, 
    expl, 
    vars,
    partial = TRUE, 
    partial_covariates = mems_col
  )
  
  selection <- ordistep(
  model,
  direction = "both",
  Pin = 0.05,
  Pout = 0.1, 
  trace = FALSE
  )

  model <- rda(
    selection$call$formula,
    expl
  )

  labels <- RDA.data$label %>% merge(
    abund_only_genera %>% dplyr::filter(Basin == basin) %>% dplyr::select(id, Date, Season),
    by = c("id", "Date")
    )
  
  p <- make_triplot(model, labels, basin)
  
  ggsave(
      file.path(rda_dir, basin),
      device = IMAGE_FORNMAT,
      plot = p,
      width = 10,
      height = 10
    )
  return(model)
  }, 
  simplify = FALSE
)

models <- sapply(
  names(all_data),
  function(basin) {
    selected_genera <- pull_taxa(all_data[[basin]], 0.5)
    RDA.data <- create_resp_expl_data(
  env_data = chem_phys,
  sites_taxa = abund_only_genera %>% dplyr::filter(Basin == basin),
  mems = mems[[basin]],
  covariates = vars,
  species = selected_genera,
  log_transform = TRUE, 
  db =  NULL, 
  hellinger = FALSE
)

  mems_col <- RDA.data$expl %>% dplyr::select(dplyr::contains("MEM")) %>% colnames()
  model <- compute_rda_model(
    RDA.data$resp, 
    RDA.data$expl, 
    vars,
    partial = TRUE, 
    partial_covariates = mems_col
  )

  
}, 
  simplify = FALSE
)


do.call(
  rbind, 
  sapply(
    names(models), 
    function(basin) {
      scores(models[[basin]], display = "bp") %>%
      as.data.frame() %>%
      rownames_to_column(var = "Var") %>%
      mutate(
        Basin = basin
      )
    }, 
    simplify = FALSE
  )
) %>% pivot_longer(
  cols = c("RDA1", "RDA2"), 
  names_to = "Axis", 
  values_to = "Value"
) %>% ggplot() + 
geom_tile(aes(x = Var, y = Basin, fill = Value)) +
facet_wrap(~Axis, scales = "free")
