library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(ggvegan)
library(MASS)
library(dplyr)
library(openxlsx)
library(rjson)
library(patchwork)
library(ggpp)
library(ggrepel)

IMAGE_FORNMAT <- "pdf"
HOME_ <- "."
source(file.path(HOME_, "utils.r"))

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
    expl, 
    scale = TRUE
    )
  } else {
    mod <- rda(resp ~ ., expl)
  }
  return(mod)
}

compute_cca_model <- function(resp, expl, covariates, partial = FALSE, partial_covariates = NULL) {
  if (partial) {
    mod <- cca(
      formula(
        paste(
          "resp ~", paste(covariates, collapse = " + "), 
          "+ Condition(", 
          paste(partial_covariates, collapse = " + "),
          ")",
          sep = " "
        )
      ), 
    expl, 
    scale = TRUE
    )
  } else {
    mod <- cca(resp ~ ., expl)
  }
  return(mod)
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

log_trans <- function(x) {
    eps <- x[x != 0] %>% na.omit() %>% min()
    return(as.numeric(log10(x + eps)))
}

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
    expl <- RDA.data  %>% dplyr::select(all_of(c(covariates, mems_cols))) #%>% decostand(method = "standardize")
    return(list(resp = resp, expl = expl, label = RDA.data %>% dplyr::select(all_of(c("id", "Date")))))
 }

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))

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


abund_only_genera <- phyto_abund %>%
        dplyr::filter(Det_level %in% c("Genus", "Species")) %>%
        group_by(id, Date, Genus) %>% 
        summarise(
            abund = sum(Num_cell_l),
            Region = first(Region),
            Season = first(Season), 
            .groups = "drop"
        )


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





chem_phys %>% dplyr::filter(Basin == "NA") %>% ggplot() + 
geom_boxplot(aes(x = Season, y = PO4, fill = Season)) + 
facet_wrap(~Region, ncol = 1, scales = "free")

vars_to_transform <- c("NH4", "NO3","PO4", "Salinity", "SiO4", "TN", "TP", "NP_tot", "pH")
chem_phys <- chem_phys %>% mutate(across(all_of(vars_to_transform), boxcox_transform)) %>% 
dplyr::select(-c(O_sat, NP)) %>%
na.omit() %>% 
dplyr::filter(if_all(where(is.numeric), ~ is.finite(.))) 


sheets <- getSheetNames(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"))
IndVal <- sapply(sheets, function(sheet) {
    data <- read.xlsx(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), sheet = sheet)
    colnames(data)[1] <- "Taxon"
    return(data)
}, simplify = FALSE)
names(IndVal) <- sheets

mems <- read.csv(file.path(HOME_, "morans.csv"))
mems_cols <- colnames(mems)[-1][mems[1, -1] > 0]
positive_mems <- mems %>% dplyr::select(
    id, all_of(mems_cols)
)


vars <- c("DO", "NH4", "NO3", "PO4", "Salinity", "T","pH", "SiO4", "Season")
vars_to_transform <- c("NH4", "NO3","PO4", "SiO4", "TN", "TP", "pH")
covariates <- chem_phys


plot(variation_partitioning)

showvarparts(4)

run_variation_partitioning <- function(threshold, log_transform = TRUE, transfo = "hellinger") {
    relevant_genera <- sapply(
        IndVal,
        function(data) {
            data %>% dplyr::mutate(across(where(is.numeric), ~ .^2)) %>% 
                dplyr::filter(apply(dplyr::select(., where(is.numeric)), 1, max, na.rm = TRUE) > threshold) %>% 
                pull(Taxon)
        }
    ) %>% unlist() %>% unique()
        
    sites_genera <- abund_only_genera %>% dplyr::filter(Genus %in% relevant_genera) %>% pivot_wider(
        id_cols = c(Region, Season, id, Date),
        names_from = Genus,
        values_from = abund,
        values_fill = 0
    )
    
   sufficient_genera <- sites_genera %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
      rowSums() %>%
      `>`(floor(length(relevant_genera) / 10))
    
    sites_genera <- sites_genera[sufficient_genera, ]

    ordination.data <- create_resp_expl_data(
        covariates, 
        sites_genera, 
        positive_mems, 
        vars,
        relevant_genera, 
        log_transform = log_transform, 
        hellinger = FALSE
    )
    
    variation_partitioning <- vegan::varpart(
        ordination.data$resp, 
        ordination.data$expl[, mems_cols], 
        ordination.data$expl[, setdiff(colnames(ordination.data$expl), c(mems_cols,"Season"))], 
        as.factor(ordination.data$expl[, "Season"]),
        transfo = transfo, 
        chisquare = ifelse(transfo == "chi.square", TRUE, FALSE),
        scale = FALSE
    )

    result <-  data.frame(
        Adj.R.squared = variation_partitioning[["part"]][["indfract"]][, "Adj.R.square"], 
        Log = ifelse(log_transform, "log", "no_log"),
        N_genera = length(relevant_genera),
        Transfo = ifelse(transfo == "hellinger", "RDA", "CCA"),
        Term = c("MeMs", "Env", "Season","MeMs+ Env", "Envs + Season", "MeMs + Season", "Mems + Env + Season", "Residuals")
    )
    return(result)
}


results <- mapply(
    run_variation_partitioning,
    threshold = rep(rep(c(0.25, 0.5), each = 2), times = 2),
    log_transform = rep(c(TRUE, FALSE), each = 4),
    transfo = rep(c("hellinger", "chi.square"), times = 4),
    SIMPLIFY = FALSE
) %>% bind_rows()



results %>% dplyr::mutate(
    Adj.R.squared = ifelse(Adj.R.squared < 0, 0, Adj.R.squared), 
    Term = factor(Term, levels = c("MeMs", "Env", "Season", "MeMs+ Env", "Envs + Season", "MeMs + Season", "Mems + Env + Season", "Residuals"))
) %>% pivot_wider(
    id_cols = Term,
    names_from = c("Log", "N_genera", "Transfo"),
    values_from = Adj.R.squared
) %>% write.csv(
    file = file.path(HOME_, paste("variation_partitioning_results.csv", sep = "")), 
    row.names = FALSE
)

p <- results %>% dplyr::mutate(
    model = paste(Log, N_genera, Transfo, sep = "_"), 
    Adj.R.squared = ifelse(Adj.R.squared < 0, 0, Adj.R.squared)
) %>% ggplot() + 
    geom_bar(aes(x = model, y = Adj.R.squared, fill = Term), stat = "identity", position = "fill") + 
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
        axis.text.y = element_text(angle = 0, hjust = 1, size = 15), 
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18), 
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 12)
    ) +
    labs(
        x = "Model", 
        y = "Adj. R-squared", 
        fill = "Term", 
        title = "Variation Partitioning Results"
    ) 
p
ggsave(
    filename = file.path(HOME_, paste("variation_partitioning_results.", IMAGE_FORNMAT, sep = "")), 
    plot = p,
    width = 12,
    height = 8,
    dpi = 300
)


relevant_genera <- sapply(
        IndVal,
        function(data) {
            data %>% dplyr::mutate(across(where(is.numeric), ~ .^2)) %>% 
                dplyr::filter(apply(dplyr::select(., where(is.numeric)), 1, max, na.rm = TRUE) > 0.5) %>% 
                pull(Taxon)
        }
    ) %>% unlist() %>% unique()



sites_genera <- abund_only_genera %>% dplyr::filter(Genus %in% relevant_genera) %>% pivot_wider(
        id_cols = c(Region, Season, id, Date),
        names_from = Genus,
        values_from = abund,
        values_fill = 0
    )

sufficient_genera <- sites_genera %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
      rowSums() %>%
      `>`(floor(length(relevant_genera) / 10))
    
sites_genera <- sites_genera[sufficient_genera, ]


ordination.data <- create_resp_expl_data(
  covariates, 
  sites_genera, 
  positive_mems, 
  vars,
  relevant_genera, 
    log_transform = TRUE, 
    hellinger = TRUE
)


CCA <- compute_cca_model(
    ordination.data$resp, 
    ordination.data$expl, 
    covariates = vars, 
    partial = TRUE, 
    partial_covariates = mems_cols
  )


RDA <- compute_rda_model(
    ordination.data$resp, 
    ordination.data$expl, 
    covariates = vars, 
    partial = TRUE, 
    partial_covariates = mems_cols
  )



mems_vif <- vif.cca(RDA)[mems_cols]
selected_mems_cols <- names(mems_vif[mems_vif < 10])
RDA <- compute_rda_model(
    ordination.data$resp, 
    ordination.data$expl, 
    covariates = vars, 
    partial = TRUE, 
    partial_covariates = selected_mems_cols
  )

null_model <- rda(
    as.formula(paste("ordination.data$resp ~ 1 + Condition(", 
          paste(selected_mems_cols, collapse = " + "),
          ")", sep = "")), 
    ordination.data$expl
)


step_mod <- ordiR2step(null_model, scope = formula(RDA), direction = "backward", R2scope = TRUE, permutations = 999)

model <- step_mod
sites <- cbind(
  ordination.data$label,
  scores(model, display = "sites", scaling = 2)
) %>% merge(
  phyto_abund %>% dplyr::select(id, Date, Season) %>% distinct(),
  by = c("id", "Date")
)
env_arrows <- scores(model, display = "bp", scaling = 2) %>% as.data.frame()
centroids <- scores(model, display = "cn", scaling = 2) %>% as.data.frame()
rownames(centroids) <- gsub("^Season", "", rownames(centroids))
species <- scores(model, display = "species", scaling = 2) %>% as.data.frame() %>% 
rownames_to_column("Genus") %>%
dplyr::filter(Genus %in% important_taxa) %>%
 mutate(
            Genus = case_when(
                Genus == "Thalassionema" ~ "Tham", 
                Genus == "Thalassiosira" ~ "Thar",
                TRUE ~ Genus
            )
        )
species$Genus <- substr(species$Genus, 1, 4)
perc <- round(100*(summary(model)$cont$importance[2, 1:2]), 2)
important_taxa <- species %>%
  dplyr::mutate(dist = sqrt((.[[1]])^2 + (.[[2]])^2)) %>%
  dplyr::arrange(desc(dist)) %>%
  head(15) %>%
  rownames()
axis1 <- grep("A1", colnames(env_arrows), value = TRUE)
axis2 <- grep("A2", colnames(env_arrows), value = TRUE)

spe_factor <- seq(1, 1, length.out = length(important_taxa))
p_spe <- sites %>% ggplot() + 
geom_point(
    aes(x = !!as.symbol(axis1), 
      y = !!as.symbol(axis2), 
      fill = Season), 
    shape = 21, 
    color = "black", 
    size = 3, 
    alpha = 0.8
    ) + 
geom_segment(
    data = species, 
    aes(x = 0, y = 0, xend = !!as.symbol(axis1) * spe_factor, yend = !!as.symbol(axis2) * spe_factor), 
    arrow = arrow(length = unit(0.2, "cm")), 
    col = "black", 
    linewidth = 1
    ) +
geom_label_repel(
    data = species, 
    aes(x = !!as.symbol(axis1) * spe_factor, y = !!as.symbol(axis2) * spe_factor), 
    position = position_nudge_center(x = 0.01, y = 0.01,
                                                    center_x = 0, center_y = 0),
    label = species %>% pull(Genus),
    size = 6,  
    color = "black", 
    fill = "white"
) +
labs(x = paste(axis1, " (",perc[1], "%)", sep = "") , y = paste(axis2, " (",perc[2], "%)", sep = "")) + 
theme_bw() + 
theme(legend.position = "bottom")


env_factor <- rep(1, nrow(env_arrows) - 3)
env_arrows["T", ] <- env_arrows["T", ] / 2
p_env <- sites %>% ggplot() + 
geom_point(
    aes(x = !!as.symbol(axis1), 
      y = !!as.symbol(axis2), 
      fill = Season), 
    shape = 21, 
    color = "black", 
    size = 3, 
    alpha = 0.8
    ) + 
geom_segment(
  data = env_arrows[seq(4, nrow(env_arrows)), ], 
  aes(x = 0, y = 0, xend = !!as.symbol(axis1) * env_factor, yend = !!as.symbol(axis2) * env_factor), 
  arrow = arrow(length = unit(0.2, "cm")), 
  col = "black", 
  linewidth = 1
  ) +
geom_label_repel(
    data = rbind(env_arrows[seq(4, nrow(env_arrows)), ], centroids),
    aes(
      x = !!as.symbol(axis1) * c(env_factor, c(1, 1, 1, 1)), 
      y = !!as.symbol(axis2) * c(env_factor, c(1, 1, 1, 1)),
      label = c(rownames(env_arrows)[seq(4, nrow(env_arrows))], rownames(centroids))
      ),
    position = position_nudge_center(
      x = 0.005, y = 0.005,
      center_x = 0, center_y = 0
      ),  
    size = 7, fontface = "bold",
    color = "black",
    fill = "white", 
    force = 50
) + 
geom_point(
    data = centroids, 
    aes(x = !!as.symbol(axis1), y = !!as.symbol(axis2)),  
    size = 5, 
    color = "black"
) + 
labs(x = paste(axis1, " (",perc[1], "%)", sep = "") , y = paste(axis2, " (",perc[2], "%)", sep = "")) + 
theme_bw() + 
theme(legend.position = "bottom")
p_env
p <-  p_env / p_spe + plot_layout(ncol = 1) + 
plot_layout(guides = "collect") + 
plot_annotation(title = "RDA of phytoplankton communities in Italy") & 
plot_annotation(
    title = "RDA of phytoplankton communities in Italy", 
    tag_levels = "A", tag_suffix = ")") &
theme(
  legend.title = element_text(face = "bold", size = 15),
  legend.position = "bottom",
  legend.text = element_text(size = 15), 
  plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  ) & 
coord_cartesian(xlim = c(-0.6, 0.5)) 
p
ggsave(
    filename = file.path(HOME_, paste("ordination_all_italy.", IMAGE_FORNMAT, sep = "")), 
    plot = p,
    width = 10,
    height = 12,
    dpi = 300
)

p_all <- sites %>% ggplot() + 
geom_point(
    aes(x = !!as.symbol(axis1), 
      y = !!as.symbol(axis2), 
      fill = Season), 
    shape = 21, 
    color = "black", 
    size = 3, 
    alpha = 0.8
    ) + 
geom_segment(
    data = species[important_taxa, ], 
    aes(x = 0, y = 0, xend = !!as.symbol(axis1) * 5, yend = !!as.symbol(axis2) * 5), 
    arrow = arrow(length = unit(0.2, "cm")), 
    col = "grey", 
    linewidth = 1
    ) +
geom_label(
    data = species[important_taxa, ], 
    aes(x = !!as.symbol(axis1) * 5, y = !!as.symbol(axis2) * 5, 
    label = rownames(species[important_taxa, ])),
    size = 6,  
    color = "black", 
    fill = "white"
) +
geom_segment(data = env_arrows[seq(4, nrow(env_arrows)), ], aes(x = 0, y = 0, xend = !!as.symbol(axis1) * 5, yend = !!as.symbol(axis2) * 5), 
            arrow = arrow(length = unit(0.2, "cm")), col = "black") +
geom_label(
    data = env_arrows[seq(4, nrow(env_arrows)), ],
    aes(x = !!as.symbol(axis1) * 5, y = !!as.symbol(axis2) * 5, label = rownames(env_arrows)[seq(4, nrow(env_arrows))]),
    size = 7, fontface = "bold",
    color = "black",
    fill = "white"
) + 
labs(title = paste("Ordination entire Italy"), x = paste(axis1, " (",perc[1], "%)", sep = "") , y = paste(axis2, " (",perc[2], "%)", sep = "")) + 
theme_bw() + 
theme(legend.position = "bottom")

p_all


#in previous plot was modified
site_scores <- scores(step_mod, display = "sites")

# Correlate environmental variables with RDA axes
cors <- cor(ordination.data$expl[, c("T", "Salinity", "PO4", "pH", "NO3", "DO", "SiO4")], site_scores[, 1:2])

# Normalize each row of the correlation matrix so that the sum of absolute values in each row is 1
cors_normalized <- t(apply(cors, 1, function(x) x / sqrt(sum(x^2))))
cors_normalized
p_load <- cors_normalized %>% as.data.frame() %>%
rownames_to_column("Variable") %>%
mutate(
    Variable = factor(Variable, levels = c(rownames(env_arrows)[seq(4, nrow(env_arrows))], rownames(centroids)), ordered = TRUE)
) %>%
pivot_longer(
    cols = c(axis1, axis2),
    names_to = "Axis", 
    values_to = "Value"
) %>% ggplot() + 
geom_bar(
    aes(y = Variable, x = Value), 
    stat = "identity", position = "dodge"
) + facet_wrap(~Axis, scales = "free_y", ncol = 1) + 
scale_y_discrete(limits = rev) +
labs(
    title = "Corr. of env variables with RDA axes", 
    y = "Variable", 
    x = "Correlation"
) + 
theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 15), 
    axis.text.y = element_text(angle = 0, hjust = 1, size = 15), 
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18), 
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
)
p_load
ggsave(
    filename = file.path(HOME_, paste("ordination_all_italy_env_arrows.", IMAGE_FORNMAT, sep = "")), 
    plot = p_load,
    width = 8,
    height = 9,
    dpi = 300
)

vars <- c("T")

fit <- envfit(step_mod, ordination.data$expl[, c("T", "Salinity", "PO4", "pH", "NO3", "DO", "SiO4", "Season")], permutations = 999)

fit
dev.off()
plot(step_mod, display = c("sites", "species"))
plot(fit, p.max = 0.05, col = "red")
scores(fit, display = "vectors")

0.9833 ^ 2 + 0.1812 ^ 2



fit$vectors$arrows
