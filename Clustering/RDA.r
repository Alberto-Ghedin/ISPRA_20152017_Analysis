library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(tibble)
library(zoo)
library(vegan)
library(jsonlite)
library(MASS)
library(tidyverse)
HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

#CHOOSING THE ENV_DATA  
env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) %>% dplyr::select(-c(Secchi_depth))
env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Cop_variables_on_sampling_sites.csv", sep = "/")) %>% dplyr::select(-c(Latitude, Longitude, Closest_coast)) 
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/sites_taxa_matrix.csv", sep = "/"))
sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
names(indval_list) <- sheets
morans <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/morans.csv", sep = "/"), row.names = 1)
closest_coast <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Stations_info.csv", sep = "/")) %>% dplyr::select(all_of(c("id", "Closest_coast")))
#aems <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/aems.csv", sep = "/"))
#aems_season_year <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/aems_season_year.csv", sep = "/"))


#read json file 
params <- fromJSON(paste(HOME_ , "/ISPRA_20152017_Analysis/params.json", sep = "/"))
seasons <- params[["seasons"]]
seasons <- rep(names(seasons), each = 3)

#create season column 
env_data$Date <- as.Date(env_data$Date, format = "%Y-%m-%d")
env_data <- env_data %>% mutate(Season = factor(seasons[as.integer(format(Date, "%m"))]), Year = factor(format(Date, "%Y"))) %>% mutate(Season_year = paste(Season, Year, sep = "_"))
env_data <- merge(env_data, closest_coast, by = "id")


for (sheet in sheets){
    print(sheet)
    print(indval_list[[sheet]] %>% filter(stat > 0.5) %>% group_by(index) %>% summarise(n = n()))
}

boxcox_transform <- function(values, min_val = 0.001) {
  
  if (any(is.na(values))) {
    clean_values <- as.numeric(na.omit(values))
  } else {
     clean_values <- as.numeric(values)
  }
  if (min(clean_values) == 0.0) {
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



filter_merge_data <- function(env_data, sites_taxa, morans, covariates, species, n_moran_vec, threshold_obs_species = 0) {
  temp_env <- env_data %>% dplyr::select(all_of(c("Date", "id", "Year", "Season", "Season_year", covariates))) %>% na.omit()
  temp_taxa <- sites_taxa %>% dplyr::select(all_of(c("Date", "id", species)))
  temp_taxa <- temp_taxa %>% filter(rowSums(temp_taxa[, sapply(temp_taxa, is.numeric)] > 0) > threshold_obs_species) 
  temp_morans <- morans[-1, ] 

  mem_vec <- sapply(C(1:n_moran_vec), function(ith) {paste("MEM", ith, sep = "")})

  merged_data <- merge(
    merge(
      temp_env,
      temp_taxa,
      by = c("Date", "id")
    ), 
    morans[-1, ] %>% dplyr::select(all_of(mem_vec)) %>% rownames_to_column(var = "id"), 
    by = "id"
    )

  return(merged_data)
}

species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.5)%>% rownames()}, simplify = TRUE)))

species_obs <- sites_taxa %>% dplyr::select(where(is.numeric)) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
species <- names(species_obs[species_obs > 170])




### POSSIBLE CHOICES ###
# 1) all covariates
# 2) physical + pH + Chla 
# 3) physical + pH + Chla + TN + SiO4

covariates <- c(
  "T",
  "Salinity",
  "O_sat",
  "pH",
  "Chla",
  "NO3",
  "NO2",
  "NH4",
  "TN",
  "PO4",
  "TP",
  "SiO4", 
  "Closest_coast"
)

all_nutrients <- c(
  "Chla",
  "NO3",
  "NO2",
  "NH4",
  "TN",
  "PO4",
  "TP",
  "SiO4"
)

names(env_data) 
covariates_cop <- c(
  "temp",
  "sal",
  "o2",
  "ph",
  "chl",
  "no3",
  "nppv",
  "mld",
  "phyc",
  "dic",
  "nh4",
  "po4",
  #"TP",
  "Closest_coast"
)


heatmap_data <- env_data %>% dplyr::select(c(all_of(covariates_cop))) %>% cor() %>% as.data.frame() %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to = "Var2", values_to = "value")


pca <- rda(env_data %>% dplyr::select(all_of(covariates_cop)), scale = TRUE)

env <- pca$CA$eig
cumsum(env / sum(env) * 100)
biplot(pca, scaling = 2)
ggplot(data = heatmap_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Variable 1", y = "Variable 2", title = "Correlation Heatmap") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  theme_minimal()


#count rows where at least one covariaets is NA
nrow(env_data) - env_data %>% is.na(dplyr::select(., all_of(covariates))) %>% rowSums() %>% sum()

sum(env_data %>% dplyr::select(all_of(covariates)) %>% is.na() %>% rowSums() > 0 ) / nrow(env_data)

env_data %>% dplyr::select(all_of(covariates)) %>% head(30)

variable_selection <- function(resp, expl) {
  full_model <- rda(resp ~ ., expl)
  print(attr(full_model$terms, "term.labels"))
  print("Forward selection...")
  mod0 <- rda(resp ~ 1, expl)
  forward <- ordiR2step(mod0, full_model, perm.max = 200, trace = FALSE)
  print(attr(forward$terms, "term.labels"))

  print("Backward selection...")
  backward <- ordistep(full_model, direction = "backward", perm.max = 200, trace = FALSE)
  print(attr(backward$terms, "term.labels"))
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

create_resp_expl_data <- function(
  env_data,
  sites_taxa,
  morans,
  covariates, 
  species,
  n_moran_vec, 
  temporal_factor, 
  boxcox = FALSE, 
  rem_multiv_out = FALSE, 
  threshold_obs_species = 1, 
  log_transform = FALSE
 ) {
  print(paste("using as covariates: ", paste(covariates, collapse = ", ")))
  print(paste("number of spatial eigenvectors: ", n_moran_vec))
  print(paste("using as temporal factor: ", temporal_factor))

  print(paste("applying boxcox", boxcox))
  if (boxcox) {
    env_data <- env_data %>% mutate(across(-c(id, Date, Season, Year, Season_year, Closest_coast), ~ boxcox_transform(.))) 
  }

  RDA.data <- filter_merge_data(env_data, sites_taxa, morans, covariates, species, n_moran_vec = n_moran_vec, threshold_obs_species = 10)
  print(paste("Dim after filtering: ", paste(RDA.data %>% dim(), collapse = "x")))
 
  print(paste("removing multivariate outliers:", rem_multiv_out))
  if (rem_multiv_out) {
    D2 <- stats::mahalanobis(RDA.data %>% dplyr::select(all_of(covariates)), colMeans(RDA.data %>% dplyr::select(all_of(covariates))), cov(RDA.data %>% dplyr::select(all_of(covariates))))
    cut_off <- qchisq(0.95, df = length(covariates))
    RDA.data <- RDA.data %>% filter(D2 < cut_off) 
    print(paste("Dim after removing multivariate outliers: ", paste(RDA.data %>% dim(), collapse = "x")))
  }

  mem_vec <- sapply(C(1:n_moran_vec), function(ith) {paste("MEM", ith, sep = "")})

  resp <- RDA.data %>% dplyr::select(all_of(species)) %>% mutate(across(everything(), ~ if (log_transform) log1p(.) else .))  %>% vegan::decostand("hellinger") %>% mutate(across(everything(), ~ . - mean(.)))

  
  expl <- RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(c(covariates, mem_vec, temporal_factor)))
  return(list(resp = resp, expl = expl))
 }




test_models <- function(
  env_data,
  sites_taxa,
  morans,
  covariates, 
  species,
  n_moran_vec, 
  temporal_factor, 
  boxcox = FALSE, 
  rem_multiv_out = FALSE, 
  threshold_obs_species = 10,
  partial = FALSE, 
  partial_covariates = NULL, 
  log_transform = FALSE
  ) {
  
  data <- create_resp_expl_data(
    env_data,
    sites_taxa,
    morans,
    covariates, 
    species,
    n_moran_vec, 
    temporal_factor, 
    boxcox = boxcox, 
    rem_multiv_out = rem_multiv_out, 
    threshold_obs_species = threshold_obs_species, 
    log_transform = log_transform
  )

  if (partial) {
    print(paste("using as partial covariates: ", paste(partial_covariates, collapse = ", ")))
  }
  RDA.model <- compute_rda_model(data$resp, data$expl, partial = partial, partial_covariates = partial_covariates)
  tot.chi  <- RDA.model$CCA$tot.chi + RDA.model$CA$tot.chi
  print(paste("constrained :" , RDA.model$CCA$tot.chi / tot.chi, seq = ""))
  print(paste("unconstrained :" , RDA.model$CA$tot.chi / tot.chi, seq = ""))

  return(RDA.model)
}

ggplot_rda <- function(mod, scaling = 2, fill_row = NULL) {
  scores.sites <- scores(mod, display = "sites", choices = c(1, 2), scaling = scaling)
  scores.species <- scores(mod, display = "species", choices = c(1, 2), scaling = scaling)
  scores.covariates <- scores(mod, display = "bp", choices = c(1, 2), scaling = scaling)

  p <- ggplot() + 
  geom_point(
    data = scores.sites,
    aes(x = RDA1, y = RDA2, fill = fill_row), shape = 21, size = 3
    ) +
    geom_segment(
      data = scores.species,
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2), colour = "black"
    ) +
    geom_segment(
      data = scores.covariates,
      aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = arrow(length = unit(0.2, "cm")), colour = "blue"
    ) +
    geom_text(
      data = scores.covariates,
      aes(x = RDA1, y = RDA2, label = rownames(scores.covariates)),
      hjust = 1.5, vjust = 1.5, colour = "blue"
    ) +
    theme_minimal() + 
    theme(legend.position = "none")
  return(p)
}





##VARIABLE SELECTION
n_moran_vec <- 23 
mem_vec <- sapply(C(1:n_moran_vec), function(ith) {paste("MEM", ith, sep = "")})
species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.5)%>% rownames()}, simplify = TRUE)))
RDA.data <- filter_merge_data(env_data, sites_taxa, morans, covariates_cop, species, n_moran_vec = n_moran_vec, threshold_obs_species = 10)
variable_selection(
  RDA.data %>% dplyr::select(all_of(species)) %>% vegan::decostand("hellinger") %>% mutate(across(all_of(species), ~ . - mean(.))), 
  RDA.data %>% dplyr::mutate(across(all_of(covariates_cop), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(covariates_cop))
  )
variable_selection(
  RDA.data %>% dplyr::select(all_of(species)) %>% vegan::decostand("hellinger") %>% mutate(across(all_of(species), ~ . - mean(.))), 
  RDA.data %>% dplyr::select(all_of(mem_vec))
  )


species_obs <- sites_taxa %>% dplyr::select(where(is.numeric) & !matches("Other.phytoplankton")) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
quantile(species_obs, probs = seq(0, 1, 0.1))
species <- names(species_obs[species_obs > 10 & species_obs < 211])
species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.65)%>% rownames()}, simplify = TRUE)))
length(species)


#BUILDING THE FINAL MODEL
data <- create_resp_expl_data(env_data, sites_taxa, morans, covariates_cop, species, n_moran_vec = 23, temporal_factor = "Season", boxcox = FALSE, rem_multiv_out = FALSE, threshold_obs_species = 5, log_transform = FALSE)
RDA.model <- compute_rda_model(data$resp, data$expl, covariates = covariates_cop, partial = TRUE, partial_covariates = c(mem_vec, "Season"))

RDA.model
vegan::RsquareAdj(RDA.model)
anova.cca(RDA.model, by = "margin")


p <- ggplot_rda(RDA.model, scaling = 2, fill_row = as.factor(data$expl$Season))

p <- ggplot() + 
  geom_point(
    data = scores(RDA.model, display = "sites", choices = c(1, 2), scaling = 2),
    aes(x = RDA1, y = RDA2, fill = as.factor(data$expl$Season)) , shape = 21, size = 3
    )
plot(p)


scores_sp <- scores(RDA.model, display = "species", scaling = 2)
#compute norm rowwise
scores_sp <- scores_sp %>% as.data.frame() %>% mutate(norm = sqrt(RDA1^2 + RDA2^2))
scores_sp[, "norm"]
scores_sp %>% arrange(desc(norm))
barplot(scores_sp %>% arrange(desc(norm)) %>% pull(norm))

#select first 4 rownames 


resid <- cbind(residuals(RDA.model), 
data$expl)

##ggplot scateerplot first three columns againts T 
high_score_species <- scores_sp %>% arrange(desc(norm)) %>% rownames() %>% head(7)
resid_long <- resid %>% dplyr::select(all_of(c("Chla", high_score_species))) %>% 
tidyr::pivot_longer(cols = high_score_species, names_to = "Variable", values_to = "Value")


ggplot(resid_long, aes(x = Chla, y = Value, color = Variable)) +
  geom_point() +
  labs(x = "Chla", y = "Value", color = "Variable") +
  theme_minimal()


only_mem <- rda(data$resp ~ ., data$expl %>% dplyr::select(all_of(mem_vec)))
mem_resid <- residuals(only_mem)
resid[,c(1:3)] %>% head()
mem_resid[,c(1:3)]  %>% head()
plot(data$expl[["T"]], resid[,1])
points(data$expl[["T"]], resid[,2])
points(data$expl[["T"]], resid[,3])

data$expl[["T"]] %>% sd()



for (season in names(params$"seasons")) {
  print(season)
  data <- create_resp_expl_data(
    env_data %>% filter(Season == season), 
    sites_taxa, 
    morans, 
    covariates, 
    species, 
    n_moran_vec = 23, 
    temporal_factor = "Season", 
    boxcox = FALSE, 
    rem_multiv_out = FALSE, 
    threshold_obs_species = 10
    )

  RDA.model <- compute_rda_model(data$resp, data$expl , covariates = covariates, partial = TRUE, partial_covariates = mem_vec)
  
  p <- ggplot_rda(RDA.model, scaling = 2)
  plot(p)
}


names(sites_taxa)

env_data %>% head()
species_obs <- sites_taxa %>% dplyr::select(where(is.numeric)) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
species <- names(species_obs[species_obs > 120])


variation_partitioning <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Season", "Env", "Shared"))))
for (region in unique(env_data$Region)) {
  print(region)
  RDA.data <- filter_merge_data(env_data %>% dplyr::filter(Region == region), sites_taxa, morans, covariates, species, n_moran_vec = n_moran_vec, threshold_obs_species = 5)
  print(paste("Dim. of set", paste(RDA.data %>% dim(), collapse = "x")))
  if (nrow(RDA.data) < 30) {
    next
  }
  varp <- vegan::varpart(
  RDA.data %>% dplyr::select(all_of(species)), 
  RDA.data %>% dplyr::select(Season), 
  RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(covariates)),
  transfo = "hel"
  )
  var_expl <- setNames(data.frame(t(varp$part$indfract[c(1,2,3), "Adj.R.squared"]), row.names = region), c("Season", "Env", "Shared"))
  variation_partitioning <- rbind(var_expl, variation_partitioning)
}


ggplot(variation_partitioning, aes(x = rownames(variation_partitioning),)) +
  geom_col(aes(y = Env + Season + Shared, fill = "green")) +
  geom_col(aes(y = Env + Season, fill = "red")) +
  geom_col(aes(y = Env, fill = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


plot(varp)
summary(varp)

plot(varp, title = "AAA")


as.data.frame(varp$indfract)


varp$part$frac %>% data.frame()

vpar <- varpart(
  data$resp, 
  data$expl %>% dplyr::select(Season), 
  data$expl %>% dplyr::select(all_of(covariates)),
  data$expl %>% dplyr::select(all_of(mem_vec))
)

plot(
  vpar, 
  Xnames = c("S", "C", "M"), 

)



#dbRDa 
species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.5)%>% rownames()}, simplify = TRUE)))

pecies_obs <- sites_taxa %>% dplyr::select(where(is.numeric) & !matches("Other.phytoplankton")) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
species <- names(species_obs[species_obs > 120])
RDA.data <- filter_merge_data(env_data, sites_taxa, morans, covariates_cop, species, n_moran_vec = 23, threshold_obs_species = 10)  
mem_vec <- sapply(C(1:n_moran_vec), function(ith) {paste("MEM", ith, sep = "")})

resp <- RDA.data %>% dplyr::select(all_of(species)) %>% vegan::vegdist("bray", binary = FALSE)  
expl <- RDA.data %>% dplyr::mutate(across(all_of(covariates_cop), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(c(covariates_cop, mem_vec, "Season")))
partial_covariates <- c(mem_vec, "Season")
RDA.model <- vegan::dbrda(formula(
        paste(
          "resp ~", paste(covariates_cop, collapse = " + "), 
          "+ Condition(", 
          paste(partial_covariates, collapse = " + "),
          ")",
          sep = " "
        )
      ), 
    expl
    )

RDA.model
