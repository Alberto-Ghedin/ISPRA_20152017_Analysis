library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(tibble)
library(zoo)
library(vegan)
library(jsonlite)
library(tidyverse)
source("Clustering/RDA_functions.r")

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

#CHOOSING THE ENV_DATA  
env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) %>% dplyr::select(-c(Secchi_depth)) %>% dplyr::filter(!Region %in% c("Lazio", "Basilicata", "Calabria"))
env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Cop_variables_on_sampling_sites.csv", sep = "/")) %>% dplyr::select(-c(Latitude, Longitude, Closest_coast)) 

#LOADING OTHER DATASETS
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
  "SiO4"#, 
  #"Closest_coast"
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
quantile(species_obs, probs = seq(0, 1, 0.05))
sort(species_obs)
species <- names(species_obs[species_obs > 120 & species_obs < 600])
species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.65)%>% rownames()}, simplify = TRUE)))
length(species)


#BUILDING THE FINAL MODEL
data <- create_resp_expl_data(env_data, sites_taxa, morans, covariates, species, n_moran_vec = 23, temporal_factor = "Season", boxcox = FALSE, rem_multiv_out = FALSE, threshold_obs_species = 5, log_transform = FALSE)
RDA.model <- compute_rda_model(data$resp, data$expl, covariates = covariates, partial = TRUE, partial_covariates = c(mem_vec, "Season"))

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


ggplot(variation_partitioning, aes(x = rownames(variation_partitioning),)) +
  geom_col(aes(y = Env + Season + Shared, fill = "green")) +
  geom_col(aes(y = Env + Season, fill = "red")) +
  geom_col(aes(y = Env, fill = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


selected_Regions <- c("Friuli-Venezia Giulia", "Veneto", "Emilia-Romagna", "Marche", "Abruzzo", "Molise", "Puglia")
selected_Regions <- c("Liguria", "Toscana", "Campania", "Sardegna")
data <- create_resp_expl_data(env_data %>% dplyr::filter(Region %in% selected_Regions), sites_taxa, morans, covariates, species, n_moran_vec = 23, temporal_factor = "Season", boxcox = FALSE, rem_multiv_out = FALSE, threshold_obs_species = 5, log_transform = FALSE)
RDA.model <- compute_rda_model(data$resp, data$expl, covariates = covariates, partial = TRUE, partial_covariates = c(mem_vec, "Season"))
RDA.model





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



specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
data

#VARIATION PARTITIONING OF ALL MODELS
adj_r2 <- list()
effect <- list()
model <- list()
#species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.5)%>% rownames()}, simplify = TRUE)))
species_obs <- sites_taxa %>% dplyr::select(where(is.numeric) & !matches("Other.phytoplankton")) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
species <- names(species_obs[species_obs > 170])
n_moran_vec <- 23 
mem_vec <- sapply(C(1:n_moran_vec), function(ith) {paste("MEM", ith, sep = "")})


##FULL MODEL###
RDA.data <- filter_merge_data(env_data, sites_taxa, morans, covariates, species, n_moran_vec = n_moran_vec, threshold_obs_species = 5)

varp <- vegan::varpart(
  RDA.data %>% dplyr::select(all_of(species)), 
  RDA.data %>% dplyr::select(all_of(mem_vec)),
  RDA.data %>% dplyr::select(Season), 
  RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(covariates)),
  transfo = "hel"
)

adj_r2 <- unlist(c(adj_r2,  varp$part$indfract[, "Adj.R.square"]))
effect <- unlist(c(effect, c("MEM", "Season", "Env", "MEM + Season", "MEM + Env", "Season + Env", "MEM + Season + Env", "Residual")))
model <- unlist(c(model, rep("Full", length(varp$part$indfract[, "Adj.R.square"]))))



##MODEL PER SEASON###
for (season in names(params$"seasons")) {
  print(season)
  RDA.data <- filter_merge_data(env_data %>% dplyr::filter(Season == season), sites_taxa, morans, covariates, species, n_moran_vec = n_moran_vec, threshold_obs_species = 5)
  varp <- vegan::varpart(
    RDA.data %>% dplyr::select(all_of(species)), 
    RDA.data %>% dplyr::select(all_of(mem_vec)),
    RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(covariates)),
    transfo = "hel"
  )
  adj_r2 <- unlist(c(adj_r2,  varp$part$indfract[, "Adj.R.squared"]))
  effect <- unlist(c(effect, c("MEM", "Env", "MEM + Env", "Residual")))
  model <- unlist(c(model, rep(season, length(varp$part$indfract[, "Adj.R.squared"]))))
}

##MODEL PER REGION###
regions <- unique(env_data$Region)
for (region in regions) {
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
  adj_r2 <- unlist(c(adj_r2,  varp$part$indfract[, "Adj.R.squared"]))
  effect <- unlist(c(effect, c("Season", "Env", "Season + Env", "Residual")))
  model <- unlist(c(model, rep(region, length(varp$part$indfract[, "Adj.R.squared"]))))
}

##MODEL WITH BRAY CURTIS DISTANCE###
varp <- vegan::varpart(
  RDA.data %>% dplyr::select(all_of(mem_vec)),
  RDA.data %>% dplyr::select(Season), 
  RDA.data %>% dplyr::mutate(across(all_of(covariates), ~ decostand(., "standardize"))) %>% dplyr::select(all_of(covariates))
)

adj_r2 <- unlist(c(adj_r2,  varp$part$indfract[, "Adj.R.squared"]))
effect <- unlist(c(effect, c("Season", "Env", "Season + Env", "Residual")))
model <- unlist(c(model, rep("Bray-Curtis", length(varp$part$indfract[, "Adj.R.squared"]))))


variation <- data.frame(adj_r2 = adj_r2, effect = effect, model = model)

variation <- variation %>% mutate(adj_r2 = ifelse(adj_r2 < 0, 0, adj_r2))

ggplot(variation, aes(fill=effect, y=adj_r2, x=model)) + 
    geom_bar(position="fill", stat="identity")
