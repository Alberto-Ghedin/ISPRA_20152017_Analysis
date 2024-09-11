library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(readxl)
library(tibble)
library(jsonlite)
library(zoo)
library(klaR)
library(vegan)
library(jsonlite)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) 
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/sites_taxa_matrix.csv", sep = "/"))
sheets <- getSheetNames(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"))
indval_list <-  lapply(sheets, function(sheet) read.xlsx(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/IndVal_method_2.xlsx", sep = "/"), sheet = sheet, rowNames = TRUE))
morans <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/morans.csv", sep = "/"), row.names = 1)


barplot(morans[1, ])
#read json file 
params <- fromJSON(paste(HOME_ , "/ISPRA_20152017_Analysis/params.json", sep = "/"))
seasons <- params[["seasons"]]
seasons <- rep(names(seasons), each = 3)

#create season column 
env_data$Date <- as.Date(env_data$Date, format = "%Y-%m-%d")
env_data <- env_data %>% mutate(Season = factor(seasons[as.integer(format(Date, "%m"))]), Year = factor(format(Date, "%Y"))) 


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

  lambda <- boxcox(clean_values ~ 1, plotit = FALSE)
  max_lambda <- lambda$x[which(lambda$y == max(lambda$y))]
  if (max_lambda == 0) {
    clean_values <- log(clean_values)
  } else {
    clean_values <- (clean_values^max_lambda - 1) / max_lambda
  }
  values[which(!is.na(values))] <- clean_values
  return(values)
}


filter_merge_data <- function(env_data, sites_taxa, covariates, species, threshold_obs_species = 0) {
  temp_env <- env_data %>% dplyr::select(all_of(c("Date", "id", "Season", "Year", covariates))) %>% na.omit()
  temp_taxa <- sites_taxa %>% dplyr::select(all_of(c("Date", "id", species)))
  temp_taxa <- temp_taxa %>% filter(rowSums(temp_taxa[, sapply(temp_taxa, is.numeric)] > 0) > threshold_obs_species) 

  merged_data <- merge(
    temp_env,
    temp_taxa,
    by = c("Date", "id")
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
  "SiO4"
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

test_models <- function(covariates, 
species,
boxcox = FALSE, 
rem_multiv_out = FALSE, 
threshold_obs_species = 10, 
env_data. = env_data,
sites_taxa. = sites_taxa,
partial = FALSE
) {
  print(paste("using as covariates: ", paste(covariates, collapse = ", ")))

  print(paste("applying boxcox", boxcox))
  if (boxcox) {
    env_data. <- env_data. %>% mutate(across(-c(Region, id, Date, Secchi_depth, Season, Year), ~ boxcox_transform(.))) 
  }

  RDA.data <- filter_merge_data(env_data., sites_taxa., covariates, species, threshold_obs_species = 10)
  print(paste("Dim after filtering: ", paste(RDA.data %>% dim(), collapse = "x")))
 
  print(paste("removing multivariate outliers:", rem_multiv_out))
  if (rem_multiv_out) {
    D2 <- stats::mahalanobis(RDA.data %>% dplyr::select(all_of(covariates)), colMeans(RDA.data %>% dplyr::select(all_of(covariates))), cov(RDA.data %>% dplyr::select(all_of(covariates))))
    cut_off <- qchisq(0.95, df = length(covariates))
    RDA.data <- RDA.data %>% filter(D2 < cut_off) 
    print(paste("Dim after removing multivariate outliers: ", paste(RDA.data %>% dim(), collapse = "x")))
  }

  print("model with no log transformation")
  resp <- RDA.data %>% dplyr::select(all_of(species)) %>% vegan::decostand("hellinger") %>% mutate(across(all_of(species), ~ . - mean(.)))
  expl <- RDA.data %>% dplyr::select(all_of(covariates)) %>% decostand("standardize")
  RDA.model <- rda(resp ~ ., expl)
  tot.chi  <- RDA.model$CCA$tot.chi + RDA.model$CA$tot.chi
  print(paste("constrained :" , RDA.model$CCA$tot.chi / tot.chi, seq = ""))
  print(paste("unconstrained :" , RDA.model$CA$tot.chi / tot.chi, seq = ""))

  print("model with log transformation")
  resp <- RDA.data %>% dplyr::select(all_of(species)) %>% apply(2, function(x) {log(x+1)}) %>% vegan::decostand("hellinger") %>% data.frame() %>% mutate(across(all_of(species), ~ . - mean(.)))
  expl <- RDA.data %>% dplyr::select(all_of(covariates)) %>% decostand("standardize")
  RDA.model <- rda(resp ~ ., expl)
  tot.chi  <- RDA.model$CCA$tot.chi + RDA.model$CA$tot.chi
  print(paste("constrained :" , RDA.model$CCA$tot.chi / tot.chi, seq = ""))
  print(paste("unconstrained :" , RDA.model$CA$tot.chi / tot.chi, seq = ""))
}


species <- unique(unlist(sapply(sheets, function(sheet) {indval_list[[sheet]] %>% filter(stat > 0.5)%>% rownames()}, simplify = TRUE)))
test_models(covariates, species, boxcox = TRUE, rem_multiv_out = TRUE)

species_obs <- sites_taxa %>% dplyr::select(where(is.numeric)) %>% mutate(across(where(is.numeric), ~ .>0)) %>% colSums()
species <- names(species_obs[species_obs > 170])

test_models(covariates, species, boxcox = TRUE, rem_multiv_out = TRUE)



mod0 <- rda(resp ~ 1, expl)
### forward selection: Salinity + T + TN + SiO4 + O_sat + PO4 + Chla + NH4 +      pH + NO2 + TP
ordiR2step(mod0, RDA.model, perm.max = 200)
### backward selection:  T + Salinity + O_sat + pH + Chla + NO3 + NO2 + NH4 + TN + PO4 + TP + SiO4
ordistep(RDA.model, direction = "backward", perm.max = 200)






mod0 <- rda(resp ~ 1, expl)
### forward selection: TN + T + PO4 + Salinity + O_sat + SiO4 + TP + Chla + pH + NO2 + NH4 + NO3
ordiR2step(mod0, RDA.model, perm.max = 200)
### backward selection:  T + Salinity + O_sat + pH + Chla + NO3 + NO2 + NH4 + TN + PO4 + TP + SiO4
ordistep(RDA.model, direction = "backward", perm.max = 200)





## PARTIAL RDA
resp <- RDA.data %>% dplyr::select(all_of(species)) %>% vegan::decostand("hellinger") %>% mutate(across(all_of(species), ~ . - mean(.)))
expl <- RDA.data %>% dplyr::select(all_of(c(covariates, "Year", "Season"))) %>% dplyr::mutate(across(all_of(covariates), ~  decostand(.,"standardize"))) 
RDA.pmodel <- rda(formula(paste("resp ~", paste(covariates, collapse = " + "), "+ Condition(Season)"), sep = " "), expl)
summary(RDA.pmodel)
RsquareAdj(RDA.pmodel)
anova.cca(RDA.model, steps = 1000)


resp <- RDA.data%>% dplyr::select(all_of(species)) %>% apply(2, function(x) {log(x+1)}) %>% vegan::decostand("hellinger") %>% data.frame() %>% mutate(across(all_of(species), ~ . - mean(.)))
expl <- RDA.data %>% dplyr::select(all_of(c(covariates, "Year", "Season"))) %>% dplyr::mutate(across(all_of(covariates), ~  decostand(.,"standardize"))) 
RDA.pmodel <- rda(formula(paste("resp ~", paste(covariates, collapse = " + "), "+ Condition(Season)"), sep = " "), expl)
RsquareAdj(RDA.pmodel)
summary(RDA.pmodel)


summary(RDA.model)
names(RDA.model)

summary.rda

RDA.model$tot.chi
