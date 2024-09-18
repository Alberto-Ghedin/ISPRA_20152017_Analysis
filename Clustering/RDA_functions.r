library(MASS)
library(dplyr)


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
