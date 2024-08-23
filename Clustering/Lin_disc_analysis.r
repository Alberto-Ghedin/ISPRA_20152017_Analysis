library(MASS)
library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(tibble)
library(jsonlite)
library(zoo)
library(klaR)
library(DAAG)
library(rstatix)
library(ggpubr)
library(vegan)
library(tidyr)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) 
cluster_indices <- read.xlsx(paste(HOME_, "/ISPRA_20152017_Analysis/Clustering/Unknown_effect/cluster_indices.xlsx", sep = ""),sheet = "method_2", detectDates = TRUE)


cluster_indices$Date <- openxlsx::convertToDate(cluster_indices$Date)
env_data$Date <- as.Date(env_data$Date)
training_data <- merge(env_data, cluster_indices %>% select(Date, id, ward_4, ward_11), by = c("Date", "id"))

clustering <- "ward_11"

### POSSIBLE CHOICES ###
# 1) all covariates
# 2) physical + pH + Chla 
# 3) physical + pH + Chla + TN + SiO4

#c("Chla", "pH", "NH4", "TN", "TP")
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


training_data %>% dplyr::select(all_of(c(clustering, covariates))) %>% summary()


training_data %>% 
  select(all_of(c(clustering, covariates))) %>%
  na.omit() %>% 
  group_by(!!as.symbol(clustering)) %>%
  summarise(total_rows = n()) / total_samples


# Assuming training_data is your dataframe and covariates is a vector of column names
# Reshape the data to long format
long_data <- training_data %>%
  select(all_of(covariates)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Create histograms for each variable
ggplot(long_data, aes(x = value)) +
  geom_histogram(binwidth = 30, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Histograms of Covariates", x = "Value", y = "Frequency")




training_data.clean <- training_data %>%
  select(all_of(c(clustering, covariates))) %>%
  na.omit() %>%
  mutate(across(all_of(covariates), ~ decostand(boxcox_transform(.), "standardize"))) 


ggplot(training_data.clean %>%
  select(all_of(covariates)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value"), aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Histograms of Covariates", x = "Value", y = "Frequency")

## MANOVA TEST ##

# testing for normality of each covariates within each group
training_data %>% 
  group_by(!!as.symbol(clustering)) %>% 
  shapiro_test(covariates) %>% arrange(variable) %>% filter(p > 0.01)

ggqqplot(training_data.clean, "Chla", facet.by = clustering,
         ylab = "Chla", ggtheme = theme_bw())
# No multinormality

#testing for homogeneity
box_m(training_data.clean[, covariates] %>% as.matrix(), training_data.clean[[clustering]] %>% as.matrix())

training_data.clean %>% 
  gather(key = "variable", value = "value", all_of(covariates)) %>% 
  group_by(variable) %>%
  levene_test(as.formula(paste("value ~ as.factor(", clustering, ")", sep = "")))

# No homogeneity



## LOOKING FOR STATISTICALLY DIFERENT VARIABLES ##
res.man <- manova(as.formula(paste0("cbind(",paste0(covariates, collapse = ","), ") ~ factor(", clustering, ")", sep = "")), training_data.clean )
summary(res.man)

post_hoc <- summary.aov(res.man)
for (i in seq_along(covariates)) {
  print(paste("p-val of", covariates[i], "is", post_hoc[[i]][["Pr(>F)"]][[1]]))
}


summary(manova(as.formula(paste0("cbind(",paste0(covariates, collapse = ","), ") ~ ", clustering, sep = "")), training_data.clean ))
summary(manova(as.formula(paste0("cbind(",paste0(covariates_wrong, collapse = ","), ") ~ ", clustering, sep = "")), training_data.clean ))

for (cov in covariates_wrong) {
  print(paste("p-val of", cov, "is", summary(manova(as.formula(paste0(cov, " ~ ", clustering, sep = "")), training_data.clean ))$stats[1, 4]))
}

names(aov(O_sat ~ ward_11, data = training_data.clean))

oneway.test(NO2~ ward_11, data = training_data.clean)

summary(aov(O_sat ~ factor(ward_11), data = training_data.clean))

## LDA ## 

covariates <- c("NO3","NO2",  "Salinity", "Chla", "O_sat", "SiO4", "PO4", "T")
lda_formula <- as.formula(
  paste(clustering, "~", paste(covariates, collapse = " + "))
)
LDA <- MASS::lda(lda_formula, data = training_data.clean, CV = FALSE)

sum(LDA$svd^2)

df_lda <- data.frame(training_data.clean, predict(LDA)$x)
lda_scaling <- data.frame(LD1 = LDA$scaling[,1], LD2 = LDA$scaling[,2])
centroids <- df_lda %>%
  group_by(!!sym(clustering)) %>%
  summarize(LD1 = mean(LD1), LD2 = mean(LD2))

explained_var <- LDA$svd^2 / sum(LDA$svd^2)

LDA$scaling %>% apply(1, norm, type = "2") %>% sort(decreasing = TRUE) %>% plot()
LDA$scaling %>% apply(1, norm, type = "2") %>% sort(decreasing = TRUE)


## PLOT WITH CENTROIDS AND ARROWS
p <- ggplot() +
    geom_segment(data = lda_scaling, 
                             aes(x = 0, y = 0, xend = LD1, yend = LD2), 
                             arrow = arrow(length = unit(0.3, "cm")), 
                             color = "black") + 
    geom_text(data = lda_scaling, 
                      aes(x = LD1, y = LD2, label = rownames(lda_scaling)), 
                      hjust = -0.5, vjust = -0.5) +
    geom_text(data = centroids, 
                         aes(x = LD1, y = LD2, label = rownames(centroids), colour = rownames(centroids)), 
                         size = 6, 
              show.legend = FALSE) +
    theme_classic() + 
    labs(
      x = paste("LD1 (", round(explained_var[1] * 100, 2), "%)", sep = ""),
      y = paste("LD2 (", round(explained_var[2] * 100, 2), "%)", sep = "")
    ) +
    theme(axis.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 18))

p

var(training_data.clean %>% select(all_of(covariates)))
df_lda %>% head()
ggplot() +
    geom_point(data = df_lda,
                         aes(x = LD1, 
                                 y = LD2, 
                                 col = factor(!!sym(clustering))), 
                         size = 3) +
    stat_ellipse(data = df_lda, aes(x = LD1, y = LD2, fill = factor(ward_11)), level = 0.95, alpha = 0.4) +
    labs(color = "Groups") +
    geom_segment(data = lda_scaling, 
                             aes(x = 0, y = 0, xend = LD1, yend = LD2), 
                             arrow = arrow(length = unit(0.3, "cm")), 
                             color = "black") + 
    geom_text(data = lda_scaling, 
                      aes(x = LD1, y = LD2, label = rownames(lda_scaling)), 
                      hjust = -0.5, vjust = -0.5) +
    geom_text(data = centroids, 
                         aes(x = LD1, y = LD2, label = rownames(centroids), colour = rownames(centroids)), 
                         size = 6, 
              show.legend = FALSE) +
    theme_classic() + 
    labs(
      x = paste("LD1 (", round(explained_var[1] * 100, 2), "%)", sep = ""),
      y = paste("LD2 (", round(explained_var[2] * 100, 2), "%)", sep = "")
    ) +
    theme(axis.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 18))

pca <- vegan::rda(training_data.clean %>% select(all_of(covariates)), SCALE = FALSE)

biplot(pca, scaling = 1, type = "text")

pca$CA$eig / sum(pca$CA$eig)

vegan::cleanplot.pca(pca)




## F-ratio test 
for (removed_cov in covariates) {
  # Create a new list of covariates excluding the current one
  covariates_subset <- setdiff(covariates, removed_cov)

  # Define the LDA formula with the updated list of covariates
  lda_formula <- as.formula(
    paste(clustering, "~", paste(covariates_subset, collapse = " + "))
  )
  LDA <- MASS::lda(lda_formula, data = training_data.clean %>% select(all_of(c(clustering, covariates))) %>% na.omit(), CV = FALSE)
  print(paste(
    "F-ratio with", removed_cov, "removed is:", 
    sum(LDA$svd^2), "with", length(LDA$svd), "components"
  ),
  sep = " " 
  )
}


backward_elimination <- function(dataset, covariates, model = NULL, covariates_removed = NULL, f_ratios = NULL) {
  if (length(covariates) == 1) {
    return(c(model = list(model), f_ratios = list(setNames(f_ratios, covariates_removed))))
  }

  if (is.null(f_ratios)) {
    lda_formula <- as.formula(
      paste(clustering, "~", paste(covariates, collapse = " + "))
    )
    LDA <- MASS::lda(lda_formula, data = dataset, CV = FALSE)
    f_ratios <- sum(LDA$svd^2)
    covariates_removed <- "full_model"

    backward_elimination(dataset, covariates, covariates_removed = covariates_removed, f_ratios = f_ratios)
  }

  results <- vector("double", length = length(covariates))
  names(results) <- covariates

  for (removed_cov in covariates) {
    # Create a new list of covariates excluding the current one
  covariates_subset <- setdiff(covariates, removed_cov)

  
  # Define the LDA formula with the updated list of covariates
  lda_formula <- as.formula(
    paste(clustering, "~", paste(covariates_subset, collapse = " + "))
  )
  LDA <- MASS::lda(lda_formula, data = dataset, CV = FALSE)
  results[[removed_cov]] <- sum(LDA$svd^2)
  }

  cov_to_remove <- names(results)[which(results == max(results))]
  covariates <- setdiff(covariates, cov_to_remove)
  f_ratios <- c(f_ratios, results[[cov_to_remove]])
  covariates_removed <- c(covariates_removed, cov_to_remove)
  
  backward_elimination(dataset, covariates, model = LDA, covariates_removed = covariates_removed, f_ratios = f_ratios)
}



test <- backward_elimination(training_data.clean %>% select(all_of(c(clustering, covariates))), covariates)


plot(x = test$covariates_removed, y = test$f_ratios,xlab = "Number of covariates removed", ylab = "F-ratio", main = "F-ratio vs Number of Covariates Removed")

test$f_ratios
plot(test$f_ratios)
plot(test$f_ratios[-length(test$f_ratios)] - test$f_ratios[-1], type = "o", xlab = "Number of covariates removed", ylab = "F-ratio", main = "F-ratio vs Number of Covariates Removed")

test$f_ratios

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
covariates <- c("NO2",  "Salinity", "Chla", "O_sat", "SiO4", "PO4", "T", "TN")
lda_formula <- as.formula(
  paste(clustering, "~", paste(covariates, collapse = " + "))
)
LDA <- MASS::lda(lda_formula, data = training_data.clean, CV = FALSE)
sum(LDA$svd^2)
