library(MASS)
library(ggplot2)
library(dplyr)
library(lubridate)
library(openxlsx)
library(readxl)
library(tibble)
library(jsonlite)
library(zoo)
library(klaR)
library(DAAG)
library(rstatix)
library(ggpubr)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

training_data <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/training_data_LDA.csv", sep = "/")) 

clustering <- "ward_4"

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

training_data %>% group_by(ward_4) %>% summarise(n = n())

total_samples <- training_data %>% group_by(!!as.symbol(clustering)) %>% summarise(n = n())

training_data %>% 
  select(all_of(c(clustering, covariates))) %>%
  na.omit() %>% 
  group_by(!!as.symbol(clustering)) %>%
  summarise(total_rows = n()) / total_samples

training_data.clean <- training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit()

training_data.clean %>% select(all_of(covariates)) %>% base::sum()
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


model <- lm(formula(paste0("cbind(",paste0(covariates, collapse = ","), ") ~ ", clustering, sep = "")), training_data.clean )
Manova(model, test.statistic = "Pillai") 


## LDA ## 
lda_formula <- as.formula(
  paste(clustering, "~", paste(covariates, collapse = " + "))
)
LDA <- MASS::lda(lda_formula, data = training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit(), CV = FALSE)
confusion(LDA$class, training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit() %>% pull(clustering))

LDA$svd

heatmap(confusion(LDA$class, training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit() %>% pull(clustering))$confusion, Colv=NA, Rowv=NA, scale='none')


qda_formula <- as.formula(
  paste(clustering, "~", paste(covariates, collapse = " + "))
)
QDA <- MASS::qda(qda_formula, data = training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit(), CV = TRUE)


confusion(QDA$class, training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit() %>% pull(clustering))

heatmap(confusion(QDA$class, training_data %>% select(all_of(c(clustering, covariates))) %>% na.omit() %>% pull(clustering))$confusion, Colv=NA, Rowv=NA, scale='none')




df_lda <- data.frame(groups = na.omit(training_data %>% select(all_of(c(clustering, covariates)))) %>% pull(clustering), lda = LDA$x)
lda_scaling <- data.frame(LD1 = 5 * LDA$scaling[,1], LD2 = 5 * LDA$scaling[,2])
centroids <- df_lda %>%
  group_by(groups) %>%
  summarize(centroid_LD1 = mean(lda.LD1), centroid_LD2 = mean(lda.LD2))

p <- ggplot() +
    geom_point(data = df_lda,
                         aes(x = lda.LD1, 
                                 y = lda.LD2, 
                                 col = factor(groups)), 
                         size = 4) +
    labs(color = "Groups") +
    scale_color_manual(values = c("blue", "red", "green", "orange")) + 
    geom_segment(data = lda_scaling, 
                             aes(x = 0, y = 0, xend = LD1, yend = LD2), 
                             arrow = arrow(length = unit(0.3, "cm")), 
                             color = "black") + 
    geom_point(data = centroids, 
                         aes(x = centroid_LD1, y = centroid_LD2, col = factor(groups)), 
                         size = 6, shape = 4) +
    theme_classic() + # formatting the plot to make it pretty
    theme(axis.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 18))



confusion(predict(LDA))


table(predict(LDA)$class, training_data$ward_4)

predict(LDA)$class
training_data$ward_4
p <- klaR::partimat(as.factor(ward_4) ~ T + Salinity + O_sat + pH + Chla + NO3 + NO2 + NH4 + TN + PO4 + TP + SiO4, data = training_data, method = "lda")

warnings()
