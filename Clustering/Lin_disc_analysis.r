library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(labdsv)
library(openxlsx)
library(readxl)
library(tibble)
library(jsonlite)
library(zoo)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

training_data <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/training_data_LDA.csv", sep = "/")) 

names(training_data)

LDA <- MASS::lda(ward_4 ~ T + Salinity + O_sat + pH + Chla + NO3 + NO2 + NH4 + TN + PO4 + TP + SiO4, data = training_data)

df_lda <- data.frame(groups = training_data$ward_4, lda = predict(LDA)$x)

lda_scaling <- data.frame(LD1 = 5 * LDA$scaling[,1], LD2 = 5 * LDA$scaling[,2])

centroids <- df_lda %>%
  group_by(groups) %>%
  summarize(centroid_LD1 = mean(lda.LD1), centroid_LD2 = mean(lda.LD2))

ggplot() +
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




LDA$scaling

library(klaR)

klaR::partimat(as.factor(ward_4) ~ T + Salinity + O_sat + pH + Chla + NO3 + NO2 + NH4 + TN + PO4 + TP + SiO4, data = training_data, method = "lda")
