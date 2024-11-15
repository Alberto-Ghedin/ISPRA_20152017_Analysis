library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS) 
library(lme4)
library(glmmTMB)

HOME_ <- "."
phyto_abund <- read.csv(paste(HOME_, "phyto_abund.csv", sep = "/"))

sample_abund <- phyto_abund %>% group_by(Date, id) %>% summarise(sample_abund = as.integer(sum(Num_cell_l)), Region = first(Region), Season = first(Season), Basin = first(Basin)) 
phyto_abund %>% head()
model_nb <-  MASS::glm.nb(sample_abund ~ Region, data = sample_abund)
model_lmer <- lme4::lmer(sample_abund ~ Region + (1| Season), data = sample_abund)
model_nbmx <- lme4::glmer.nb(sample_abund ~ Season + (1| Region), data = sample_abund)

summary(model_nb)
summary(model_lmer)
summary(model_nbmx)
qqnorm(resid(model_nb))
qqline(resid(model_nb))
model_nb
res <- residuals(model_nb)
fitted <- fitted(model_nb)

plot(log(fitted), res, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, col = "red")

dev.off()
plot(model_lmer)

qqnorm(res)
abline(a = 0, b = 1, col = "red")
taxon_abund <- phyto_abund %>% group_by(Taxon) %>% summarise(mean = mean(Num_cell_l), sd = sd(Num_cell_l))


ggplot(taxon_abund) + geom_point(aes(x = log(mean), y = log(sd))) + theme_minimal() + labs(x = "Mean", y = "SD")
