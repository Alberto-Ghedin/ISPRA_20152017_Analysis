library(dplyr)
library(adespatial)
library(lubridate)

HOME_ <- paste(path.expand("~"), "PHD", sep = "/")

env_data <- read.csv(paste(HOME_, "/ISPRA_20152017_Analysis/Create_dataset/df_chem_phys_mod_data_cleaned_long_format.csv", sep = "/")) 
sites_taxa <- read.csv(paste(HOME_, "ISPRA_20152017_Analysis/sites_taxa_matrix.csv", sep = "/"))

env_data$Date <- as.Date(env_data$Date, format = "%Y-%m-%d")
sites_taxa$Date <- as.Date(sites_taxa$Date, format = "%Y-%m-%d")

common_dates <- sort(as.Date(intersect(env_data$Date, sites_taxa$Date)))


# Function to compute the difference in months between two dates
diff_in_months <- function(date1, date2) {
  time_length(interval(date2, date1), "months")
}

diff_in_months(common_dates[3], common_dates[1])
# Create a matrix where each entry is the difference in months between dates
date_diff_matrix <- outer(common_dates, common_dates, diff_in_months)

corr_threshold <- 6

binary.mat <- ifelse(date_diff_matrix > 0 & date_diff_matrix < corr_threshold, 1, 0)

max(date_diff_matrix)
date_diff_matrix[c(1:6), c(1:6)]
binary.mat[c(1:7), c(1:7)]


indices <- which(binary.mat == 1, arr.ind = TRUE)
n <- length(common_dates)
years <- as.numeric(format(common_dates, format = "%Y"))
months <- as.numeric(format(common_dates, format = "%m"))

aem_bn <- aem.build.binary(coords=cbind(c(1:n), years - min(years), months),link=indices[, c(2,1)], rm.same.y = FALSE, rot.angle = 90)
weight <- numeric()

date_diff_matrix
for(i in 1:length(indices) / 2){
weight[i] <- 1 - (date_diff_matrix[indices[i,1],indices[i,2]] / corr_threshold) ^2
}


aem_res <- aem(aem_bn, weight = weight, rm.link0 = TRUE, print.binary.mat = TRUE)

df_aem <- data.frame(Date = c("eigen", common_dates), rbind(aem_res$values, aem_res$vectors))

names(df_aem) <- c("Date", sapply(c(1:(length(common_dates) -1)), function(ith) {paste("AEM", ith, sep = "")}))

write.csv(df_aem, paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/aems.csv", sep = "/"), row.names = FALSE)




#AEM for season-year
aems <- aem.time(10)

Seasons <- c(c("Summer", "Autumn"), rep(c("Winter", "Spring", "Summer", "Autumn"), 2))
Seasons
years <- c(rep(2015, 2), rep(2016, 4), rep(2017, 4))
df_aem <- data.frame(Season = c("eigen", Seasons), Year = c("eigen", years), rbind(aems$values, aems$aem))

names(df_aem) <- c("Season", "Year", sapply(c(1:(length(Seasons) -1)), function(ith) {paste("AEM", ith, sep = "")}))

write.csv(df_aem, paste(HOME_, "ISPRA_20152017_Analysis/Clustering/Unknown_effect/aems_season_year.csv", sep = "/"), row.names = FALSE)
