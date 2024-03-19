library(iNEXT)
library(dplyr)
library(tidyr)
library(vegan)
library(tibble)
file_path <- "/mnt/d/PHD/ISPRA_20152017_Analysis/Tests/Olli_2019"

samp <- read.table(file = paste(file_path, '/sample.txt', sep = ""), header = T, sep = '\t', stringsAsFactors = F)
samp$date <- as.Date(samp$date)

# parse Date
samp$year <- as.numeric(format(samp$date, '%Y'))
samp$month <- as.numeric(format(samp$date, '%m'))
samp$jul <- as.numeric(format(samp$date, '%j'))


# load count table
count <- read.table(file = paste(file_path, '/count.txt', sep = ""), header = T, sep = '\t', stringsAsFactors = F)



# make community matrix
xtab <- tapply(count$biomass, list(as.factor(count$sampleID), as.factor(count$name)), sum)

xtab <- count %>% 
    group_by(sampleID, name) %>% 
    summarise(biomass = sum(biomass)) %>% 
    pivot_wider(names_from = name, values_from = biomass, values_fill = 0) %>%
    column_to_rownames(var = "sampleID")

dim(xtab) # should be: 15640 samples x 1383 taxa

xtab[is.na(xtab)] <- 0 # change all NA's to zero

xtab %>% filter(sampleID %in% samp$sampleID) %>% as.matrix()

# add taxon richness to sample table
samp$nsp <- rowSums(decostand(xtab[samp$sampleID, ], 'pa'))

# separate sample tables for the Baltic Seaa (bs) and Chesapeake Bay (chesa)

samp.bs <- subset(samp, basin == 'bs') # 7731  samples
samp.chesa <- subset(samp, basin == 'chesa') # 7896 samples

taxa <-  read.table(file = paste(file_path, '/taxa.txt', sep = ""), header = T, sep = '\t', stringsAsFactors = F)


# salinity windows
s.bs <- seq(1, 26, by = 0.5) # salinity windows for the Baltic Sea
s.chesa <- seq(1, 27, by = 0.5) # salinity windows for the Chesapeake Bay

# frequency distribution of samples within salinity windows
x.bs <- x.chesa <- numeric()

for(i in 1:length(s.bs)){ x.bs[i] <- length(which(samp.bs$Salin > (s.bs[i]-1) & samp.bs$Salin < s.bs[i]+1))}
for(i in 1:length(s.chesa)){ x.chesa[i] <- length(which(samp.chesa$Salin > s.chesa[i]-1 & samp.chesa$Salin < s.chesa[i]+1))}

minx <- min(c(x.bs, x.chesa))

for (i in 1:length(s.bs)) {lst.bs[[i]] <- t(xtab[samp.bs[which(samp.bs$Salin > s.bs[i]-1 & samp.bs$Salin < s.bs[i]+1), 'sampleID'], ])}
names(lst.bs) <- paste('Sal', s.bs, sep = '')

lst.bs[[1]][c(1:2),c(1:4)]
