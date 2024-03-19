# This R script can be run by sourcing, provided the data tables,sample.txt, count.txt, and taxa.txt, are in the R working directory *and* the required packages are installed. The figures and tables will be saved in the working directory.

# The data is available from: https://datadryad.org/resource/doi:10.5061/dryad.q6m40r0

# The data for the Chesapeake Bay were kindly provided by Kevin G. Sellner (Chesapeake Research Consortium, USA) 
# The Baltic Sea data were provided by eight academic institutions: 
#   Finnish Environment Institute, 
#   Finnish Institute of Marine Research, 
#   City of Helsinki Environmental Centre (Finland), 
#   Institute of Aquatic Sciences (Latvia), Stockholm University (Sweden), 
#   Institut für Ostseeforschung Warnemünde (Germany), 
#   National Environmental Research Institute (Denmark),
#   Estonian Marine Institute.


# source('./SalinityDiversityCode.R')

# The figures and tables were produced with R version 3.4.1 (2017-06-30) on x86_64-apple-darwin15.6.0 platform, using the package versions specified below.
# If using the script on ohter platforms, change the quartz() command to platform specific one. 

# The taxon table code chunk reveals how the taxon table was constructed from the count table, using the WORMS web sevices through the taxiezoop library. It is not needed to run, as the taxon table is provides for convenience.

# Load required packages
library(vegan) # vegan version 2.4-3 
library(mgcv) # mgcv version 1.8-17
library(iNEXT) # iNEXT version 2.0.15 
library(plotrix) # plotrix version 3.6-6
library(colorRamps) # colorRamps version 2.3
library(rworldmap) # rworldmap version 1.3-6
library(rworldxtra) # rworldxtra version 1.01
library(plyr) # plyr version 1.8.4
library(latexpdf) # latexpdf version 0.1.6 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LOAD DATA               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# the data files should be in the working directory.

# load sample table
samp <- read.table(file = './sample.txt', header = T, sep = '\t', stringsAsFactors = F)
samp$date <- as.Date(samp$date)

# parse Date
samp$year <- as.numeric(format(samp$date, '%Y'))
samp$month <- as.numeric(format(samp$date, '%m'))
samp$jul <- as.numeric(format(samp$date, '%j'))

# load count table
count <- read.table(file = './count.txt', header = T, sep = '\t', stringsAsFactors = F)

# make community matrix
xtab <- tapply(count$biomass, list(as.factor(count$sampleID), as.factor(count$name)), sum)
dim(xtab) # should be: 15640 samples x 1383 taxa
xtab[is.na(xtab)] <- 0 # change all NA's to zero

# add taxon richness to sample table
samp$nsp <- rowSums(decostand(xtab[samp$sampleID,], 'pa'))

# separate sample tables for the Baltic Seaa (bs) and Chesapeake Bay (chesa)

samp.bs <- subset(samp, basin == 'bs') # 7731  samples
samp.chesa <- subset(samp, basin == 'chesa') # 7896 samples

taxa <-  read.table(file = './taxa.txt', header = T, sep = '\t', stringsAsFactors = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIG 1 MAP               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stns.chesa <- unique(samp.chesa[, c('lon','lat'),])
stans.bs <- unique(samp.bs[, c('lon','lat'),])

map <- getMap(resolution = 'high')
lmap <- getMap(resolution = 'low')

quartz(width = 8, height = 4)
par(mfrow = c(1,2), mar = c(1,4,1,1), oma = c(2,0,0,0))

plot(lmap, xlim=c(9,30), ylim=c(54,66), col=grey(0.95), lwd=.75, axes=T, las = 1) 
points(stans.bs, cex = 0.7, pch = 19)
mtext(' A', side = 3, line = -1.5, adj = 0, cex = 1.5)

plot(map,xlim=c(-78, -75),ylim=c(36.5, 39.8), col=grey(0.95), lwd=.75, axes=T, las = 1) # chesa
points(stns.chesa, cex = 0.7, pch = 19)
mtext(' B', side = 3, line = -1.5, adj = 0, cex = 1.5)
dev.copy2pdf(file='./Fig1.pdf')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIG 2 SALINITY vs DIVERSITY ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# GAM models
gm.bs <- gam(nsp~s(Salin), data = samp.bs, method = 'REML') 
gm.chesa <- gam(nsp~s(Salin), data = samp.chesa, method = 'REML')

# predicted values
gm.bs$pred <- predict(gm.bs, se = T,  newdata = data.frame(Salin=seq(0, 31, by = 0.1)))
gm.chesa$pred <- predict(gm.chesa, se = T,  newdata = data.frame(Salin=seq(0, 31, by = 0.1)))
gm.bs$pred$Salin <- gm.chesa$pred$Salin <- seq(0, 31, by = 0.1)

quartz('Fig 2', width = 8, height = 4)
par(mfrow = c(1,2), mar = c(4,4,1,1), oma = c(0,0,0,0))

# Baltic Sea panel
plot(nsp~Salin, data = samp.bs, ylab = 'Species richness', xlab='Salinity', ylim = c(0, 70), pch = 19, col = rgb(0,0,0, 0.2), las = 1)
lines(gm.bs$pred$Salin, gm.bs$pred$fit, col = 2, lwd= 3)
lines(gm.bs$pred$Salin, gm.bs$pred$fit + gm.bs$pred$se, col=2, lty=2)
lines(gm.bs$pred$Salin, gm.bs$pred$fit - gm.bs$pred$se, col=2, lty=2)
v.bs <- gm.bs$pred$Salin[which(gm.bs$pred$fit==min(gm.bs$pred$fit))]
lines(c(v.bs, v.bs), c(0, 65), lty = 2, col = 2, lwd = 2)
text( v.bs, 65, as.character(v.bs), pos = 4, cex = 1.3)
mtext('A ', side = 3, adj = 1, line = -1.5, cex = 1.5)

# Chesapeake Bay panel
plot(nsp~Salin, data = samp.chesa, ylab = 'Species richness', xlab='Salinity', ylim = c(0, 70), pch = 19, col = rgb(0,0,0, 0.2), las = 1)
lines(gm.chesa$pred$Salin, gm.chesa$pred$fit, col = 2, lwd= 3)
lines(gm.chesa$pred$Salin, gm.chesa$pred$fit + gm.chesa$pred$se, col=2, lty=2)
lines(gm.chesa$pred$Salin, gm.chesa$pred$fit - gm.chesa$pred$se, col=2, lty=2)
v.chesa <- gm.chesa$pred$Salin[which(gm.chesa$pred$fit == min(gm.chesa$pred$fit))]
lines(c(v.chesa, v.chesa), c(0, 65), lty = 2, col = 2, lwd = 2)
text( v.chesa, 65, as.character(v.chesa), pos = 4, cex = 1.3)
mtext('B ', side = 3, adj = 1, line = -1.5, cex = 1.5)
dev.copy2pdf(file='./Fig2.pdf')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIG 3 RAREFACTION       ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# salinity windows
s.bs <- seq(1, 26, by = 0.5) # salinity windows for the Baltic Sea
s.chesa <- seq(1, 27, by = 0.5) # salinity windows for the Chesapeake Bay

# frequency distribution of samples within salinity windows
x.bs <- x.chesa <- numeric()
for(i in 1:length(s.bs)){ x.bs[i] <- length(which(samp.bs$Salin > s.bs[i]-1 & samp.bs$Salin < s.bs[i]+1))}
for(i in 1:length(s.chesa)){ x.chesa[i] <- length(which(samp.chesa$Salin > s.chesa[i]-1 & samp.chesa$Salin < s.chesa[i]+1))}

# minimum salinity window, fixing the sample based rarefaction
(minx <- min(c(x.bs, x.chesa))) # 108

## iNEXT rarefaction objects ###
lst.bs <- lst.chesa <- list()

#  Baltic Sea
for (i in 1:length(s.bs)) {lst.bs[[i]] <- t(xtab[samp.bs[which(samp.bs$Salin > s.bs[i]-1 & samp.bs$Salin < s.bs[i]+1), 'sampleID'], ])}
names(lst.bs) <- paste('Sal', s.bs, sep = '')

# Chesapeake Bay
for (i in 1:length(s.chesa)) {lst.chesa[[i]] <- t(xtab[samp.chesa[which(samp.chesa$Salin > s.chesa[i]-1 & samp.chesa$Salin < s.chesa[i]+1), 'sampleID'], ])}
names(lst.chesa) <- paste('Sal', s.chesa, sep = '')

# Calculate rarefaction with iNEXT for sample sizes from 1 to 500
print('Rarefaction calculations for the Baltic Sea; takes time')
suppressWarnings(iNEXTse.bs <- iNEXT(lst.bs,  datatype = 'incidence_raw', size = 1:500, se = TRUE))

print('Rarefaction calculations for the Chesapeake Bay; takes time')
suppressWarnings(iNEXTse.chesa <- iNEXT(lst.chesa,  datatype = 'incidence_raw', size = 1:500, se = TRUE))

## iNEXT wrapper functions ##

# coverage based rarefaction curves from iNEXT rarefaction objects
slopeCoverage <- function(inext.lst = iNEXTse.bs, delta = seq(0.3, 2, by = 0.1) ) {
  SC <- abundance <- ucl <- lcl <- matrix(NA, length(inext.lst$iNextEst), length(delta), dimnames = list(names(inext.lst$iNextEst), delta))
  for(i in 1:length(delta)){
    tmp <- unlist(sapply(inext.lst$iNextEst, function(x){which(abs(diff(x$qD) - delta[i]) == min (abs(diff(x$qD) - delta[i])))}))
    a <- b <- u <- l <- numeric()
    for (j in 1:length(inext.lst$iNextEst)){
      a[j] <- inext.lst$iNextEst[[j]]$qD[tmp[j]]     # abundance at the specified slope coverage
      u[j] <- inext.lst$iNextEst[[j]]$qD.UCL[tmp[j]] # upper CL the specified slope coverage
      l[j] <- inext.lst$iNextEst[[j]]$qD.LCL[tmp[j]] # lower CL at the specified slope coverage
      b[j] <- inext.lst$iNextEst[[j]]$SC[tmp[j]]     # iNEXT coverage at the specified slope coverage
    }
    abundance[,i] <- a
    ucl[,i] <- u
    lcl[,i] <- l
    SC[,i] <- b
  }
  return(list(abundance = abundance, ucl = ucl, lcl = lcl, SC = SC))
} # returns richness at predefined (delta) slopeCoverages

# sample based rarefaction curves from iNEXT rarefaction objects
sampleRaref <- function(inext.lst = iNEXTse.bs){
  SC <- abundance <- ucl <- lcl<- matrix(NA, length(inext.lst$iNextEst), length(inext.lst$iNextEst[[1]]$t), dimnames = list(names(inext.lst$iNextEst), inext.lst$iNextEst[[1]]$t))
  for(i in 1:length(inext.lst$iNextEst)){
    abundance[i,] <- inext.lst$iNextEst[[i]]$qD
    ucl[i,] <- inext.lst$iNextEst[[i]]$qD.UCL
    lcl[i,] <- inext.lst$iNextEst[[i]]$qD.LCL
    SC[i,] <- inext.lst$iNextEst[[i]]$SC
  }
  return(list(abundance = abundance, ucl = ucl, lcl = lcl, SC = SC))
} # returns richness at pre defined sample size

sRar.bs <- sampleRaref(iNEXTse.bs)
sRar.chesa <- sampleRaref(iNEXTse.chesa)

# surface salinity frequency distributions for figure background
bs.salinityFrequency <- c(470, 2652, 4009, 5259, 6532, 6502, 7866, 10784, 15575, 20299, 25182, 28376, 25847, 19729, 12607, 5530, 2057, 1459, 1017, 759, 635, 621, 613, 570, 525, 477, 437, 430, 446, 480, 556, 640, 708, 769, 782, 743, 710, 668, 626, 632, 634, 596, 552, 491, 435, 439, 462, 495, 530, 546, 571)
chesa.salinityScale <- c(0.25, 1.50, 3.75, 6.25, 8.75, 11.25, 13.75, 16.50, 19.50, 22.50, 25.50, 28.50)
chesa.salinityFrequency <- c(855, 1097, 1392, 2668, 3483, 8754, 13437, 14297, 8029, 4737, 1232, 16)

quartz('Fig 3', width = 8, height = 4)
par(mfrow=c(1,2), oma = c(1,1,1,1), mar = c(4,4,1,1), yaxs = 'i')

## Baltic Sea panel #
bs.salinityFrequency.scaled <- rescale(c(0, bs.salinityFrequency, 0), c(175, 405))
plot(s.bs, s.bs, type = 'n', las = 1,  xlab = 'Salinity', ylab = 'Rarefied species richness', ylim = c(175, 450))
polygon(c(s.bs[1], s.bs, s.bs[length(s.bs)]), bs.salinityFrequency.scaled, col = gray(0.95) ) # salnity polygon

# sample based rarefaction
lines(s.bs, sRar.bs$abundance[, minx], lwd = 2)
plotCI(s.bs, sRar.bs$abundance[, minx] , ui = sRar.bs$ucl[, minx], li = sRar.bs$lcl[, minx], pch = 19, add = T)

# coverage based rarefaction
lines(s.bs, slopeCoverage(iNEXTse.bs, delta = 0.6)$abundance, lwd=2) 
plotCI(s.bs, slopeCoverage(iNEXTse.bs, delta = 0.6)$abundance, ui = slopeCoverage(iNEXTse.bs, delta = 0.6)$ucl, li = slopeCoverage(iNEXTse.bs, delta = 0.6)$lcl, add = T, pch = 21, pt.bg = 'white')

# complete the plot
lines(c(v.bs, v.bs),c(0, 405), lty = 3, col = 1, lwd = 1) # minimum diversity vertical line
text( v.bs, 405, as.character(v.bs), pos = 4, cex = 1.2) # text the minimum diversity salinity
mtext('A ', side = 3, adj = 1, line = -1.5, cex = 1.5)
legend(x=12, y = 380, legend=c('sample based','coverage based'), pch = c(19, 21), col  = 1, bty ='n')

## Chesapeake Bay panel #
chesa.salinityFrequency.scaled <- rescale(c(0, chesa.salinityFrequency, 0), c(175, 270))
plot(s.chesa, s.chesa, type = 'n', las = 1,  xlab = 'Salinity', ylab = 'Rarefied species richness', ylim = c(175, 300))
polygon(c(chesa.salinityScale[1], chesa.salinityScale, chesa.salinityScale[length(chesa.salinityScale)]), chesa.salinityFrequency.scaled, col = gray(0.95) )

# sample based rarefaction
lines(s.chesa, sRar.chesa$abundance[, minx], lwd = 2)
plotCI(s.chesa, sRar.chesa$abundance[, minx], ui = sRar.chesa$ucl[, minx], li = sRar.chesa$lcl[, minx], pch = 19, add = T)

# coverage based rarefaction
lines(s.chesa, slopeCoverage(iNEXTse.chesa, delta = 0.6)$abundance, lwd = 2) 
plotCI(s.chesa, slopeCoverage(iNEXTse.chesa, delta = 0.6)$abundance, ui = slopeCoverage(iNEXTse.chesa, delta = 0.6)$ucl, li = slopeCoverage(iNEXTse.chesa, delta = 0.6)$lcl, add = T, pch = 21, pt.bg = 'white')

# complete the plot
lines(c(v.chesa, v.chesa),c(0, 300*0.9), lty = 3, col = 1, lwd = 1) # minimum diversity vertical line
text(v.chesa, 300*0.9,as.character(v.chesa), pos = 4, cex = 1.3) # plot the minimum diversity salinity
mtext(' B', side = 3, adj = 0, line = -1.5, cex = 1.5)

dev.copy2pdf(file='./Fig3.pdf')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIG 4 NMDS plots         ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Due to seasonality in phytoplankton community composition, NMDS was restricted to summer season (May to September).

# Salinity regions with high sampling frequency are downsampled

# Baltic Sea
samp0.bs <- subset(samp.bs, month  %in% 5:9 & samp.bs$nsp > 2)
dens.bs <- density(samp0.bs$Salin, from = min(samp0.bs$Salin), to = max(samp0.bs$Salin))
salW.bs <- approx(dens.bs$x, dens.bs$y, xout = samp0.bs$Salin)$y 
idx.bs <- sample(1:nrow(samp0.bs), 3000, replace = F, prob = 1/salW.bs)
sampR.bs <- samp0.bs[idx.bs, ]

# Chesapeake Bay
samp0.chesa <- subset(samp.chesa, month  %in% 5:9)
dens.chesa <- density(samp0.chesa$Salin, from = min(samp0.chesa$Salin), to = max(samp0.chesa$Salin))

salW.chesa <- approx(dens.chesa$x, dens.chesa$y, xout = samp0.chesa$Salin)$y 
idx.chesa <- sample(1:nrow(samp0.chesa), 3000, replace = F, prob = 1/salW.chesa)
sampR.chesa <- samp0.chesa[idx.chesa,]

xtabs.bs <- xtab[rownames(xtab) %in% sampR.bs$sampleID,]; xtabs.bs <- xtabs.bs[, which(colSums(xtabs.bs) > 0)] 
xtabs.chesa <- xtab[rownames(xtab) %in% sampR.chesa$sampleID,]; xtabs.chesa <- xtabs.chesa[, which(colSums(xtabs.chesa) > 0)]

mds.chesa <- metaMDS(xtabs.chesa)  # [takes time]
mds.bs <- metaMDS(xtabs.bs)

scores.bs <- scores(mds.bs)
scores.chesa <- scores(mds.chesa)

quartz('Fig 4', width=8, height = 4); par(mfrow = c(1,2), mar = c(4,4,1,1))

# Baltic Sea
idx <- match(rownames(scores.bs), samp.bs[,'sampleID'])
ord <- order(samp0.bs[idx, 'nsp'])
plot(samp.bs[idx, 'Salin'][ord], scores.bs[ord, 1], las = 1, pch = 21, col = 1, bg = matlab.like(40)[cut(samp.bs[idx, 'nsp'][ord], 40)], ylim = quantile(scores.bs[,1], c(0.001, .999)), xlab = 'Salinity', ylab = 'NMDS1 scores')

txt <- paste(" n = ", nrow(scores.bs), ' ', sep = '')
mtext(' A ', side = 3, adj = 0, line = -1.5, cex = 1.5)
# mtext(txt, side = 3, adj = 1, line = -3, cex = 1.3)
mtext(txt, side = 1, adj = 1, line = -1.5, cex = 1.3)

# add Baltic Sea GAM line
mdsBS.x <- samp.bs[idx, 'Salin']
mdsBS.y <- scores.bs[,1]
mdsBS.gm <- gam(mdsBS.y ~ s(mdsBS.x))
mdsBS.xpred <- seq(0, 31, by = 0.1) 
mdsBS.ypred <- predict(mdsBS.gm, newdata = data.frame(mdsBS.x = mdsBS.xpred))
lines(mdsBS.xpred, mdsBS.ypred, col = 2, lwd = 3)

abline(v = v.bs, lty  = 2, col = 1, lwd = 1) #  minimum line
text( v.bs, .83, as.character(v.bs), pos = 2, cex = 1.3) 

color.legend(15, -0.75, 30, -.55,  as.character(round(seq(min(samp.bs[idx, 'nsp']), max(samp.bs[idx, 'nsp']), length = 5))), 
             matlab.like(40), align = 'lt', gradient = 'x')

# Chesapeake Bay panel
idx <- match(rownames(scores.chesa), samp.chesa[,'sampleID'])
ord <- order(samp0.chesa[idx, 'nsp'])
plot(samp.chesa[idx, 'Salin'][ord], scores.chesa[ord, 1], las = 1, pch = 21, col = 1, bg = matlab.like(40)[cut(samp.chesa[idx, 'nsp'][ord], 40)], ylim = quantile(scores.chesa[,1], c(0.001, .999)), xlab = 'Salinity', ylab = 'NMDS1 scores')

txt <- paste(" n = ", nrow(scores.chesa), ' ', sep = '')
mtext(' B ', side = 3, adj = 0, line = -1.5, cex = 1.5)
mtext(txt, side = 1, adj = 1, line = -1.5, cex = 1.3)

# add Chesapeake Bay GAM line
mdsChesa.x <- samp.chesa[idx, 'Salin']
mdsChesa.y <- scores.chesa[,1]
mdsChesa.gm <- gam(mdsChesa.y ~ s(mdsChesa.x))
mdsChesa.xpred <- seq(0, 31, by = 0.1) 
mdsChesa.ypred <- predict(mdsChesa.gm, newdata = data.frame(mdsChesa.x = mdsChesa.xpred))
lines(mdsChesa.xpred, mdsChesa.ypred, col = 2, lwd = 3)
abline(v=v.chesa, lty  = 2, col = 1, lwd = 1) #  minimum line
text( v.chesa, .77, as.character(v.chesa), pos = 2, cex = 1.3)

color.legend(15, -1, 30, -.79,  as.character(round(seq(3, max(samp.chesa[idx, 'nsp']), length = 5))), 
             matlab.like(40), align = 'lt', gradient = 'x')
dev.copy2pdf(file='./Fig4.pdf')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIG 5 Cumulative Likelihood plot   ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# predicted salinity ranges for Baltic Sea and Chesapeake Bay
predSalin.bs  <-  data.frame(Salin=seq(0, 28, by = 0.1))
predSalin.chesa  <-  data.frame(Salin=seq(0, 33, by =0.1))

xtab.bs <- xtab[match(samp.bs$sampleID, rownames(xtab)),]
xtab.chesa <- xtab[match(samp.chesa$sampleID, rownames(xtab)),]

# turn basin community matrices to binary presence absence tables
xtab.bs.bin <- decostand(xtab.bs, 'pa') # vegan::decostand, presence absence
xtab.chesa.bin <- decostand(xtab.chesa, 'pa')

# leave out rare species
xtab.bs.bin <- xtab.bs.bin[, colSums(xtab.bs.bin) > 3]
xtab.chesa.bin <- xtab.chesa.bin[, colSums(xtab.chesa.bin) > 3]

# fit binary GAM along salinity gradient
res.bs <- matrix(0, nrow = nrow(predSalin.bs), ncol = ncol(xtab.bs.bin))
res.chesa <- matrix(0, nrow = nrow(predSalin.chesa), ncol = ncol(xtab.chesa.bin))

# Species likelihood matrices ##
# Baltic Sea [takes time]
for(i in 1:ncol(xtab.bs.bin)){
  res.bs[,i] <- predict(gam(xtab.bs.bin[, i]~s(Salin), method = 'REML', data = samp.bs, family = binomial), predSalin.bs, se=F, type = 'response')
  print(i)
}
colnames(res.bs) <- colnames(xtab.bs.bin)

# Chesapeake Bay
for(i in 1:ncol(xtab.chesa.bin)){
  pred <- try(predict(gam(xtab.chesa.bin[,i]~s(Salin), method = 'REML', data = samp.chesa, family = binomial), predSalin.chesa, se=F, type = 'response'), silent = TRUE)
  if (class(pred) == "try-error")
  {pred <-  predict(gam(xtab.chesa.bin[,i]~s(Salin), data = salin.chesa, family = binomial), predSalin.chesa, se=F, type = 'response')}
  res.chesa[,i] <- pred
  print(i)
}
colnames(res.chesa) <- colnames(xtab.chesa.bin) # transfer taxon names

# species optimal salinities (mean)
salin.opt.bs <- salin.opt.chesa <- array()
for(i in 1:ncol(res.bs)){
  pred <- res.bs[,i]
  salin.opt.bs[i] <- predSalin.bs[which(pred == max(pred)),]
}
names(salin.opt.bs) <- colnames(xtab.bs.bin) # add taxon names
for(i in 1:ncol(res.chesa)){
  pred <- res.chesa[,i]
  salin.opt.chesa[i] <- predSalin.chesa[which(pred==max(pred)),]
}
names(salin.opt.chesa) <- colnames(xtab.chesa.bin) # add taxon names

# re-order results according to increasing salinity optimum

# Baltic Sea
ord.bs <- order(salin.opt.bs)
res.bs.ord <- res.bs[,ord.bs]
salin.opt.bs.ord <- round(salin.opt.bs[ord.bs], 1)

# Chesapeake Bay
ord.chesa <- order(salin.opt.chesa)
res.chesa.ord <- res.chesa[,ord.chesa]
salin.opt.chesa.ord <- round(salin.opt.chesa[ord.chesa], 1)

# matlab.like color lookup tables

# Baltic Sea
col.bs <-  matlab.like(max(salin.opt.bs)*10+1) 
col.bs <- col.bs[match(round(unique(salin.opt.bs.ord), 1), round(seq(0, max(predSalin.bs), by = 0.1), 1))] 

# Chesapeake Bay
col.chesa <-  matlab.like(max(salin.opt.chesa)*10+1) 
col.chesa <- col.chesa[match(round(unique(salin.opt.chesa.ord), 1), round(seq(0, max(predSalin.chesa), by = 0.1), 1))]

# cumulative likelihood
res.cum.bs.ord <- apply(res.bs.ord, 1, cumsum) # Baltic Sea
res.cum.chesa.ord <- apply(res.chesa.ord, 1, cumsum) # Chesapeake Bay

# unique salinity differences

# Baltic Sea
dif.salin.opt.bs <- which(diff(salin.opt.bs.ord) > 0) 
dif.salin.opt.bs <- c(1, dif.salin.opt.bs, length(salin.opt.bs.ord))

# Chesapeake Bay
dif.salin.opt.chesa <- which(diff(salin.opt.chesa.ord) > 0) 
dif.salin.opt.chesa <- c(1, dif.salin.opt.chesa, length(salin.opt.chesa.ord))

## plot Fig 5 cumulative likelihood #

quartz('Fig 5', width = 8, height = 5)
par(mfrow=c(1,2), oma = c(1,1,1,1), mar = c(4,4,1,1))

# Baltic Sea panel
res.cum.bs.ord[1,] <- 0 # to start the first polygon with zero
plot((predSalin.bs$Salin), res.cum.bs.ord[nrow(res.cum.bs.ord),], type = 'l', xaxs = 'i', yaxs = 'i', ylim = c(0,54),las=1, ylab = 'Cumulative likelihood of presence', xlab = 'Salinity')
xx <- c(predSalin.bs$Salin, rev(predSalin.bs$Salin))
for(i in 1:(length(dif.salin.opt.bs)-1)){
  yy <-   c(res.cum.bs.ord[dif.salin.opt.bs[i],], rev(res.cum.bs.ord[dif.salin.opt.bs[i+1],]))
  polygon(xx, yy, col=col.bs[i], border = col.bs[i], density = NA)
}

# finish Baltic Sea panel
lines(c(v.bs, v.bs ),c(0,35), lty = 2, col = 1, lwd = 2) # minimum diversity vertical line
text( v.bs, 35, as.character(v.bs), pos = 4, cex = 1.3) # plot the minimum diversity salinity
tmp <- abs(round(salin.opt.bs.ord,1) - v.bs)
i  <- which(tmp==min(tmp))[1]
# cumulative boundary - below are taxa with salnity optima below the diversity minimum, and vice versa 
lines(predSalin.bs$Salin, res.cum.bs.ord[i,], col = 1, lwd = 2) 
mtext(' A', side = 3, adj = 0, line = -1.5, cex = 1.5)
color.legend(6, 45, 24, 50,  as.character(round(seq(min(salin.opt.bs), max(salin.opt.bs), length = 5))), 
             matlab.like(40), align = 'lt', gradient = 'x')

# Chesapeake Bay panel
res.cum.chesa.ord[1,] <- 0 # to start the first polygon with zero
plot(predSalin.chesa$Salin, res.cum.chesa.ord[nrow(res.cum.chesa.ord),], type = 'l', ylim = c(0,54),las=1, xaxs = 'i', yaxs = 'i', ylab = 'Cumulative likelihood of presence', xlab = 'Salinity')
xx <- c(predSalin.chesa$Salin, rev(predSalin.chesa$Salin))
for(i in 1:(length(dif.salin.opt.chesa)-1)){
  yy <-   c(res.cum.chesa.ord[dif.salin.opt.chesa[i],], rev(res.cum.chesa.ord[dif.salin.opt.chesa[i+1],]))
  polygon(xx, yy, col=col.chesa[i], border = col.chesa[i], density = NA)
}

lines(c(v.chesa, v.chesa),c(0,35), lty = 2, col = 1, lwd = 2)
text( v.chesa, 35, as.character(v.chesa), pos = 4, cex = 1.3) 
tmp <- abs(round(salin.opt.chesa.ord, 1) - v.chesa)
i  <- which(tmp==min(tmp))[1]
lines(predSalin.chesa$Salin, res.cum.chesa.ord[i,], col = 2, lwd = 2)
mtext(' B', side = 3, adj = 0, line = -1.5, cex = 1.5)
color.legend(7, 45, 27, 50,  as.character(round(seq(min(salin.opt.chesa.ord), max(salin.opt.chesa.ord), length = 5))), 
             matlab.like(40), align = 'lt', gradient = 'x')

dev.copy2pdf(file='./Fig5.pdf')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TABLE 1;  GAM models                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Baltic Sea  ###
modI.bs  <- gam(nsp ~ s(Salin), method='REML', data = samp.bs) 
modII.bs <- gam(nsp ~ s(Salin) + s(jul) + s(Temp) + s(year), method='REML', data = samp.bs)

# detrend
modD.bs <- gam(nsp ~  s(jul) + s(Temp) + s(year), method='REML', data = samp.bs)
tmp.bs <- modD.bs$model # 
tmp.bs$resid <- modD.bs$residuals
tmp.bs$fitted <- modD.bs$fitted.values
tmp.bs$Salin <- samp.bs[rownames(tmp.bs),'Salin']
modIII.bs <- gam(resid ~ s(Salin), method = 'REML', data = tmp.bs)

# Chesapeake Bay  ###
modI.chesa  <- gam(nsp ~ s(Salin), method='REML', data = samp.chesa) 
modII.chesa <- gam(nsp ~ s(Salin) + s(jul) + s(Temp) + s(year), method='REML', data = samp.chesa) 

# detrend
modD.chesa <- gam(nsp ~  s(jul) + s(Temp) + s(year), method='REML', data = samp.chesa)
tmp.chesa <- modD.chesa$model # 
tmp.chesa$resid <- modD.chesa$residuals
tmp.chesa$fitted <- modD.chesa$fitted.values
tmp.chesa$Salin <- samp.chesa[rownames(tmp.chesa),'Salin']
modIII.chesa <- gam(resid ~ s(Salin), method = 'REML', data = tmp.chesa) # dev.exp 18.3

summary(modI.chesa)
summary(modII.chesa)
summary(modIII.chesa)

# prepare Table 1

gmtbl <- function(x){x <- summary(x)
a <- c(paste('n = ', x$n, ';',sep =''),'edf', round(x$s.table[,'edf'], 1) )
b <- c(paste('dev.expl =', sprintf("%1.0f%%", 100*x$dev.expl)), 'F', round(x$s.table[,'F'],1) )
return(data.frame(a = a, b = b))
}

tb1 <- rbind(cbind(gmtbl(modI.chesa), gmtbl(modI.bs)),
             cbind(gmtbl(modII.chesa), gmtbl(modII.bs)),
             cbind(gmtbl(modIII.chesa), gmtbl(modIII.bs)))
tb1$d <- c('Model I','Variable','Salinity','Model II', 'Variable','Salinity','Day of Year','Temperature','Year','Model III','Variable','Salinity')

tbl1 <- as.ltable(tb1[, c(5,1:4)], colbreaks=c(1,0,1,0), rowbreaks=c(1,0,2,1,0,0,0,0,2,1,0), grid = T, colgroups=c('','Chesapeake Bay','Chesapeake Bay','Baltic Sea','Baltic Sea'),
                  justify=c('left','center','center','center','center'),
                  caption = "\\textbf{Generalized additive model (GAM) outputs of the regionally explicit data sets.}" , source.label = "Note: ", source = "Model I -- species richness as a function of salinity; Model II -- species richness as a function of salinity, seasonal trend (day of year), temperature, and long term trend (year); Model III -- species richness detrended with respect to temperature, day of year and year, as a function of salinity. Due to the high sample size (n) all effects are highly significant (p < 0.001). The explanatory power is shown by the deviance explained by each model.", footnote.size = "footnotesize")

tbl1[7] <- "\\hline "

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TABLE 2;  CLAMTEST                  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtab.bs <- xtab[match(samp.bs$sampleID, rownames(xtab)),]
xtab.bs <- xtab.bs[, colSums(xtab.bs) > 0]# 7744  963
xtab.chesa <- xtab[match(samp.chesa$sampleID, rownames(xtab)),]
xtab.chesa <- xtab.chesa[, colSums(xtab.chesa) > 0] # 7896  800

# metadata in the Results section

ncol(xtab) # combined taxon pool in both basins: 1366
ncol(xtab.chesa) # taxon pool in the Chesapeake Bay: 806 taxa
ncol(xtab.bs) # taxon pool in the Baltic Sea: 963 taxa
length(intersect(colnames(xtab.chesa), colnames(xtab.bs))) # 397 common taxa
length(setdiff(colnames(xtab.chesa), colnames(xtab.bs))) # 403 taxa unique to the Chesapeake Bay
length(setdiff(colnames(xtab.bs), colnames(xtab.chesa))) # 566 taxa unique to the Baltic Sea

## CLAMTEST ##

# clamtest solutions on presence-absence data for the Chesapeake Bay and the Baltic Sea

if(!exists('v.chesa')){v.chesa = 7.8}
if(!exists('v.bs')){v.chesa = 8.9}

table(samp.chesa$Salin > v.chesa) # 3295 freshwater and 4601 marine samples in the Chesapeake bay
table(samp.bs$Salin > v.bs) # 6073 freshwater and 1671 marine samples in the Baltic Sea

# run the vegan::clamtest
sol.chesa <- clamtest(decostand(xtab.chesa, method = 'pa'), samp.chesa$Salin > v.chesa)
sol.bs <- clamtest(decostand(xtab.bs, method = 'pa'), samp.bs$Salin > v.bs)

# for clarity change clamtest default column names: Total_FALSE and Total_TRUE to Fresh and Marine 
names(sol.chesa)[2:3] <- names(sol.bs)[2:3] <- c('Fresh', 'Marine')
sol.bs$Species <- as.character(sol.bs$Species)
sol.chesa$Species <- as.character(sol.chesa$Species)

# add WORMS AphiaId
mch <- match(sol.bs$Species, count$name)
sol.bs <- cbind(count[mch, 'AphiaID'], sol.bs)
names(sol.bs)[1] <- 'AphiaID'

mch <- match(sol.chesa$Species, count$name)
sol.chesa <- cbind(count[mch, 'AphiaID'], sol.chesa)
names(sol.chesa)[1] <- 'AphiaID'


if(FALSE){ # make taxa table
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Taxon table code chunk               ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # make taxon table from count table using WORMS web services
  library(taxizesoap)
  
  # taxizesoap::worms_hierarchy wrapper function
  wormscl <- function(x){
    # takes AphiaID as input
    cl <- list()
    for(i in 1:length(x)){
      print(i)
      tmp <- try(worms_hierarchy(ids = x[i]), silent = TRUE)
      if(class(tmp) != 'try-error'){
        cl[[i]]  <- tmp
      } else {
        cl[[i]]  <- NA
      }
    }
    return(cl) }
  
  taxa <- unique(count[, c('AphiaID', 'name')]) # 1378 x 2
  cl <- wormscl(taxa$AphiaID) # classification
  
  cl[[which(taxa$AphiaID == 1)]] <- data.frame(wormsid = 1, rank = 'Biota', scientificname = 'Biota', stringsAsFactors = F) # worms has no true entry for 'Biota'; Biota stands for unidentified taxa
  
  # We do not distinct between botanial and zoological nomenclature codes
  cl1 <- lapply(cl, function(x){x$rank <- sub(' \\(.*\\)','', x$rank); return(x)})
  
  tmp <- data.frame(rank = unique(unlist(lapply(cl1, function(x){x$rank}))), stringsAsFactors = F)
  
  for(i in 1:length(cl1)){
    suppressWarnings(tmp <- merge(tmp, cl1[[i]][,2:3], all = T, by = 'rank'))
    print(i)
  }
  
  tmp1 <- as.data.frame(t(tmp[,-1]), stringsAsFactors = F)
  colnames(tmp1) <- tmp[,'rank']
  
  taxa <- cbind(taxa, tmp1)
  taxa$rank <- unlist(lapply(cl, function(x){tail(x[, 2], 1)})) # add taxonomic resolution
  # taxon table ready
  write.table(taxa, file = './dryad/taxa.txt', quote = FALSE, sep = '\t', row.names= FALSE)
  
  rm(tmp, tmp1) # clean garbage
} # make taxa table

taxa <- read.table(file = './taxa.txt', header = T, sep = '\t', stringsAsFactors = F)

# wrapper function
clamtb <- function(sol){
  # match sol with classification and add class code
  mch <- match(sol$AphiaID, taxa$AphiaID)
  sol$Class <- taxa[mch, 'Class']
  sol$Phylum <- taxa[mch, 'Phylum']
  sol$rank <- taxa[mch, 'rank']
  
  # add spp to genus only
  sol$Species <- ifelse(sol$rank == 'Genus', paste(sol$Species, 'spp.'), sol$Species) 
  
  # aggregate green algae from chlorophyte and streptophyte clades
  
  idx <- which(sol$Phylum %in% c('Chlorophyta','Charophyta'))
  sol[idx, 'Class'] <- 'Chlorophyceae'
  
  # aggregate Ohters
  idx <- which(!sol$Class %in% c('Bacillariophyceae','Dinophyceae','Chlorophyceae','Cyanophyceae'))
  sol[idx, 'Class'] <- 'Others'
  # make 5 categories - marine exclusive, fresh exclusive, marine mix, fresh mix, general mix, rare
  sol$Habitat <- 'Rare'
  sol$Habitat <- ifelse(sol$Fresh == 0, 'Marine excl', sol$Habitat) # unique marine taxa; no freshwater occurrence
  sol$Habitat <- ifelse(sol$Marine == 0, 'Fresh excl', sol$Habitat) # unique fresh taxa; no marine occurrence
  sol$Habitat <- ifelse(sol$Marine >  0 & sol$Fresh >  0 & sol$Classes == 'Specialist_FALSE', 'Fresh mix', sol$Habitat) 
  sol$Habitat <- ifelse(sol$Marine >  0 & sol$Fresh >  0 & sol$Classes == 'Specialist_TRUE', 'Marine mix', sol$Habitat) 
  sol$Habitat <- ifelse(sol$Marine >  0 & sol$Fresh >  0 & sol$Classes == 'Generalist', 'General mix', sol$Habitat)
  sol$Habitat <- ifelse(sol$Classes == 'Too_rare', 'Rare', sol$Habitat) # 
  
  # change Total_FALSE Total_TRUE  to Fresh Marine
  
  # top 10 of each category
  top10 <- list()
  
  tmp <- sol
  tmp <- tmp[which(tmp$rank %in% c('Species','Genus','Subspecies','Variety')),]  # only at least genus level taxa
  
  tmp <- tmp[order(tmp$Habitat, tmp$Marine, decreasing = TRUE),]
  top10[['Marine excl']] <- head(tmp[tmp$Habitat == 'Marine excl', c('Species','AphiaID','Fresh','Marine')], 10)
  top10[['Marine mix']] <- head(tmp[tmp$Habitat == 'Marine mix', c('Species','AphiaID','Fresh','Marine')], 10)
  
  tmp <- tmp[order(tmp$Habitat, tmp$Fresh, decreasing = TRUE),]
  top10[['Fresh excl']] <- head(tmp[tmp$Habitat == 'Fresh excl', c('Species','AphiaID','Fresh','Marine')], 10)
  top10[['Fresh mix']] <- head(tmp[tmp$Habitat == 'Fresh mix', c('Species','AphiaID','Fresh','Marine')], 10)
  
  tmp <- tmp[order(tmp$Habitat, tmp$Marine + tmp$Fresh, decreasing = TRUE),]
  top10[['Generalists']] <- head(tmp[tmp$Habitat == 'General mix', c('Species','AphiaID','Fresh','Marine')], 10)
  
  return(list(sumtb = table(sol$Habitat, sol$Class)[c('Fresh excl','Fresh mix','General mix','Marine mix','Marine excl','Rare'),], top10 = top10))
} # frequency of clamtest categories with taxonomic classes + top10 list

if(TRUE){ # clamtest metadata in Results text and Table 2
  # Chesapeake Bay
  table(sol.chesa$Fresh==0) # 205 unique to marine
  table(sol.chesa[sol.chesa$Fresh==0, 'Classes'])  # of these 17 are marine specialists, 188 are too rare to be classified
  table(sol.chesa$Marine==0) # 141 unique to fresh
  table(sol.chesa[sol.chesa$Marine==0, 'Classes']) # of these 14 are fresh specialists, 127 are too rare to be classified
  
  table(sol.chesa[sol.chesa$Marine > 0 & sol.chesa$Fresh > 0 , 'Classes'])
  # mixed: 454 taxa sperad on both, marine and fresh sides
  # of these 90 are generalists, 157 fresh and  107 marine are marine specialists 100 are too rare to be classified.
  
  # Baltic Sea
  table(sol.bs$Fresh==0) # 111 unique to marine
  table(sol.bs[sol.bs$Fresh==0, 'Classes'])  # of these 45 are marine specialists, 66 are too rare to be classified
  table(sol.bs$Marine==0) # 565 unique to fresh
  table(sol.bs[sol.bs$Marine==0, 'Classes']) # of these 168 are fresh specialists, 397 are too rare to be classified
  
  table(sol.bs[sol.bs$Marine > 0 & sol.bs$Fresh > 0 , 'Classes'])
  # mixed: 287 taxa sperad on both, marine and fresh sides
  # of these 84 are generalists, 80 fresh and  90 marine are marine specialists 33 are too rare to be classified.
} # clamtest metadata

bstb <- clamtb(sol.bs)
chesatb <- clamtb(sol.chesa)

tb.chesa <- as.matrix(chesatb$sumtb)
tb.bs <- as.matrix(bstb$sumtb)

# calculate totals:
tb.chesa <- cbind(tb.chesa, rowSums(tb.chesa))
tb.chesa <- as.data.frame.matrix(rbind(tb.chesa, colSums(tb.chesa)))

tb.bs <- cbind(tb.bs, rowSums(tb.bs))
tb.bs <- as.data.frame.matrix(rbind(tb.bs, colSums(tb.bs)))

(tb2 <- (rbind(tb.chesa, tb.bs)))

colnames(tb2) <- c('Diat', 'Chloro', 'Cyano', 'Dino', 'Others', 'Total')
tb2$Habitat <- rep(c('Fresh excl', 'Fresh', 'Generalists', 'Marine', 'Marine excl', 'Rare', 'Total' ), 2)

tbl2 <- as.ltable(tb2[,c(7, 1:6)], colbreaks=c(0,0,0,0,0,1), rowbreaks=c(0,0,0,0,0,1,2,0,0,0,0,0,1), rowgroups = rep(c('Chesapeake Bay', 'Baltic Sea'), each = 7), grid = T, caption = "\\textbf{Frequency table showing the classification of the observed taxa according to habitat preference.}", source.label = "Note: ", footnote.size = "footnotesize", source = "The rows designate generalists, and habitat specialists of either the fresh water or marine side, and rare taxa. Habitat specialists observed exclusively on the marine or fresh side are shown separately (Fresh excl, Marine excl). The columns distinguish major algal lineages --- diatoms, chlorophytes (pooling classes form the Streptophyta and Chlorophyta clades), cyanobacteria, dinoflagellaes, and other taxa." )

# clean up some latex code
for(i in c('Chesapeake Bay', 'Baltic Sea')){tbl2 <- sub(paste('^.*(',i,').*}', sep = ''),'\\1', tbl2) }

# print Table 1 and Table 2 as pdf files [needs a working latex engine]

as.pdf(c(tbl1, tbl2), stem = 'Tables', preamble = makePreamble(geometryPackage = command("usepackage", options = list(left = "10mm", top = "5mm", bottom = "1mm", right = "10mm"), args = "geometry")))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# APPENDIX A                         ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

top10.bs <- rbind.fill(bstb$top10)
top10.chesa <- rbind.fill(chesatb$top10)

top10.bs$group <- rep(names(bstb$top10), each = 10)
top10.chesa$group <- rep(names(chesatb$top10), each = 10)

# prepare Appendix S1 for latex

ltable.bs <- as.ltable(top10.bs[, -5], caption = "Baltic Sea (n = 7744)", rowgroups = as.character(top10.bs$group), rowgrouplabel = 'Habitat', grid = T, rowgrouprule = 1, source = "The numbers in columns \\textit{Fresh} and \\textit{Marine} designate the number of samples where the species was observed in respective habitats. \\textit{AphiaID} stands for the WoRMS id code of the taxon (http://www.marinespecies.org). n – the total number of samples in a water body. Taxa above genus resolution have not been considered.", source.label = "Note: ", footnote.size = "footnotesize") 

ltable.chesa <- as.ltable(top10.chesa[, -5], caption = "Chesapeake Bay (n = 7896)",rowgroups = as.character(top10.chesa$group), rowgrouplabel = 'Habitat', grid = T, rowgrouprule = 1, colrouprule = 1)

# Appendix caption text
appcap <- c("\\textbf{Appendix A: Supplementary Tables} \\\\ ", "Table A1. Most frequently observed species classified as habitat specialists and generalists in the Chesapeake Bay and the Baltic Sea.")

# concatenate table caption and the two top 10 tables for Chesapeake Bay and the Baltic Sea
appendixA <- c(appcap, ltable.chesa, ltable.bs)

# clean up some latex code
for(i in names(bstb$top10)){appendixA <- sub(paste('^.*(',i,').*}', sep = ''),'\\1', appendixA) }

# compile Appendix
as.pdf(appendixA, stem = 'AppendixA', preamble = makePreamble(
  geometryPackage = command("usepackage", options = list(left = "10mm", top = "5mm", bottom = "1mm", right = "10mm"), args = "geometry"),
  geometry = command("geometry", args = list(paste0("papersize=", paste0("{", 200, "mm", ",", 270, "mm}")))),
  morePreamble = c(command("usepackage", args = "caption"), 
                   command("captionsetup", options = "table", args = "labelformat = empty, justification = raggedright, singlelinecheck = false, textfont = {large}, skip = 0pt"))
))

## end of file ##