library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(rjson)
library(vegan)

IMAGE_FORNMAT <- "pdf"
HOME_ <- "."
source(file.path(HOME_, "utils.r"))

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% dplyr::select(id, Transect,SeaDepth)
)
phyto_abund$Region <- from_region_to_abreviation[as.character(phyto_abund$Region)]
phyto_abund$Transect <- factor(phyto_abund$Transect, levels = ordered_transect, ordered = TRUE)
phyto_abund <- phyto_abund %>% mutate(
    Basin = case_when(
        Region %in% c("FVG", "VEN", "EMR") ~ "NA",
        Region %in% c("MAR", "ABR") ~ "CA", 
        Region == "MOL" ~ "SA",
        Transect %in% c("FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO") ~ "SA",
        Transect %in% c("PORTO_CESAREO", "PUNTA_RONDINELLA") ~ "SM",
        Region == "BAS" ~ "SM",
        Transect %in% c("Villapiana", "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche") ~ "SM",
        Region == "SIC" ~ "SIC", 
        Transect %in% c("Vibo_marina", "Cetraro") ~ "ST",
        Region == "CAM" ~ "ST",
        Transect %in% c("m1lt01", "m1lt02") ~ "ST",
        Transect %in% c("m1rm03", "m1vt04") ~ "NT",
        Region == "TOS" ~ "NT",
        Region == "LIG" ~ "LIG",
        Region == "SAR" ~ "SAR"
    )
)
phyto_abund$Basin <- factor(phyto_abund$Basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)


abund_only_genera <- phyto_abund %>%
        dplyr::filter(Det_level %in% c("Genus", "Species")) %>%
        group_by(id, Date, Genus) %>% 
        summarise(
            abund = sum(Num_cell_l),
            Region = first(Region),
            Season = first(Season), 
            .groups = "drop"
        ) %>% dplyr::filter(Genus != "") 

abund_only_genera %>% pull(Genus) %>% unique() %>% length()

#Testign for presence of spacetime structure
genera_abund.hel <- abund_only_genera %>%
tidyr::pivot_wider(
    id_cols = c(id, Date), 
    names_from = Genus,
    values_from = abund,
    values_fill = 0
) %>% dplyr::select(-c(id, Date)) %>% as.matrix() %>% log1p() %>% vegan::decostand(method = "hellinger") 

rank_freq <- colSums(genera_abund.hel > 0) %>% sort() %>% rev()
plot(c(1:length(rank_freq)), rank_freq, xlim = c(0, 100))

genera_high <- names(rank_freq[rank_freq >= 160])


space_time.mat <- abund_only_genera %>%
    tidyr::pivot_wider(
        id_cols = c(id, Date), 
        names_from = Genus,
        values_from = abund,
        values_fill = 0
    ) %>% dplyr::select(id, Date) %>% merge(
        phyto_abund %>% 
            dplyr::select(id, Date, Region, Season, Longitude, Latitude) %>% arrange(id, Date) %>% distinct(),
        by = c("id", "Date")
    ) 


genera.abund.hel.res <- resid(lm(
    genera_abund.hel ~ Longitude * Latitude * Season,
    data = space_time.mat
    )
)


MEMs <- read.csv(
    file.path(HOME_, "MEMs_all_sites.csv")
) %>% dplyr::select(-c(Longitude, Latitude))

space_time.mat %>% head()
genera.dbmems <- rda(
    genera.abund.hel.res  ~ .,
    data = merge(
        space_time.mat, 
        MEMs,
        by = "id"
    ) %>% dplyr::select(-c(id, Longitude, Latitude, Date, Region, Season))
)

RsquareAdj(genera.dbmems)
genera.dbmems.fwd <- adespatial::forward.sel(
    Y = genera.abund.hel.res,
    X =  merge(
        space_time.mat, 
        MEMs,
        by = "id"
    ) %>% dplyr::select(-c(id, Longitude, Latitude, Date, Region, Season)) %>% as.matrix(),
    Xscale = FALSE,
    adjR2thresh = RsquareAdj(genera.dbmems)$adj.r.squared,
    nperm = 999
)
genera.dbmems.fwd
sort(genera.dbmems.fwd[, 2])

dbmem.sign <- genera.dbmems.fwd %>% pull(variables)
genera.dbmems.sign <- rda(
    genera.abund.hel.res  ~ .,
    data = merge(
        space_time.mat, 
        MEMs,
        by = "id"
    ) %>% dplyr::select(all_of(dbmem.sign))
)

anova(genera.dbmems.sign, by = "axis", permutations = 999)
