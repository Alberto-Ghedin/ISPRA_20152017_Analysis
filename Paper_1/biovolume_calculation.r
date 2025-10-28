library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
#library(ggh4x)
#library(legendry)
library(rjson)
library(openxlsx)
library(colorBlindness)


HOME_ <- "."
IMAGE_FORMAT <- "svg"
source(file.path(HOME_, "utils.r"))

sea_depth <- read.csv(file.path(HOME_, "transects_info.csv"))
params <- fromJSON(file = file.path(HOME_, "params.json"))

phyto_abund <- read.csv(file.path(HOME_, "phyto_abund.csv")) %>% dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% 
merge(
    sea_depth %>% select(id, Transect,SeaDepth)
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


BV_BM_Classes <- read.csv(file.path(HOME_, "biovolumes_classes.csv")) 
BV_BM_genera <- read.csv(file.path(HOME_, "biovolumes_genera.csv")) 


BV_BM_genera %>% mutate(
    radius = (3 * meanBV / (4 * pi))^(1/3)
) %>% dplyr::select(Genus, radius) %>% arrange(radius) %>% head(30)

genus_abund_biovol <- phyto_abund %>% dplyr::filter(Det_level == "Genus") %>% 
    group_by(Date, id, Genus) %>% 
    summarise(
        abund = sum(Num_cell_l, na.rm = TRUE),
        Season = first(Season),
        Transect = first(Transect),
        Region = first(Region),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% 
    merge(
        BV_BM_genera %>% select(Genus, meanBV, meanBM_pg),
        by = "Genus"
    ) %>% mutate(
        BV = meanBV * abund,
        BM_pg = meanBM_pg * abund
    ) 

class_abund_biovol <- phyto_abund %>% 
mutate(
    Class = case_when(
        Taxon == "Noctilucea" ~ "Dinophyceae",
        Class == "Dinoflagellata incertae sedis" ~ "Dinophyceae", 
        Taxon == "Dinoflagellata" ~ "Dinophyceae",
        TRUE ~ Class
    )
) %>% dplyr::filter(Class != "nan" & Det_level == "Higher cat.") %>% 
    group_by(Date, id, Class) %>% 
    summarise(
        abund = sum(Num_cell_l, na.rm = TRUE),
        Season = first(Season),
        Transect = first(Transect),
        Region = first(Region),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% 
    merge(
        BV_BM_Classes %>% select(Class, meanBV, meanBM_pg),
        by = "Class"
    ) %>% mutate(
        BV = meanBV * abund,
        BM_pg = meanBM_pg * abund
    )


unk_abund_biovol <- phyto_abund %>% dplyr::filter(Taxon == "Other phytoplankton (inf. 20µm)") %>% 
mutate(
    BV = Num_cell_l * 10^3 * (4/3) * pi,  # Assuming a spherical shape with diameter 20µm, 
    BM_pg = exp(-0.6573) * BV ^ 0.933 #extracted from linear model based on genera abundances biomass
)



ordered_latitude <- phyto_abund %>% dplyr::select(id, Latitude) %>% distinct() %>% 
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Latitude) %>% round(2)
ordered_longitude <- phyto_abund %>% dplyr::select(id, Longitude) %>% distinct() %>%
arrange(id = factor(id, levels = params$ordered_id, ordered = TRUE)) %>% pull(Longitude) %>% round(2)

rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    unk_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% group_by(Date, id) %>% 
    summarise(
        BV = log10(sum(BV, na.rm = TRUE)),
        BM_pg = log10(sum(BM_pg, na.rm = TRUE)),
        Region = first(Region),
        Season = first(Season),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% plot_variable_along_coast(
        var = "BM_pg", 
        group = "Region",
        title = "Biovolume (µm³/l)", 
        ylab = "Biovolume (µm³/l)", 
        ordered_latitude = ordered_latitude,
        ordered_longitude = ordered_longitude
    )

genus_abund_biovol %>%
ggplot() + 
geom_point(aes(x = log10(BV), y = log10(BM_pg)), size = 0.5)

lm(
    log10(BM_pg) ~ log10(BV), data = genus_abund_biovol
) %>% summary()


rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    unk_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% group_by(Date, id) %>% 
    summarise(
        BV = sum(BV, na.rm = TRUE),
        BM_pg = sum(BM_pg, na.rm = TRUE),
        Region = first(Region),
        Season = first(Season),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% merge(
    phyto_abund %>% group_by(Date, id) %>% 
    summarise(
        sample_abund = sum(Num_cell_l, na.rm = TRUE),
        .groups = "drop"
    ) %>% dplyr::select(Date, id , sample_abund)
) %>% ggplot() + 
geom_point(aes(y = log10(BM_pg), x = log10(sample_abund), col = Region), size = 2)# + 
stat_ellipse(aes(x = log10(BM_pg), y = log10(sample_abund), col = Region), level = 0.95, linewidth = 0.5)

rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% group_by(Date, id) %>% 
    summarise(
        BV = sum(BV, na.rm = TRUE),
        BM_pg_known = sum(BM_pg, na.rm = TRUE),
        Region = first(Region),
        Season = first(Season),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% merge(
        unk_abund_biovol %>% rename(BM_pg_unk = BM_pg) %>% dplyr::select(Date, id, BM_pg_unk), 
        by = c("Date", "id")
    ) %>% ggplot() + 
geom_histogram(aes(x = log10(BM_pg_unk / BM_pg_known), fill = Region)) + 
facet_wrap(~Season) 

rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    unk_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% group_by(Date, id) %>% 
    summarise(
        BV = sum(BV, na.rm = TRUE),
        BM_pg = sum(BM_pg, na.rm = TRUE),
        Region = first(Region),
        Season = first(Season),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% ggplot() + 
geom_histogram(aes(x = log10(BM_pg), fill = Region)) #+ 
facet_wrap(~Season) 

rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    unk_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% group_by(Date, id) %>% 
    summarise(
        BV = sum(BV, na.rm = TRUE),
        BM_pg = sum(BM_pg, na.rm = TRUE),
        Region = first(Region),
        Season = first(Season),
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% merge(
    phyto_abund %>% group_by(Date, id) %>% 
    summarise(
        sample_abund = sum(Num_cell_l, na.rm = TRUE),
        .groups = "drop"
    ) %>% dplyr::select(Date, id , sample_abund)
) %>% ggplot() + 
geom_boxplot(aes(x = Season, y = log10(BM_pg))) + 
facet_wrap(~Basin, scale = "free_y")


genus_abund_biovol %>% head()
rbind(
    genus_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    class_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin), 
    unk_abund_biovol %>% dplyr::select(Date, id, BV, BM_pg, Region, Season, Basin)
) %>% head()

#NOT INSIGHTFUL         
sheets <- getSheetNames(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"))
IndVal <- sapply(sheets, function(sheet) {
    data <- openxlsx::read.xlsx(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), sheet = sheet, colNames = TRUE)
    colnames(data)[1] <- "Taxon"
    data <- data# %>% dplyr::filter(Taxon != "Other phytoplankton")
    return(data)
}, simplify = FALSE)
names(IndVal) <- sheets



basin <- "NA"
selected_species <- order_species(IndVal[[basin]], c("Winter", "Spring", "Summer", "Autumn"), threshold = 0) 

IndVal[[basin]] %>% dplyr::filter(Taxon %in% selected_species) %>%
mutate(Genus = factor(Taxon, levels = selected_species, ordered = TRUE)) %>%
merge(
    BV_BM_genera %>% dplyr::select(Genus, meanBV, meanBM_pg, trophy) %>% distinct(Genus, .keep_all = TRUE),
    by = "Genus"
) %>% arrange(Genus) %>% 
mutate(dummy_y = 1) %>%
ggplot() + 
geom_tile(aes(y = Genus, x = dummy_y, fill = trophy)) +
scale_y_discrete(limits = rev) + 
theme_minimal() + 
labs(x = "Genus", y = "Mean Biomass (pg)", fill = "Mean Biomass (pg)")



