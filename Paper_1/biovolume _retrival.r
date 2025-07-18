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


HOME_ <- "./Paper_1/"
source(file.path(HOME_, "utils.r"))
biovolumes <- openxlsx::read.xlsx(
    file.path(HOME_, "bvol_nomp_version_2024.xlsx")
)

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

all_genera <- phyto_abund %>% pull(Genus) %>% unique() %>% sort()
biovolumes_genera <- biovolumes %>% pull(Genus)


#Most of genera have an associated biovolume 
phyto_abund %>% dplyr::filter(Det_level == "Genus") %>%
mutate(
    has_biovolume = ifelse(Genus %in% biovolumes_genera, TRUE, FALSE)
) %>% group_by(Season, Basin, has_biovolume) %>%
summarise(
    Abund = mean(Num_cell_l, na.rm = TRUE),
    .groups = "drop"
) %>% group_by(Season, Basin) %>% 
mutate(Rel_abund = Abund / sum(Abund)) %>% 
ggplot(aes(y = Season, x = Rel_abund, fill = has_biovolume)) +
geom_bar(stat = "identity") +
facet_wrap(~Basin, ncol = 1) +
scale_y_discrete(limits = rev)


biovol_columns <- c(
            "Genus", 
            "Trophy",
            "SizeClassNo", 
            "SizeRange", 
            "Calculated_volume_µm3/counting_unit", 
            "Calculated_Carbon_pg/counting_unit"
        )

phyto_abund <- phyto_abund %>% 
merge(
    biovolumes %>% dplyr::select(all_of(biovol_columns)),
    by = "Genus", all.x = TRUE
)

BV_BM_genera <- biovolumes %>% dplyr::filter(Genus %in% unique(phyto_abund$Genus)) %>% 
dplyr::select(all_of(biovol_columns)) %>% 
group_by(Genus) %>% 
summarise(
    meanBV = mean(`Calculated_volume_µm3/counting_unit`, na.rm = TRUE),
    meanBM_pg = mean(`Calculated_Carbon_pg/counting_unit`, na.rm = TRUE), 
    trophy = first(Trophy)
)

BV_BM_Classes <- biovolumes %>% dplyr::filter(Class %in% unique(phyto_abund$Class)) %>%
dplyr::select(all_of(c("Class", biovol_columns))) %>% 
group_by(Class) %>% 
summarise(
    meanBV = mean(`Calculated_volume_µm3/counting_unit`, na.rm = TRUE),
    meanBM_pg = mean(`Calculated_Carbon_pg/counting_unit`, na.rm = TRUE)
)

BV_BM_Classes %>% write.csv(
    file = file.path(HOME_, "biovolumes_classes.csv"), 
    row.names = FALSE
)
BV_BM_genera %>% write.csv(
    file = file.path(HOME_, "biovolumes_genera.csv"), 
    row.names = FALSE
)
