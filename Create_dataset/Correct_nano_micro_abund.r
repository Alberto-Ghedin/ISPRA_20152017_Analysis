library(tidyr)
library(tibble)
library(dplyr)
library(lubridate)
HOME_ <- path.expand("~")

stations_info <- read.csv(
    paste(HOME_, "PHD/ISPRA_20152017_Analysis/Stations_info.csv", sep = "/")
) 
abund_raw <- read.csv(
    paste(HOME_, "PHD/ISPRA_20152017_Analysis/phyto_nano_micro.csv", sep = "/")
) %>% rename(
    id = NationalStationID
)

abund_raw <- abund_raw %>% mutate(
    id = case_when(
        id == "0004-M000200_TR06" ~ "0004-MS00200_TR06",
        id == "0004-M000200_TR12" ~ "0004-MS00200_TR12",
        TRUE ~ as.character(id)
    )
)

abund_raw <- abund_raw %>% mutate(
    Date = as.Date(paste(Year, Month, Day, sep = "-"), format = "%Y-%m-%d"),
)

abund_raw <- abund_raw %>% merge(
    stations_info %>% select(id, Region),
    by = "id", all.x = TRUE
)

abund_raw <- abund_raw %>% 
mutate(
    Abbondanza_CD = as.numeric(
        case_when(
        Abbondanza_CD == "" ~ "0",
        Abbondanza_CD == "<120" ~ "120",
        TRUE ~ Abbondanza_CD
    )
    )
)

phyto_Abund <- read.csv(
    paste(HOME_, "PHD/ISPRA_20152017_Analysis/Paper_1/phyto_abund.csv", sep = "/")
)

phyto_Abund %>% dplyr::filter(id == "M1T1B" & Date == "2017-05-31")

abund_raw <- abund_raw %>% mutate(
    SampleDepth = as.numeric(
        case_when(
            SampleDepth == "" ~ NA,
            SampleDepth == "0-25" ~ NA,
            SampleDepth == "0-30" ~ NA,
            SampleDepth == "0-50" ~ NA,
            TRUE ~ SampleDepth
        )
    )
)

abund_raw %>% dplyr::filter(SampleDepth < 2) %>% group_by(Date, id) %>% 
summarise(
    n_samples = n_distinct(SampleDepth)
)

colnames(abund_raw)

abund_raw <- abund_raw %>% mutate(
    Date = case_when(
        Region == "Emilia-Romagna" & (Date == as.Date("2015-10-05") | Date == as.Date("2015-10-09")) ~ as.Date("2015-09-30"),
        Region == "Marche" & (Date == as.Date("2016-03-08") | Date == as.Date("2016-03-02")) ~ as.Date("2016-02-27"),
        Region == "Marche" & Date == as.Date("2017-04-04") ~ as.Date("2017-03-30"),
        Region == "Marche" & Date == as.Date("2017-05-02") ~ as.Date("2017-04-30"),
        Region == "Veneto" & Date == as.Date("2016-05-05") ~ as.Date("2016-04-30"),
        TRUE ~ Date
    )
)

abund_raw %>% dplyr::filter(FitoZoo == "F" & SampleDepth < 1) %>% group_by(id, Date) %>% summarise(
    n = n()
) %>% dplyr::filter(n > 2)

abund_raw %>% dplyr::filter(FitoZoo == "F" & SampleDepth < 1) %>% dplyr::filter(
    id == "1E_MS_CH_6" & Date == "2016-01-20"
)

abund_raw %>% dplyr::filter(FitoZoo == "F") %>% group_by(id, Year, Month, ClasseDim) %>% 
arrange(SampleDepth) %>% 
summarise(
    Date = first(Date), 
    SampleDepth = first(SampleDepth),
    Abbondanza_CD = first(Abbondanza_CD), 
    .groups = "drop"
) %>% dplyr::select(id, Date, ClasseDim, Abbondanza_CD) %>%
mutate(
    Date = lubridate::ceiling_date(Date, unit = "month") - days(1),
) %>% 
write.csv(
    paste(HOME_, "PHD/ISPRA_20152017_Analysis/phyto_nano_micro_corrected.csv", sep = "/"),
    row.names = FALSE
)
