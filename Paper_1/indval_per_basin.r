library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(openxlsx)
library(parallel)
library(indicspecies)
library(readxl)
library(patchwork)
library(ggplot2)

HOME_ <- "."

order_species <- function(df, basins, threshold = 0.5) {
    characteristic_species <- df$Taxon[apply(df[, basins], 1, function(x) any(x >= threshold))]
    df_long <- reshape2::melt(df[, c("Taxon", basins)] %>% dplyr::filter(Taxon %in% characteristic_species), id.vars = "Taxon", variable.name = "Basin", value.name = "Value")
    df_long$Basin <- factor(df_long$Basin, levels = basins, ordered = TRUE)
    ordered_ids <- df_long %>%
    group_by(Taxon) %>%
    summarise(Max_Basin = Basin[which.max(Value)], Max_Value = max(Value)) %>%
    arrange(Max_Basin, desc(Max_Value)) %>% pull(Taxon)
    return(ordered_ids)
}

phyto_abund <- read.csv(paste(HOME_, "phyto_abund.csv", sep = "/")) 
id_basin <- phyto_abund %>% dplyr::select(id, Basin) %>% distinct() %>% column_to_rownames("id")
ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")

result <- multipatt(
    sites_taxa %>% select(-c(Date, id)), 
    id_basin[sites_taxa %>% pull(id), "Basin"], 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5,6), #considering also combination of basis 14 = Lig+STyr, 15 = Lig+WMed, 16 = NA+SA
    control = how(nperm = 999) #duleg = TRUE
    )
taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "indval_per_basin.xlsx", sep = "/"), rowNames = TRUE)



result <- multipatt(
    sites_taxa %>% select(-c(Date, id)) %>% mutate(across(everything(), ~ log1p(.))), 
    basins[sites_taxa %>% pull(id), "Basin"], 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5,6, 14, 15, 16), #considering also combination of basis 14 = Lig+STyr, 15 = Lig+WMed, 16 = NA+SA
    control = how(nperm = 999) #duleg = TRUE
    )
taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "ISPRA_20152017_Analysis/Description/indval_per_basin_log_trasf.xlsx", sep = "/"), rowNames = TRUE)


#indval only genera 
abund_only_genera <- phyto_abund %>% mutate(
    Det_level = case_when(
        Class == "nan" ~ Taxon,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
) %>% group_by(Date, id, Det_level) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    basin = first(Basin),
    .groups = "drop"
) %>% pivot_wider(names_from = Det_level, values_from = Abund, values_fill = 0)


result <- multipatt(
    abund_only_genera %>% select(-c(Date, id, basin)), 
    abund_only_genera %>% pull(basin), 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5,6), 
    control = how(nperm = 999) #duleg = TRUE
    )
taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "indval_only_genera.xlsx", sep = "/"), rowNames = TRUE)


abund_only_genera_single_tyr <- phyto_abund %>% mutate(
    Det_level = case_when(
        Class == "nan" ~ Taxon,
        Genus == "" ~ Class,
        TRUE ~ Genus
    )
) %>% group_by(Date, id, Det_level) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    basin = first(Basin),
    .groups = "drop"
)
abund_only_genera_single_tyr <- abund_only_genera_single_tyr %>% mutate(
    basin = case_when(
        basin == "NorthTyr" ~ "Tyr",
        basin == "SouthTyr" ~ "Tyr",
        TRUE ~ basin
    )
) %>% pivot_wider(names_from = Det_level, values_from = Abund, values_fill = 0)


result <- multipatt(
    abund_only_genera_single_tyr %>% select(-c(Date, id, basin)), 
    abund_only_genera_single_tyr %>% pull(basin), 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4,5), 
    control = how(nperm = 999) #duleg = TRUE
    )
taxa <- result$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()

list_df <- list(
    "Indval" = result$str %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% taxa) %>% 
    mutate(across(everything(), ~ round(., 3))), 
    "Indval_A" = result$A %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3))),
    "Indval_B" = result$B %>%
    as.data.frame() %>%
    filter(rownames(.) %in% taxa) %>%
    mutate(across(everything(), ~ round(., 3)))
)
write.xlsx(list_df, paste(HOME_, "indval_only_genera_single_tyr.xlsx", sep = "/"), rowNames = TRUE)


#PLOTS#
IndVal <- read_excel("./indval_per_basin.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
IndVal <- IndVal %>% dplyr::filter(Taxon != "Other phytoplankton")
taxa_list <- order_species(IndVal, ordered_basins, threshold = 0.5)
partial_indval <- IndVal[, c("Taxon", ordered_basins)] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
pivot_longer(cols = all_of(ordered_basins), names_to = "Basin", values_to = "IndVal") 
partial_indval$Taxon <- factor(partial_indval$Taxon, levels = taxa_list)
partial_indval$Basin <- factor(partial_indval$Basin, levels = ordered_basins)
p2 <- partial_indval %>% ggplot(aes(x = Basin, y = Taxon, fill = IndVal)) +
geom_tile() + 
geom_text(aes(label = round(IndVal, 2), colour = ifelse(IndVal > 0.5, "black", "white")), size = 6) +
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22, face = "italic"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    legend.position = "bottom", 
    #axis.text.y = element_blank(), axis.ticks.y = element_blank()
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
scale_colour_manual(values=c("white"="white", "black"="black")) +
guides(colour = "none") + guides(
    fill = guide_colourbar(
        title = "IndVal", 
        title.position = "left", 
        title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
        label.theme = element_text(size = 20),
        barwidth = unit(20, "lines"),
        ticks.linewidth = 1,
        frame.linewidth = 1,
        ticks.colour = "black",
        frame.colour  ='black'
        )
    )
p2

IndVal <- read_excel("./indval_only_genera.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
IndVal <- IndVal %>% dplyr::filter(Taxon != "Other phytoplankton")
taxa_list <- order_species(IndVal, ordered_basins, threshold = 0)
complete_indval <- IndVal[, c("Taxon", ordered_basins)] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
pivot_longer(cols = all_of(ordered_basins), names_to = "Basin", values_to = "IndVal") 
complete_indval$Taxon <- factor(complete_indval$Taxon, levels = taxa_list)
complete_indval$Basin <- factor(complete_indval$Basin, levels = ordered_basins)

taxa_list <- order_species(IndVal, ordered_basins, threshold = 0.5)
partial_indval <- IndVal[, c("Taxon", ordered_basins)] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
pivot_longer(cols = all_of(ordered_basins), names_to = "Basin", values_to = "IndVal") 
partial_indval$Taxon <- factor(partial_indval$Taxon, levels = taxa_list)
partial_indval$Basin <- factor(partial_indval$Basin, levels = ordered_basins)



p1 <- complete_indval %>% ggplot(aes(x = Basin, y = Taxon, fill = IndVal)) +
geom_tile() + 
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
guides(fill = "none")

p2 <- partial_indval %>% ggplot(aes(x = Basin, y = Taxon, fill = IndVal)) +
geom_tile() + 
geom_text(aes(label = round(IndVal, 2), colour = ifelse(IndVal > 0.5, "black", "white")), size = 6) +
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20, face = "italic"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    #legend.position = "bottom", 
    #axis.text.y = element_blank(), axis.ticks.y = element_blank()
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
scale_colour_manual(values=c("white"="white", "black"="black")) +
guides(colour = "none") + guides(
    fill = guide_colourbar(
        title = "IndVal", 
        title.position = "left", 
        title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
        label.theme = element_text(size = 20),
        barwidth = unit(20, "lines"),
        ticks.linewidth = 1,
        frame.linewidth = 1,
        ticks.colour = "black",
        frame.colour  ='black'
        )
    )


p <- p1 + p2 + 
plot_annotation(
    title = "Indicator value for genera across basins", 
    tag_levels = "A"
    ) + 
plot_layout(guides = "collect") & 
theme(
    legend.position = "bottom", 
    legend.box.margin = margin(r = 40), 
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold")
    )
ggsave(
    p, 
    file = file.path(HOME_, "indval_per_basin.pdf"),
    width = 18, height = 12, dpi = 300
)

IndVal <- read_excel("./indval_only_genera_single_tyr.xlsx", sheet = "Indval", col_names = TRUE) 
colnames(IndVal)[1] <- "Taxon"
IndVal <- IndVal %>% dplyr::filter(Taxon != "Other phytoplankton")
taxa_list <- order_species(IndVal, c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed") , threshold = 0)
complete_indval <- IndVal[, c("Taxon", c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed"))] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
pivot_longer(cols = all_of(c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed")), names_to = "Basin", values_to = "IndVal") 
complete_indval$Taxon <- factor(complete_indval$Taxon, levels = taxa_list)
complete_indval$Basin <- factor(complete_indval$Basin, levels = c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed"))

taxa_list <- order_species(IndVal, c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed"), threshold = 0.5)
partial_indval <- IndVal[, c("Taxon", c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed"))] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
pivot_longer(cols = all_of(c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed")), names_to = "Basin", values_to = "IndVal") 
partial_indval$Taxon <- factor(partial_indval$Taxon, levels = taxa_list)
partial_indval$Basin <- factor(partial_indval$Basin, levels = c("NorthAdr", "SouthAdr", "Ion", "Tyr", "WestMed"))

p1 <- complete_indval %>% ggplot(aes(x = Basin, y = Taxon, fill = IndVal)) +
geom_tile() + 
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
guides(fill = "none")

p2 <- partial_indval %>% ggplot(aes(x = Basin, y = Taxon, fill = IndVal)) +
geom_tile() + 
geom_text(aes(label = round(IndVal, 2), colour = ifelse(IndVal > 0.5, "black", "white")), size = 6) +
theme_bw() +
theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 22, face = "italic"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 25),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    #legend.position = "bottom", 
    #axis.text.y = element_blank(), axis.ticks.y = element_blank()
    ) + 
scale_y_discrete(limits = rev) + 
scale_fill_continuous(type = "viridis", limits = c(0, 1)) + 
scale_colour_manual(values=c("white"="white", "black"="black")) +
guides(colour = "none") + guides(
    fill = guide_colourbar(
        title = "IndVal", 
        title.position = "left", 
        title.theme = element_text(size = 25, face = "bold", margin = margin(r = 30), vjust = 1), 
        label.theme = element_text(size = 20),
        barwidth = unit(20, "lines"),
        ticks.linewidth = 1,
        frame.linewidth = 1,
        ticks.colour = "black",
        frame.colour  ='black'
        )
    )


p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box.margin = margin(r = 40))