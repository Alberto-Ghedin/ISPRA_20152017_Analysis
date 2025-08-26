library(dplyr)
library(lubridate)
library(tibble)
library(jsonlite)
library(zoo)
library(tidyr)
library(openxlsx)
library(indicspecies)
library(patchwork)
library(ggplot2)
library(colorBlindness)
library(ComplexHeatmap)
library(circlize)

HOME_ <- "."
IMAGE_FORMAT <- "pdf"
source(file.path(HOME_, "utils.r"))

plot_indval <- function(df, basins, threshold, title) {
    taxa_list <- order_species(df, basins, threshold = 0)
    complete_indval <- df[, c("Taxon", basins)] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
    pivot_longer(cols = all_of(basins), names_to = "Basin", values_to = "IndVal") 
    complete_indval$Taxon <- factor(complete_indval$Taxon, levels = taxa_list)
    complete_indval$Basin <- factor(complete_indval$Basin, levels = basins)

    taxa_list <- order_species(df, basins, threshold = threshold)
    partial_indval <- df[, c("Taxon", basins)] %>% dplyr::filter(Taxon %in% taxa_list) %>% 
    pivot_longer(cols = all_of(basins), names_to = "Basin", values_to = "IndVal") 
    partial_indval$Taxon <- factor(partial_indval$Taxon, levels = taxa_list)
    partial_indval$Basin <- factor(partial_indval$Basin, levels = basins)



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
        title = title, 
        tag_levels = "A"
        ) + 
    plot_layout(guides = "collect") & 
    theme(
        legend.position = "bottom", 
        legend.box.margin = margin(r = 40), 
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold")
        )
    return(p)
}


plot_indval_CM <- function(df, basins, threshold, title, unique_classes) {

    taxa_list_complete <- order_species(df, basins = basins, threshold = 0)
    complete_indval <- df[, c("Taxon", basins)] %>% column_to_rownames("Taxon") %>% as.matrix()
    complete_indval <- complete_indval[taxa_list_complete, basins] 

    taxa_list_partial <- order_species(df, basins, threshold = threshold)
    partial_indval <- df[, c("Taxon", basins)] %>% column_to_rownames("Taxon") %>% as.matrix()
    partial_indval <- partial_indval[taxa_list_partial, basins]

    palette <- colorBlindness::paletteMartin
    colors <- setNames(rep(palette, length.out = length(unique_classes)), unique_classes)
    
    ha_right_partial <- rowAnnotation(
      Class = genus_class$Class[match(rownames(partial_indval), genus_class$Genus)],
      col = list(Class = colors),
      show_annotation_name = FALSE,
      annotation_legend_param = list(Class = list(
        title_gp = gpar(fontsize = 20, fontface = "bold"),
        labels_gp = gpar(fontsize = 18),
        direction = "horizontal"
      ))
    )

    indval_cell_fun <- function(j, i, x, y, width, height, fill, text_color_threshold = 0.5) {
        value = partial_indval[i, j]
        grid.text(sprintf("%.2f", value), x, y, 
        gp = gpar(fontsize = 15, 
        col = ifelse(value > text_color_threshold, "black", "white")))
    }

    # Create custom color mapping to ensure yellow is at value 1
    color_mapping <- circlize::colorRamp2(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      colors = viridis::viridis(5, begin = 0, end = 1)
    )
    
    ht_partial <- Heatmap(
    partial_indval,
    name = title,
    col = color_mapping,
    column_title = title,
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    column_names_gp = gpar(fontsize = 18),
    column_names_rot = 45,
    column_gap = unit(4, "mm"),
    cell_fun = indval_cell_fun,
    heatmap_legend_param = list(
      title = "IndVal",
      at = seq(0, 1, by = 0.2),
      labels = seq(0, 1, by = 0.2),
      legend_height = unit(10, "cm"),
      legend_width = unit(5, "cm"),
      title_gp = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 18), 
      border = "black"
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontface = "italic", fontsize = 18),
    right_annotation = ha_right_partial
    )
}

save_heatmap <- function(ht, filename, format = "pdf", width = 6, height = 5, res = 300) {
  # Extract base name and make sure extension matches format
  ext <- tolower(format)
  if (!grepl(paste0("\\.", ext, "$"), filename)) {
    filename <- paste0(filename, ".", ext)
  }

  # Open appropriate device
  switch(ext,
         pdf = pdf(filename, width = width, height = height),
         png = png(filename, width = width, height = height, units = "in", res = res),
         svg = svg(filename, width = width, height = height),
         tiff = tiff(filename, width = width, height = height, units = "in", res = res),
         stop("Unsupported format: ", ext)
  )

  # Draw the heatmap
  draw(ht)

  # Close device
  dev.off()
}


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
abund_only_genera <- phyto_abund %>% 
dplyr::filter(!(id == "VAD120" & Date == "2017-04-30")) %>% #filtering outlier from LIG
dplyr::filter(Region != "CAL") %>%
dplyr::filter(Det_level %in% c("Species", "Genus")) %>%
group_by(Date, id, Genus) %>% 
summarize(
    Abund = sum(Num_cell_l), 
    basin = first(New_basin),
    Season = first(Season),
    .groups = "drop"
) %>% pivot_wider(names_from = Genus, values_from = Abund, values_fill = 0)




result <- multipatt(
    log2(abund_only_genera %>% select(-c(Date, id, basin)) +1), 
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
write.xlsx(list_df, paste(HOME_, "indval_only_genera_new_basin_log2.xlsx", sep = "/"), rowNames = TRUE)


    
result_per_basin <- sapply(
    abund_only_genera %>% pull(basin) %>% unique(),
    function(x) {
    indval <- multipatt(
    abund_only_genera %>% dplyr::filter(basin == x) %>% 
    dplyr::select(where(~ sum(is.numeric(.)) != 0)),
    abund_only_genera %>% dplyr::filter(basin == x) %>% pull(Season), 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4), 
    control = how(nperm = 999) #duleg = TRUE
    )
    taxa <- indval$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()
    return(indval$str %>% 
        as.data.frame() %>% 
        filter(rownames(.) %in% taxa) %>% 
        mutate(across(everything(), ~ round(., 3)))
    )
    },
    simplify = FALSE
)
names(result_per_basin) <- abund_only_genera %>% pull(basin) %>% unique()
write.xlsx(result_per_basin, paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), rowNames = TRUE)

result_per_basin <- sapply(
    abund_only_genera %>% pull(basin) %>% unique(),
    function(x) {
    indval <- multipatt(
    abund_only_genera %>% dplyr::filter(basin == x) %>% 
    dplyr::select(where(~ sum(is.numeric(.)) != 0)) %>%
    mutate(across(everything(), ~ log2(. + 1))),
    abund_only_genera %>% dplyr::filter(basin == x) %>% pull(Season), 
    func = "IndVal.g", 
    restcomb = c(1,2,3,4), 
    control = how(nperm = 999) #duleg = TRUE
    )
    taxa <- indval$sign %>% filter(p.value <= 0.05) %>% select(index, stat, p.value) %>% arrange(index) %>% rownames()
    return(indval$str %>% 
        as.data.frame() %>% 
        filter(rownames(.) %in% taxa) %>% 
        mutate(across(everything(), ~ round(., 3)))
    )
    },
    simplify = FALSE
)
names(result_per_basin) <- abund_only_genera %>% pull(basin) %>% unique()
write.xlsx(result_per_basin, paste(HOME_, "indval_only_genera_log2_per_basin.xlsx", sep = "/"), rowNames = TRUE)



sheets <- getSheetNames(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"))
IndVal <- sapply(sheets, function(sheet) {
    data <- read.xlsx(paste(HOME_, "indval_only_genera_per_basin.xlsx", sep = "/"), sheet = sheet)
    colnames(data)[1] <- "Taxon"
    return(data)
}, simplify = FALSE)
names(IndVal) <- sheets

all_taxa <- unique(unlist(lapply(IndVal, function(df) df$Taxon))) 
genus_class <- phyto_abund %>% 
    dplyr::filter(Genus %in% all_taxa) %>%
    dplyr::select(Genus, Class) %>% 
    mutate(
        Genus = factor(Genus, levels = all_taxa, ordered = TRUE),
        Class = case_when(
            Class == "Dinoflagellata incertae sedis" ~ "Dinophyceae", 
            Class == "Cryptophyta incertae sedis" ~ "Cryptophyceae", 
            Class == "Bacillariophyceae" ~ Class,
            Class == "Coccolithophyceae" ~ Class,
            Class == "Dinophyceae" ~ Class,
            TRUE ~ "Else"
        )
    ) %>% 
    arrange(Genus) %>%
    distinct()
unique_classes <- unique(genus_class$Class)


plots <- mapply(
    function(df, basin) {
        p <- plot_indval_CM(
            df %>% mutate(across(where(is.numeric), ~ .^2)), 
            c("Winter", "Spring", "Summer", "Autumn"), 
            0.25, 
            title = paste("Characteristic species in ", basin, sep = ""), 
            unique_classes = unique_classes
            )
    }, 
    IndVal, 
    names(IndVal)
)


indval_path <- paste(HOME_, "IndVal", sep = "/")
dir.create(indval_path, showWarnings = FALSE)
mapply(
    function(p, basin) {
        ggsave(
            p, 
            file = file.path(indval_path, paste(paste("indval_per_basin_", basin, sep = ""), IMAGE_FORMAT, sep = ".")),
            width = 18, height = 12, dpi = 300
        )
    },
    plots,
    names(plots)
)
mapply(
    function(p, basin) {
        save_heatmap(
            p, 
            filename = file.path(indval_path, paste(paste("indval_per_basin_", basin, sep = ""), IMAGE_FORMAT, sep = ".")),
            format = IMAGE_FORMAT,
            width = 10, height = 12, res = 300
        )
    },
    plots,
    names(plots)
)




plots <- mapply(
    function(df, basin) {
        p <- plot_indval_CM(
            df %>% mutate(across(where(is.numeric), ~ .^2)), 
            c("Winter", "Spring", "Summer", "Autumn"), 
            0.25, 
            title = paste("Characteristic species in ", basin, sep = ""), 
            unique_classes = unique_classes
            )
    }, 
    all_data, 
    names(all_data)
)

genera_basin <- sapply(
    names(IndVal), 
    function(basin) {
        genera <- order_species(
            IndVal[[basin]] %>% mutate(across(where(is.numeric), ~ .^2)),
            c("Winter", "Spring", "Summer", "Autumn"),
            threshold = 0.25
        )
        data.frame(
            Genus = genera, 
            basin = basin
        )
    }, 
    simplify = FALSE
) %>% bind_rows()

genera_basin %>% 
group_by(Genus) %>% 
mutate(n = n()) %>%
arrange(desc(n), Genus) %>% write.csv(
    file = file.path(HOME_, "characteristic_genera_frequency.csv"), 
    row.names = FALSE
)



indval_path <- paste(HOME_, "IndVal_log2", sep = "/")
dir.create(indval_path, showWarnings = FALSE)
mapply(
    function(p, basin) {
        save_heatmap(
            p, 
            filename = file.path(indval_path, paste(paste("indval_per_basin_", basin, sep = ""), IMAGE_FORMAT, sep = ".")),
            format = IMAGE_FORMAT,
            width = 10, height = 12, res = 300
        )
    },
    plots,
    names(plots)
)




ha_right_complete <- rowAnnotation(
    Class = genus_class$Class[match(rownames(complete_indval), genus_class$Genus)],
    col = list(Class = colors), 
    show_annotation_name = FALSE
)
# Create custom color mapping for complete heatmap too
complete_color_mapping <- circlize::colorRamp2(
  breaks = c(0, 0.25, 0.5, 0.75, 1.0),
  colors = viridis::viridis(5, begin = 0, end = 1)
)

ht_complete <- Heatmap(
  complete_indval,
  name = "Complete IndVal",
  col = complete_color_mapping,
column_title = "Complete IndVal",
column_title_gp = gpar(fontsize = 18, fontface = "bold"),
show_heatmap_legend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  right_annotation = ha_right_complete
)
