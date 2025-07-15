one_every_n_item <- function(list, n) {
    return(
        list[seq(1, length(list), n)]
    )
}

plot_variable_along_coast <- function(data, var, group, title, ylab, ordered_latitude = ordered_latitude, ordered_longitude = ordered_longitude) {
    colors <- scales::hue_pal()(length(unique(data[[group]])))
    palette <- setNames(colors, sort(as.character(data[[group]] %>% unique())))
    p <- data %>% 
    mutate(index_id = match(id, params$ordered_id)) %>% 
    ggplot(aes(x = index_id, y = !!as.symbol(var))) + 
    geom_boxplot(aes(group = id, fill = !!as.symbol(group))) +
    scale_fill_manual(values = palette) + 
    labs(title = title, y = ylab) +
    scale_x_continuous(
        name = "Latitude",
        breaks =  seq(1, length(params$ordered_id), 3), 
        labels = one_every_n_item(ordered_latitude, 3),
        limits = c(-0.01, 162.01), 
        sec.axis = sec_axis(
          ~., 
          name = "Longitude",
          breaks = seq(1, length(params$ordered_id), 3),
          labels = one_every_n_item(ordered_longitude, 3)
        )
    ) + 
    ggplot2::theme(
            axis.text.x.bottom = element_text(angle = 45, hjust = 1, size = 17),
            axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 17),
            axis.text.y = element_text(angle = 0, hjust = 0, size = 17),
            axis.title.y = element_text(size = 22),
            axis.title.x = element_text(size = 22),
            plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18, face = "bold"),
        ) 
    return(p)
} 

ggplot_fill_scale <- function(fill_limits, fill_title) {
    z <- list(
        ggplot2::scale_fill_continuous(type = "viridis", limits = fill_limits, oob=scales::squish), 
        ggplot2::scale_colour_manual(values=c("white"="white", "black"="black")), 
        ggplot2::guides(colour = "none"), 
        ggplot2::guides(
        fill = guide_colourbar(
            title = fill_title, 
            title.position = "left", 
            title.theme = element_text(size = 25, face = "bold", margin = margin(r = 3, l = -1), vjust = 1), 
            label.theme = element_text(size = 20),
            barwidth = unit(20, "lines"),
            ticks.linewidth = 1,
            frame.linewidth = 1,
            ticks.colour = "black",
            frame.colour  ='black'
            )
        )
    )
    return(z)     
}

sciencific_notation <- function(x) {
    y <- gsub("\\+0", "", scales::scientific(x, digits = 2))
    y <- gsub("0\\.0e0", "0", y)
    return(y)
}

from_region_to_abreviation <- c(
    "Friuli-Venezia-Giulia" = "FVG",
    "Veneto" = "VEN", 
    "Emilia-Romagna" = "EMR",
    "Marche" = "MAR",
    "Abruzzo" = "ABR",
    "Molise" = "MOL",
    "Puglia" = "PUG",
    "Basilicata" = "BAS",
    "Calabria" = "CAL",
    "Sicilia" = "SIC",
    "Campania" = "CAM", 
    "Lazio" = "LAZ",
    "Toscana" = "TOS",
    "Liguria" = "LIG",
    "Sardegna" = "SAR"
)
ordered_basins <- c("NorthAdr", "SouthAdr", "Ion", "SouthTyr", "NorthTyr", "WestMed")

ordered_transect <- c(
    "SMTS", "SMLG", "VENEZIA", "ROSOLINA", "PORTO_GARIBALDI", "CESENATICO", "RIMINI", "Chienti", "Esino", "GU",
    "VA", "R14001_B2", "FOCE_CAPOIALE", "FOCE_OFANTO", "BARI_TRULLO", "BRINDISI_CAPOBIANCO", "PORTO_CESAREO", "PUNTA_RONDINELLA", "SINNI", "Villapiana",
    "Capo_Rizzuto", "Caulonia_marina", "Saline_Joniche", "Isole_Ciclopi", "Plemmirio", "Isola_Correnti", "San_Marco", "Isole_Egadi", "Capo_Gallo","Vibo_marina", 
    "Cetraro", "Cilento", "Salerno", "Napoli", "Domizio", "m1lt01", "m1lt02",  "m1rm03",  "m1vt04", "Collelungo","Carbonifera",
    "Donoratico", "Fiume_Morto", "Mesco" , "Portofino"  ,  "Voltri"  , "Quiliano" , 
    "Olbia", "Arbatax", "Villasimius", "Cagliari", "Oristano", "Alghero", "Porto_Torres"

)

log_trans <- function(x) {
    eps <- x[x != 0] %>% na.omit() %>% min()
    return(as.numeric(log10(x + eps)))
}

boxcox_transform <- function(values) {
  
  if (any(is.na(values))) {
    clean_values <- as.numeric(na.omit(values))
  } else {
     clean_values <- as.numeric(values)
  }

  if (any(is.infinite(clean_values))) {
    clean_values <- as.numeric(clean_values[is.finite(clean_values)])
  }


  if (min(clean_values) == 0.0) {
    min_val <- clean_values[clean_values != 0] %>% min()
    clean_values <- clean_values + min_val * 0.1
  }
  lambda <- MASS::boxcox(clean_values ~ 1, plotit = FALSE)
  max_lambda <- lambda$x[which(lambda$y == max(lambda$y))]
  if (max_lambda == 0) {
    clean_values <- log(clean_values)
  } else {
    clean_values <- (clean_values^max_lambda - 1) / max_lambda
  }
  values[which(!(is.na(values) | is.infinite(values)))] <- clean_values
  return(values)
}


regression_plot_region <- function(data, var, log_env = FALSE) {
    var_sym <- sym(var)
    data %>% dplyr::select(Region, New_basin, Season, sample_abund, !!var_sym) %>% na.omit() %>% 
    ggplot(aes(y = log10(sample_abund), x = !!var_sym, col = Region)) + 
    geom_point() + 
    facet_wrap(~New_basin, scales = "free") + 
    geom_smooth(aes(x = !!var_sym, y = log10(sample_abund), col = Region), method = "lm") + 
    labs(title = var)
}


compute_TRIX <- function(data) {
    max_chla <- data %>% pull(Chla) %>% max(na.rm = TRUE) 
    min_chla <- data %>% pull(Chla) %>% min(na.rm = TRUE)
    max_din <- data %>% mutate(DIN = NO3 +  NH4 ) %>% pull(DIN) %>% max(na.rm = TRUE)
    min_din <- data %>% mutate(DIN = NO3 +  NH4 ) %>% pull(DIN)
    min_din <- min_din[min_din != 0] %>% min(na.rm = TRUE)
    max_tp <- data %>% pull(TP[TP !=0]) %>% max(na.rm = TRUE) 
    min_tp <- data %>% pull(TP[TP !=0]) 
    min_tp <- min_tp[min_tp != 0] %>% min(na.rm = TRUE) 
    max_o_dev <- abs(data %>% pull(O_sat) - 100) %>% max(na.rm = TRUE)
    min_o_dev <- abs(data %>% pull(O_sat) - 100)
    min_o_dev <- min_o_dev[min_o_dev != 0] %>% min(na.rm = TRUE)
 output <- data %>% 
 mutate(
    Chla = ifelse(Chla == 0, min_chla, Chla),
    DIN = ifelse(NO3 + NH4 == 0, min_din, NO3 + NH4),
    TP = ifelse(TP == 0, min_tp, TP),
    O_dev = ifelse(abs(O_sat - 100) == 0, min_o_dev, abs(O_sat - 100))
    ) %>% 
mutate(
    TRIX = 10 / 4 * (
        (log10(Chla / min_chla)) / (log10(max_chla / min_chla)) +
        (log10(DIN / min_din)) / (log10(max_din / min_din)) +
        (log10(TP / min_tp)) / (log10(max_tp / min_tp)) + 
        (log10(O_dev / min_o_dev)) / (log10(max_o_dev / min_o_dev))
    )
    )
    return(output$TRIX)
}

compute_TRIX_simplified <- function(data) {
    min_chla <- data %>% pull(Chla) %>% min(na.rm = TRUE)
    min_din <- data %>% mutate(DIN = NO3 +  NH4 ) %>% pull(DIN)
    min_din <- min_din[min_din != 0] %>% min(na.rm = TRUE)
    min_tp <- data %>% pull(TP[TP !=0]) 
    min_tp <- min_tp[min_tp != 0] %>% min(na.rm = TRUE) 
    min_o_dev <- abs(data %>% pull(O_sat) - 100)
    min_o_dev <- min_o_dev[min_o_dev != 0] %>% min(na.rm = TRUE)
 output <- data %>% 
 mutate(
    Chla = ifelse(Chla == 0, min_chla * 893.5, Chla * 893.5),
    DIN = ifelse(NO3 + NH4 == 0, min_din, NO3 * 62 + NH4 * 18),
    TP = ifelse(TP == 0, min_tp * 31, TP * 31),
    O_dev = ifelse(abs(O_sat - 100) == 0, min_o_dev, abs(O_sat - 100))
    ) %>% 
mutate(
    TRIX = 12 / 10 * (log10(Chla * DIN * TP * O_dev) + 1.5)
    )
    return(output$TRIX)
}


fit_reg_basin <- function(group, group_col, vars) {
    if (nrow(cleaned_data) < 3) {
        return(c(NA, NA))
    }
        model <- lm(
            paste("log10(sample_abund) ~", paste(vars, collapse = " + ")), 
            data = cleaned_data %>% dplyr::filter(!!sym(group_col) == group)
        )
        return(
            summary(model)$coefficients[-1, c(1,4)] %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Var") %>% 
            rename(p_value = `Pr(>|t|)`)
        )
    }

process_abund_groups <- function(phyto_abund) {
    abund_groups <- phyto_abund %>% mutate(
        higher_group = case_when(
            Class == "Dinoflagellata incertae sedis" ~ "Dinoflagellata",
            Taxon == "Noctilucea" ~ "Dinoflagellata",
            Class == "nan" ~ "Unknown", 
            Class == "Dinophyceae" ~ "Dinoflagellata", 
            TRUE ~ Class
        )
    ) %>% group_by(Date, id, higher_group) %>%
    summarise(
        Abund = sum(Num_cell_l),
        Region = first(Region),
        Season = first(Season), 
        Basin = first(Basin), 
        .groups = "drop"
    ) %>% 
    mutate(
        Group = case_when(
            higher_group == "Dinoflagellata" ~ "DIN",
            higher_group == "Bacillariophyceae" ~ "DIA",
            higher_group == "Coccolithophyceae" ~ "COC",
            higher_group == "Cryptophyceae" ~ "CRY",
            higher_group != "Unknown" ~ "Else",
            higher_group == "Unknown" ~ "UNK"
        )
    ) %>% group_by(Date, id, Group) %>%
    summarise(
        Abund = sum(Abund),
        Region = first(Region),
        Season = first(Season),
        New_basin = first(New_basin),
        .groups = "drop"
    )
    abund_groups$Region <- factor(abund_groups$Region, levels = unname(from_region_to_abreviation), ordered = TRUE)
    abund_groups$Season <- factor(abund_groups$Season, levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
    abund_groups$New_basin <- factor(abund_groups$New_basin, levels = c("NA", "CA", "SA", "SM", "SIC", "ST", "NT", "LIG", "SAR"), ordered = TRUE)
    return(abund_groups)
}

