plot_projection_ts_SKJ <- function (object, type = "SSB", MP = NULL, n_samples = 3, probs = c(0.1, 0.5), 
                                    ribbon_colours = RColorBrewer::brewer.pal(8, "Blues")[c(2, 4, 8)], 
                                    bbmsy_zones = c(0.4, 0.8), catch_reference = 1, clip_ylim = NULL, seed = 42) 
{
  if (is.null(object@OM$CurrentYr[[1]])) {
    warning("Missing `object@OM$CurrentYr`.\n", "Please run the MSE with a newer GitHub DLMtool version\n", 
            "or set `object@OM$CurrentYr` yourself.\n", 
            "Setting CurrentYr = 0 for now.", call. = FALSE)
    this_year <- 0
  }
  else {
    this_year <- object@OM$CurrentYr[[1]]
  }
  if(is.null(MP)) MP <- object@MPs
  ts_data <- gfdlm:::get_ts(object = object, type = type, this_year = this_year) %>% filter(mp_name %in% MP)
  ts_data$value[ts_data$Type == "B_BMSY"] <- 0.4 * ts_data$value[ts_data$Type == "B_BMSY"]
  quantiles <- gfdlm:::get_ts_quantiles(ts_data, probs = probs)
  set.seed(seed)
  sampled_ids <- sample(unique(ts_data$iter), size = n_samples)
  d <- dplyr::filter(ts_data, .data$iter %in% sampled_ids)
  .type_labels <- gsub("_", "/", unique(d$Type))
  .type_labels <- gsub("MSY", "[0]", .type_labels)
  type_df <- data.frame(Type = unique(d$Type), type_labels = .type_labels, stringsAsFactors = FALSE)
  mp_names <- sort(unique(d$mp_name))
  ref_grep <- grepl("ref", mp_names)
  if (any(ref_grep)) {
    mp_names <- c(mp_names[!ref_grep], mp_names[ref_grep])
  }
  d <- dplyr::left_join(d, type_df, by = "Type")
  quantiles <- dplyr::left_join(quantiles, type_df, by = "Type")
  d$mp_name <- factor(d$mp_name, levels = mp_names)
  quantiles$mp_name <- factor(quantiles$mp_name, levels = mp_names)
  
  g <- ggplot(d, aes_string("real_year", "value", group = "iter"))
  g <- g + ggplot2::geom_ribbon(data = quantiles, aes_string(x = "real_year", ymin = "ll", ymax = "uu"), colour = NA, fill = ribbon_colours[1], 
                                inherit.aes = FALSE)
  g <- g + ggplot2::geom_ribbon(data = quantiles, aes_string(x = "real_year", ymin = "l", ymax = "u"), colour = NA, fill = ribbon_colours[2], 
                                inherit.aes = FALSE)
  g <- g + ggplot2::geom_line(data = quantiles, aes_string(x = "real_year", y = "m"), colour = ribbon_colours[3], lwd = 1, 
                              inherit.aes = FALSE)
  
  if ("C" %in% type) {
    #average_catch <- filter(d, Type == "Catch", real_year %in% 
    #                          .hist_years[(length(.hist_years) - (catch_reference - 
    #                                                                1)):length(.hist_years)]) %>% summarize(average_catch = mean(value)) %>% 
    #  pull(average_catch)
    #g <- g + ggplot2::geom_hline(yintercept = average_catch, 
    #                             alpha = 0.2, lty = 2, lwd = 0.5)
  } else {
    lines <- data.frame(value = bbmsy_zones, type_labels = "B/B[0]", stringsAsFactors = FALSE)
    lines$type_labels <- factor(lines$type_labels, levels = "B/B[0]")
    g <- g + geom_hline(data = lines, mapping = aes(yintercept = value), alpha = 0.2, lty = 2, lwd = 0.5)
  }
  g <- g + ggplot2::geom_line(alpha = 0.3, lwd = 0.4) + ggplot2::ylab("Value") + 
    ggplot2::xlab("Year") + theme_pbs() + ggplot2::facet_grid(mp_name ~ 
                                                                type_labels, labeller = ggplot2::label_parsed) + ggplot2::coord_cartesian(expand = FALSE, 
                                                                                                                                          ylim = if (is.null(clip_ylim)) 
                                                                                                                                            NULL
                                                                                                                                          else c(0, clip_ylim * max(quantiles$m))) + ggplot2::theme(panel.spacing = grid::unit(-0.1, 
                                                                                                                                                                                                                               "lines")) + geom_vline(xintercept = this_year, 
                                                                                                                                                                                                                                                      lty = 2, alpha = 0.3)
  g
}


plot_projection_SKJ <- function (object, catch_breaks = NULL, catch_labels = NULL, rel_widths = c(2, 1.18), 
                                 msy_ylim = c(0, 1.4), catch_ylim = NULL, MP = NULL) {
  suppressMessages({
    g1 <- plot_projection_ts_SKJ(object, MP = MP, bbmsy_zones = 0.4) + ggplot2::coord_cartesian(expand = FALSE, ylim = msy_ylim) + ggplot2::theme(strip.text.y = ggplot2::element_blank())
  })
  g2 <- plot_projection_ts_SKJ(object, type = "C", MP = MP, catch_reference = 1)
  if (!is.null(catch_ylim)) {
    suppressMessages({
      g2 <- g2 + ggplot2::coord_cartesian(expand = FALSE, ylim = catch_ylim)
    })
  }
  if (!is.null(catch_breaks) && is.null(catch_labels)) {
    catch_labels <- catch_breaks
  }
  if (!is.null(catch_breaks)) {
    suppressMessages({
      g2 <- g2 + ggplot2::scale_y_continuous(breaks = catch_breaks, 
                                             labels = catch_labels)
      g2 <- g2 + ggplot2::coord_cartesian(expand = FALSE, 
                                          ylim = catch_ylim)
    })
  }
  cowplot::plot_grid(g1, g2, rel_widths = rel_widths, align = "h")
}
