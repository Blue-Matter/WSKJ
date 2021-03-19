
library(MSEtool)
library(dplyr)
library(ggplot2)
library(purrr)

#devtools::install_github("pbs-assess/gfdlm")
library(gfdlm)
source("03a_MPs.R")


MSE <- lapply(1:6, function(i) readRDS(paste0("MSE/MSE_", i, ".rds")))
OM_names <- paste0("(", 1:6, ") ", c("Base_h09", "Dome_sel", "Low_h07", "HighM_age3+", "Low_t0", "Low_t0_mat"))

# Performance metrics
## Probability that B > 0.4 B0 (in years 1 - 40)
## Probability that short-term catch > current catch (in years 1 - 10)
## Probability that long-term catch > current catch (in years 11 - 20)

# Update reference
MSE <- lapply(MSE, function(x) {
  require(dplyr)
  # Set BMSY = 0.4 B0
  x@Misc$MSYRefs$ByYear$SSBMSY <- c(0.4 * x@OM$SSB0 - 1e-4) %>% array(c(x@nsim, x@nMPs, x@nyears + x@proyears))
  x@OM$SSBMSY <- x@Misc$MSYRefs$Refs$SSBMSY <- 0.4 * x@OM$SSB0 - 1e-4
  
  x@B_BMSY <- x@SSB/x@OM$SSBMSY
  
  # Set RefY = 20 kt
  x@OM$RefY <- 20000 - 1e-4
  
  # Set FMSY = 1
  x@OM$FMSY <- x@Misc$MSYRefs$Refs$FMSY <- rep(1, x@nsim)
  return(x)
})

`40% B0` <- gfdlm::pm_factory("SBMSY", 1, c(1, 40))
`40% B0` <- DLMtool::P100
`AAVY` <- DLMtool::AAVY
`STC` <- gfdlm::pm_factory("LTY", 1, c(1, 10))
`LTC` <- gfdlm::pm_factory("LTY", 1, c(11, 20))

PM <- c("40% B0", "STC", "LTC", "AAVY")
pm_df_list <- map(MSE, ~ gfdlm::get_probs(.x, PM)) %>% structure(names = OM_names)

# Average performance metric across OMs
pm_avg <- bind_rows(pm_df_list, .id = "scenario")  %>%
  group_by(MP) %>%
  summarise_if(is.numeric, mean)



#### Tradeoff plots
plot_tradeoff(pm_df_list, "40% B0", "STC") + ggrepel::geom_text_repel(aes(label = MP)) + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff_allMPs.png", height = 6, width = 8)

plot_tradeoff(pm_df_list, "LTC", "STC") + ggrepel::geom_text_repel(aes(label = MP)) + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff_catch_allMPs.png", height = 6, width = 8)


# Remove MPs that have three year interval updates
MP_subset <- pm_df_list[[1]]$MP[!grepl("3u", pm_df_list[[1]]$MP)]
plot_tradeoff(pm_df_list, "40% B0", "STC", mp = MP_subset) + 
  ggrepel::geom_text_repel(aes(label = MP)) + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff.png", height = 6, width = 8)

plot_tradeoff(pm_df_list, "LTC", "STC", mp = MP_subset) + ggrepel::geom_text_repel(aes(label = MP)) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff_catch.png", height = 6, width = 8)


plot_tradeoff(list(Averaged = pm_avg), "40% B0", "STC", mp = MP_subset) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggrepel::geom_text_repel(aes(label = MP)) + theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff_avg.png", height = 3, width = 3)

plot_tradeoff(list(Averaged = pm_avg), "40% B0", "STC") + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggrepel::geom_text_repel(aes(label = MP)) + theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
ggsave("figures/tradeoff_all.png", height = 3, width = 3)


# Plot blank tradeoff
guide <- expand.grid(x = c("High", "Low"), y = "Biomass", xx = c("High", "Low"), yy = "Catch")
label <- mapply(function(x, y) paste(x, y, sep = "\n"), 
                x = apply(guide[, 1:2], 1, paste0, collapse = " "), 
                y = apply(guide[, 3:4], 1, paste0, collapse = " "))
guide <- data.frame(label = label, x = c(1, 0, 1, 0), y = c(1, 1, 0, 0),
                    hjust = c("right", "left", "right", "left"),
                    vjust = c("top", "top", "bottom", "bottom"))

ggplot(data.frame(x = -1, y = -1), aes(x, y)) + coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_blank() + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_label(data = guide, aes(x, y, label = label, hjust = hjust, vjust = vjust)) + 
  gfdlm::theme_pbs() + labs(x = "40% B0", y = "STC") + 
  theme(panel.grid.major.y = ggplot2::element_line(colour = "grey85"), 
        panel.grid.major.x = ggplot2::element_line(colour = "grey85"))
ggsave("figures/tradeoff_blank.png", height = 3, width = 3)




#### Tigure is a figure that serves as a table. Color scheme changed so that green represents high probability.
plot_tigure_facet(pm_df_list) + 
  scale_fill_gradient2(limits = c(0, 1), high = "#31A354", low = "#FA9FB5", midpoint = 0.5)
ggsave("figures/table_allMPs.png", height = 6, width = 6)

plot_tigure_facet(lapply(pm_df_list, filter, MP %in% MP_subset)) + 
  scale_fill_gradient2(limits = c(0, 1), high = "#31A354", low = "#FA9FB5", midpoint = 0.5)
ggsave("figures/table.png", height = 6, width = 6)


plot_tigure(pm_avg) + 
  scale_fill_gradient2(limits = c(0, 1), high = "#31A354", low = "#FA9FB5", midpoint = 0.5)
ggsave("figures/table_all_avg.png", height = 3, width = 3)

plot_tigure(filter(pm_avg, MP %in% MP_subset)) + 
  scale_fill_gradient2(limits = c(0, 1), high = "#31A354", low = "#FA9FB5", midpoint = 0.5)
ggsave("figures/table_avg.png", height = 3, width = 3)




#### Radar plots
plot_radar_facet(lapply(pm_df_list, filter, MP %in% MP_subset & MP != "NFref"))
ggsave("figures/radar_all.png", height = 6, width = 8)

plot_radar(filter(pm_avg, MP %in% MP_subset & MP != "NFref"))
ggsave("figures/radar_avg.png", height = 4, width = 6)

# Radar plot for 3 MPs
plot_radar(pm_avg[c(2, 7, 12), ])
ggsave("figures/radar_3MPs.png", height = 4, width = 6)




#### Projections that show: B/B0 and catch 
source("05a_projection_plots.R")
plot_projection_SKJ(MSE[[1]], MP = c("NFref", "curE", "spict_4010_", "CC_15kt", "CC_20kt", "CC_30kt"), 
                    rel_widths = c(0.9, 1))
ggsave("figures/projection_OM1_part1.png", height = 8, width = 6)

plot_projection_SKJ(MSE[[1]], MP = c("GBslope", "Iratio_", "Islope"), 
                    rel_widths = c(0.9, 1))
ggsave("figures/projection_OM1_part2.png", height = 4, width = 6)

plot_projection_SKJ(MSE[[5]], MP = c("NFref", "curE", "spict_4010_", "CC_15kt", "CC_20kt", "CC_30kt"), 
                    rel_widths = c(0.9, 1))
ggsave("figures/projection_OM5_part1.png", height = 8, width = 6)

plot_projection_SKJ(MSE[[5]], MP = c("GBslope", "Iratio_", "Islope"), 
                    rel_widths = c(0.9, 1))
ggsave("figures/projection_OM5_part2.png", height = 4, width = 6)

