
library(MSEtool)
library(dplyr)
library(ggplot2)
library(purrr)

#devtools::install_github("pbs-assess/gfdlm")
library(gfdlm)

SRA <- readRDS("OM/SRA_07.27.2020.rds")
OM_names <- paste0("(", 1:6, ") ", c("Base_h09", "Dome_sel", "Low_h07", "HighM_age3+", "Low_t0", "Low_t0_mat"))

# Plot catch
png("figures/catch.png", width = 4.5, height = 3.5, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
matplot(seq(SRA[[1]]@OM@CurrentYr - SRA[[1]]@OM@nyears + 1, SRA[[1]]@OM@CurrentYr), 
        SRA[[1]]@data$Chist/1e3, typ = 'l', xlab = "Year", ylab = "Catch (kt)", lty = 1)
legend("topleft", c("Baitboat", "Handline"), col = 1:2, lty = 1)
abline(h = 0, col = "grey")
dev.off()

# Plot length comps
data <- readRDS("OM/SKJ_data_07.27.2020.rds")
png("figures/CAL.png", width = 4.5, height = 3.5, units = "in", res = 400)
MSEtool:::plot_composition_SRA(2003:2018, SRA = data$CAL[46:61, , ] %>% aperm(c(3, 1, 2)),
                               CAL_bins = data$length_bin, N = NULL, dat_color = 1:2)
dev.off()


# Plot selectivity vs. maturity vs. length-at-age
get_ogives <- function(x, scenario) {
  out <- lapply(c("V", "Len_age", "Mat_age"), function(xx) getElement(x@OM@cpars, xx)[1, , x@OM@nyears])
  out[[2]] <- out[[2]]/x@OM@Linf[1]
  if(is.null(out[[3]])) out[[3]] <- x@mean_fit$obj$env$data$mat[1,]
  data.frame(Value = do.call(c, out), Age = rep(1:x@OM@maxage, 3), 
             Type = rep(c("Selectivity", "Rel length-at-age", "Maturity"), each = x@OM@maxage),
             scenario = scenario)
}
g <- purrr::map2_df(SRA, OM_names, get_ogives) %>%
  ggplot(aes(Age, Value, linetype = Type, shape = Type)) + geom_line() + geom_point() + 
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Age", y = "Value") + coord_cartesian(expand = FALSE, ylim = c(0, 1.1))
ggsave("figures/ogives.png", width = 8, height = 5)

# Plot M
get_M <- function(x, scenario) {
  data.frame(M = x@OM@cpars$M_ageArray[1, , 1], Age = 1:x@OM@maxage, scenario = scenario)
}

g <- purrr::map2_df(SRA, OM_names, get_M) %>%
  ggplot(aes(Age, M)) + geom_line() + geom_point() + 
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Age", y = "Natural Mortality") + coord_cartesian(expand = FALSE, ylim = c(0, 2))
ggsave("figures/M.png", width = 8, height = 5)

# Plot yield curve

get_yield_curve <- function(x, scenario) {
  F.vector = seq(0, 10, length.out = 5e2) # see MSEtool:::plot_yield_SCA
  
  M <- x@OM@cpars$M_ageArray[1,,x@OM@nyears]
  mat <- x@mean_fit$obj$env$data$mat[x@OM@nyears, ]
  weight <- x@OM@cpars$Wt_age[1,,x@OM@nyears]
  maxage <- x@OM@maxage
  SR <- "BH"
  
  vul <- x@OM@cpars$V[1, , x@OM@nyears]
  
  SSB0 <- vapply(x@Misc, getElement, numeric(1), "E0_SR")
  
  Arec <- vapply(x@Misc, getElement, numeric(1), "Arec")
  Brec <- vapply(x@Misc, getElement, numeric(1), "Brec")
  
  solveMSY <- function(logF, Arec, Brec, pars = FALSE) {
    Fmort <- exp(logF)
    surv <- exp(-vul * Fmort - M)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <- sum(NPR * mat * weight)
    if(SR == "BH") Req <- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * Fmort/(vul * Fmort + M) * NPR * (1 - exp(-vul * Fmort - M))
    Yield <- Req * sum(CPR * weight)
    if(pars) {
      return(c(Yield = Yield, EPR = EPR, Req = Req))
    } else return(-1 * Yield)
  }
  # Calc yield curve for each Brec (indexed by i)
  yield_curve_sim <- function(i) {
    out <- lapply(log(F.vector), solveMSY, Arec = Arec[i], Brec = Brec[i], pars = TRUE)    
    data.frame(Yield = sapply(out, getElement, 1), F = F.vector, dep = sapply(out, function(xx) prod(xx[2:3]))/SSB0[i],
               iteration = i)
  }
  
  out <- do.call(rbind, lapply(1:x@OM@nsim, yield_curve_sim))
  out$scenario <- scenario
  return(out)
}
SKJ_yield_curve <- purrr::map2_df(SRA, OM_names, get_yield_curve)
g <- SKJ_yield_curve %>% 
  ggplot(aes(F, Yield, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Apical F", y = "Yield")
ggsave("figures/yield_curve_F.png", width = 8, height = 5)

g <- SKJ_yield_curve %>% 
  ggplot(aes(dep, Yield, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Spawning depletion", y = "Yield")
ggsave("figures/yield_curve_dep.png", width = 8, height = 5)


# Plot SSB
get_SSB <- function(x, scenario) {
  all_years <- (x@OM@CurrentYr - x@OM@nyears + 1):(x@OM@CurrentYr + 1)
  reshape2::melt(x@SSB) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}
g <- purrr::map2_df(SRA, OM_names, get_SSB) %>%
  #mutate(scenario = factor(scenario, levels = OM_names)) %>%
  ggplot(aes(year, value, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Spawning biomass") + coord_cartesian(expand = FALSE, ylim = c(0, 310000))
ggsave("figures/SSB.png", width = 8, height = 5)

# Plot depletion
get_depletion <- function(x, scenario) {
  depletion <- x@SSB / vapply(x@Misc, getElement, numeric(1), "E0_SR")
  all_years <- (x@OM@CurrentYr - x@OM@nyears + 1):(x@OM@CurrentYr + 1)
  reshape2::melt(depletion) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}
g <- purrr::map2_df(SRA, OM_names, get_depletion) %>%
  #mutate(scenario = factor(scenario, levels = OM_names)) %>%
  ggplot(aes(year, value, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Spawning depletion") + coord_cartesian(expand = FALSE, ylim = c(0, 1.3))
ggsave("figures/depletion.png", width = 8, height = 5)

# Plot F
get_F2 <- function(x, scenario) {
  .F1 <- purrr::map(x@Misc, "F_at_age")
  .F <- purrr::map_dfc(.F1, ~tibble(.F = apply(.x, 1, max)))
  .F <- t(as.matrix(.F))
  row.names(.F) <- NULL
  
  last_year <- dim(.F)[2]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr)
  
  reshape2::melt(.F) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}
g <- purrr::map2_df(SRA, OM_names, get_F2) %>%
  #mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, value, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Apical F") +
  coord_cartesian(ylim = c(0, 10), expand = FALSE)
ggsave("figures/apical_F.png", width = 8, height = 5)


# Plot exploitation rate: ratio of total catch to vulnerable biomass at the start of the year
# The mean abundance during the year is much lower due to high M
get_U <- function(x, scenario) {
  U <- lapply(x@Misc, function(x) rowSums(x$Cpred))
  
  V <- x@OM@cpars$V[, , 1:(x@OM@nyears+1)]
  N <- lapply(x@Misc, getElement, "N")
  Wt_age <- x@OM@cpars$Wt_age[, , 1:(x@OM@nyears+1)]
  
  VB <- lapply(1:length(N), function(xx) colSums((N[[xx]] %>% t()) * V[xx, , ] * Wt_age[xx, , ]))
  U <- Map(function(xx, yy) rowSums(xx$Cpred)/yy[-length(yy)], xx = x@Misc, yy = VB)
  U <- do.call(rbind, U)
  
  last_year <- dim(U)[2]
  all_years <- seq(x@OM@CurrentYr - x@OM@nyears + 1, x@OM@CurrentYr)
  
  reshape2::melt(U) %>%
    rename(iteration = Var1) %>%
    mutate(year = rep(all_years, each = max(iteration))) %>%
    mutate(scenario = scenario)
}
g <- purrr::map2_df(SRA, OM_names, get_U) %>%
  #mutate(scenario = factor(scenario, levels = sc$scenario_human)) %>%
  ggplot(aes(year, value, group = iteration)) +
  geom_line(alpha = 0.05) +
  facet_wrap(vars(scenario)) +
  gfplot::theme_pbs() +
  labs(x = "Year", y = "Exploitation Rate") +
  coord_cartesian(ylim = c(0, 1), expand = FALSE)
ggsave("figures/exploitation.png", width = 8, height = 5)



# Plot fit to baitboat CPUE
get_survey <- function(sra, sc_name, survey_names) {
  n_surv <- dim(sra@Misc[[1]]$Ipred)[2]
  out2 <- purrr::map(seq_len(n_surv), function(i) {
    surveys <- do.call(cbind, purrr::map(sra@Misc, ~ .$Ipred[,i,drop=FALSE]))
    out <- reshape2::melt(surveys) %>%
      rename(year = Var1, iter = Var2)
    out$year <- out$year
    out$scenario <- sc_name
    out$survey <- survey_names[i]
    out
  })
  bind_rows(out2)
}
surv <- purrr::map2_dfr(SRA, OM_names, get_survey, survey_names = "Baitboat CPUE")
all_years <- (SRA[[1]]@OM@CurrentYr - SRA[[1]]@OM@nyears + 1):SRA[[1]]@OM@CurrentYr
surv$year <- surv$year + min(all_years) - 1

I_sd <- SRA[[1]]@data$I_sd %>% structure(dimnames = list(all_years, "Baitboat CPUE")) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey", value.name = "SD")

Index <- SRA[[1]]@data$Index %>% structure(dimnames = list(all_years, "Baitboat CPUE")) %>%
  as.data.frame() %>% cbind(data.frame(Year = all_years)) %>%
  reshape2::melt(id.vars = c("Year"), variable.name = "survey") %>% left_join(I_sd, by = c("Year", "survey"))

Index$lower <- exp(log(Index$value) - 2 * Index$SD)
Index$upper <- exp(log(Index$value) + 2 * Index$SD)

# Plot Index
g <- ggplot(surv, aes(year, value, group = paste(iter))) +
  geom_pointrange(data = Index, mapping = aes(x = Year, y = value, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, pch = 16, colour = "red") +
  geom_line(alpha = 0.05) + facet_wrap(~scenario) + gfplot::theme_pbs() +
  ylab("Baitboat CPUE") + xlab("Year") 
ggsave("figures/CPUE_fit.png", width = 8, height = 5)

# Plot predicted length comps
SRA_pred <- lapply(SRA[[1]]@Misc, function(x) x$CALpred[, , 1]) %>% unlist() %>% 
  array(c(61, 41, 100)) %>% aperm(c(3, 1, 2))

png("figures/CALpred.png", width = 4.5, height = 3.5, units = "in", res = 400)
MSEtool:::plot_composition_SRA(2003:2018, SRA = SRA_pred[, 46:61, ], dat = SRA[[1]]@data$CAL[46:61, , 1],
                               CAL_bins = SRA[[1]]@data$length_bin, N = NULL, dat_color = 4)
dev.off()


SRA_pred <- lapply(SRA[[1]]@Misc, function(x) x$CALpred[, , 2]) %>% unlist() %>% 
  array(c(61, 41, 100)) %>% aperm(c(3, 1, 2))
png("figures/CALpred_HL.png", width = 4.5, height = 3.5, units = "in", res = 400)
MSEtool:::plot_composition_SRA(2012:2018, SRA = SRA_pred[, 55:61, ], dat = SRA[[1]]@data$CAL[55:61, , 2],
                               CAL_bins = SRA[[1]]@data$length_bin, N = NULL, dat_color = 4)
dev.off()