
library(MSEtool)
library(dplyr)
library(DLMextra)

source("03a_MPs.R")

setup(8)
sfExportAll()
sfLibrary(DLMextra)


###### runMSE
SRA <- readRDS("OM/SRA_07.27.2020.rds")

# 40-year projection period
for(i in 1:length(SRA)) {
  OM <- SRA[[i]]@OM
  OM@cpars$Data@Cat <- matrix(NA, 1, 1)

  OM@interval <- c(50, 50, 1, 3)
  MSE <- runMSE(OM, MPs = c("NFref", "curE", "spict_4010_", "spict_4010_3u"), parallel = TRUE, save_name = paste0("rawMSE_", i, "_spict.rds"))
  message("MSE ", i, " (spict) complete.")
  saveRDS(MSE, file = paste0("MSE_", i, "_spict.rds"))

  #OM@interval <- c(50, 50, 50, 1, 3, 1, 3, 1, 3)
  #MSE2 <- runMSE(OM, MPs = c("CC_15t", "CC_20t", "CC_30t", "GBslope", "GBslope_3u",
  #                           "Iratio_", "Iratio_3u", "Islope", "Islope_3u"), parallel = TRUE,
  #               save_name = paste0("rawMSE_", i, "_empirical.rds"))
  #message("MSE ", i, " (empirical) complete.")
  #saveRDS(MSE2, file = paste0("MSE_", i, "_empirical.rds"))
}


merge_MSE <- function(...) {
  dots <- list(...)
  if(length(dots) == 1) dots <- dots[[1]]
  
  slots_identical <- function(slotname, x = dots, is_logical = FALSE) {
    res <- lapply(x, getElement, slotname)
    is_identical <- all(vapply(res[-1], identical, logical(1), res[[1]]))
    if(is_logical) {
      return(is_identical)
    } else return(unique(do.call(c, res)))
  }
  
  slots_identical("Name")
  slots_identical("nyears")
  slots_identical("proyears")
  slots_identical("nsim")
  
  stopifnot(slots_identical("OM", is_logical = TRUE))
  stopifnot(slots_identical("Obs", is_logical = TRUE))
  stopifnot(slots_identical("SSB_hist", is_logical = TRUE))
  stopifnot(slots_identical("CB_hist", is_logical = TRUE))
  stopifnot(slots_identical("FM_hist", is_logical = TRUE))
  
  #nMPs <- vapply(dots, getElement, numeric(1), "nMPs")
  
  slotvec <- c("B_BMSY", "F_FMSY", "B", "SSB", "VB", "FM", "C", "TAC", "Effort", "PAA", "CAA", "CAL")
  res <- lapply(slotvec, function(x) do.call(abind::abind, c(lapply(dots, getElement, x), along = 2)))
  names(res) <- slotvec
  
  Misc <- lapply(dots, slot, "Misc")
  #names(Misc[[1]])
  
  Misc_identical <- function(x) all(vapply(x[-1], identical, logical(1), x[[1]]))
  
  Data <- do.call(c, lapply(Misc, getElement, "Data"))
  
  TryMP <- do.call(c, lapply(Misc, getElement, "TryMP"))
  
  Unfished <- lapply(Misc, getElement, "Unfished")
  Unfished_Refs <- lapply(Unfished, getElement, "Refs")
  stopifnot(Misc_identical(Unfished_Refs))
  
  Unfished_ByYear <- lapply(Unfished, getElement, "ByYear")
  stopifnot(Misc_identical(Unfished_ByYear))
  
  MSYRefs <- lapply(Misc, getElement, "MSYRefs")
  MSYRefs_Refs <- lapply(MSYRefs, getElement, "Refs")
  stopifnot(Misc_identical(MSYRefs_Refs))
  
  MSYRefs_ByYear <- lapply(MSYRefs, getElement, "ByYear")
  MSYRefs_ByYear2 <- lapply(names(MSYRefs_ByYear[[1]]),
                            function(x) do.call(abind::abind, c(lapply(MSYRefs_ByYear, getElement, x), along = 2)))
  names(MSYRefs_ByYear2) <- names(MSYRefs_ByYear[[1]])
  
  Misc_new <- list(Data = Data, TryMP = TryMP,
                   Unfished = list(Refs = Unfished_Refs[[1]], ByYear = Unfished_ByYear[[1]]),
                   MSYRefs = list(Refs = MSYRefs_Refs[[1]], ByYear = MSYRefs_ByYear2))
  
  slotvec_Misc <- c("LatEffort", "Revenue", "Cost", "TAE")
  Misc2 <- lapply(slotvec_Misc, function(x) do.call(abind::abind, c(lapply(Misc, getElement, x), along = 2)))
  names(Misc2) <- slotvec_Misc
  
  ## Create MSE Object ####
  MSEout <- new("MSE", Name = slots_identical("Name"), nyears = slots_identical("nyears"),
                proyears = slots_identical("proyears"), nMPs = length(slots_identical("MPs")),
                MPs = slots_identical("MPs"), nsim = slots_identical("nsim"),
                OM = dots[[1]]@OM, Obs = dots[[1]]@Obs, B_BMSY = res$B_BMSY, F_FMSY = res$F_FMSY, B = res$B, SSB = res$SSB,
                VB = res$VB, FM = res$FM, res$C, TAC = res$TAC, SSB_hist = dots[[1]]@SSB_hist, CB_hist = dots[[1]]@CB_hist,
                FM_hist = dots[[1]]@FM_hist, Effort = res$Effort, PAA = res$PAA, CAA = res$CAA, CAL = res$CAL,
                CALbins = slots_identical("CALbins"), Misc = c(Misc_new, Misc2))
  
  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version
  
  MSEout
}

for(i in 1:length(SRA)) {
  MSE <- readRDS(paste0("MSE/MSE_", i, "_spict.rds"))
  MSE2 <- readRDS(paste0("MSE/MSE_", i, "_empirical.rds"))
  out <- merge_MSE(MSE, MSE2)
  saveRDS(out, file = paste0("MSE/MSE_", i, ".rds"))
}
