library(MSEtool)
library(dplyr)

####### Data generated from 01_generate_data.R
data <- readRDS("OM/SKJ_data_07.27.2020.rds")

# Downweight the length composition by re-scaling the sample size of the length composition. Use the log(sample size)
# as the multinomial sample size.
rescale_CAL <- function(x) {
  if(all(is.na(x))) {
    return(x)
  } else {
    N <- sum(x, na.rm = TRUE)
    return(x * log(N)/ N)
  }
}
for(i in 1:dim(data$CAL)[3]) {
  for(y in 1:dim(data$CAL)[1]) {
    data$CAL[y, , i] <- rescale_CAL(data$CAL[y, , i])
  }
}



Year <- c(1958:2018)

#### Stock
stock_SKJ <- new("Stock")
stock_SKJ@Name <- "WATL SKJ"
stock_SKJ@Common_Name <- "Western Atlantic Skipjack Tuna"
stock_SKJ@Species <- "Katsuwonus pelamis"
stock_SKJ@maxage <- 6

stock_SKJ@R0 <- 50000 # Start value of R0 for SRA_scope, will be estimated by conditioning

stock_SKJ@SRrel <- 1
stock_SKJ@Perr <- c(0.4, 0.4) # Arbitrary value for now

##comments on VB model parameters
##new available information Soares et al 2019
##there is also another paper (more recent) from Dr. Castello which shows 
##changes in population structure and growth
##so for now my suggestion is to use Garbin & Castello (2014), as base case OM

#stock_SKJ@Linf <- c(92.4, 92.4) # Garbin & Castello (2014)
#stock_SKJ@K <- c(0.161, 0.161)  # Garbin & Castello (2014)
#stock_SKJ@t0 <- c(-2.9, -2.9)   # Garbin & Castello (2014)
stock_SKJ@Linf <- c(90.1, 90.1)   # Soares et al 2019
stock_SKJ@K <- c(0.24, 0.24)    # Soares et al 2019
stock_SKJ@t0 <- c(-0.54, -0.54)   # Soares et al 2019
stock_SKJ@LenCV <- c(0.10, 0.10)  # Arbitrary value for now, population variance not reported in either paper

stock_SKJ@a <- 0.004  # Soares et al 2019
stock_SKJ@b <- 3.4184 # Soares et al 2019 

get_M_at_age <- function(Linf = stock_SKJ@Linf, K = stock_SKJ@K, t0 = stock_SKJ@t0) {
  Length_at_age <- Linf * (1 - exp(-K * (c(1:stock_SKJ@maxage) - t0)))
  
  M_at_age <- 12.01 * exp(-0.08*Length_at_age + 0.0005 * Length_at_age^2) #2014 assessment
  ifelse(Length_at_age < 15, M_at_age + 1.77, M_at_age)
}


stock_SKJ@L50_95 <- c(7, 7) # my guess for Length increment from 50 percent to 95 percent maturity
#stock_SKJ@L50 <- c(51,51) # Vilella & Castello (1993): alternative 
stock_SKJ@L50 <- c(43.2,43.2) # Soares et al 2019 for females


# Comments on SDs;why this range is c(0,0)??? (no variability?)
# Follow paper from Garbin & Castello (2014) which shows changes in population structure and growth
# it might be better to allow some inter-annual variation
# also given these parameter are very poorly quantified, specially M, I replaced this slot for all parameters  
# to c (0,0.1) to allow variability in of up 10%
# not sure if is more appropriate to consider this in the slot "parameter@grad", which my understanding
# represent long-term trends on these parameters
stock_SKJ@M <- stock_SKJ@M2 <- get_M_at_age()
stock_SKJ@Linfsd <- stock_SKJ@Ksd <- stock_SKJ@Msd <- c(0, 0)

stock_SKJ@h <- c(0.9, 0.9) 
  
# 100% discard mortality but we do not have discards in the OM (all removals in catch)
stock_SKJ@Fdisc <- c(1, 1)

# Will be updated by SRA_scope, need a placeholder for now
stock_SKJ@AC <- stock_SKJ@D <- c(0, 0)
# I don't know how SRA_scope works, but I just put the output of DBSRA preliminary run (B2018/K) 
#stock_SKJ@D <- c(0.48,0.48)

# Effectively parameterize a single area model
# These parameters at ("0.5") allow some exchanges (migrations) between "W and E stocks" ?
stock_SKJ@Frac_area_1 <- stock_SKJ@Size_area_1 <- stock_SKJ@Prob_staying <- c(0.5, 0.5)


#### Fleet
fleet_SKJ <- new("Fleet")
fleet_SKJ@nyears <- length(Year) #101
fleet_SKJ@Spat_targ <- c(1, 1) # No spatial dynamics
fleet_SKJ@DR <- c(0, 0) # All removals will be accounted for in the catch
fleet_SKJ@CurrentYr <- max(Year)

# Placeholder, not used. Effort/F/selectivity will be updated by SRA_scope, qinc/qcv only affect effort MPs
#fleet_SKJ@EffYears <- c(1, fleet_SKJ@nyears)
fleet_SKJ@qinc <- fleet_SKJ@qcv <- c(0, 0)
#fleet_SKJ@Esd <- fleet_SKJ@qinc <- fleet_SKJ@qcv <- c(0, 0)

# Selectivity parameters will be starting values for HBLL survey selectivity (fixed for all other fleets/surveys)
fleet_SKJ@L5 <- c(40, 40)
fleet_SKJ@LFS <- c(55, 55)
fleet_SKJ@Vmaxlen <- c(0.5, 0.5)
fleet_SKJ@isRel <- FALSE

#### Observation and Implementation
obs_SKJ <- DLMtool::Precise_Unbiased
imp_SKJ <- DLMtool::Perfect_Imp

OM <- new("OM", stock_SKJ, fleet_SKJ, obs_SKJ, imp_SKJ)
OM@Name <- "WATL SKJ"
OM@proyears <- 40

#### cpars
OM@nsim <- 100


# Run SRA for base model (high steepness)
map_log_rec_dev <- rep(NA, 16) %>% c(1:45)

#vul_par_high <- matrix(c(55, 45, 1), 3, 2)
#vul_par_low <- matrix(c(40, 35, 1), 3, 2)
#map_vul_par <- matrix(NA, 3, 2)
#map_vul_par <- matrix(c(1, 2, NA), 3, 2)

SRA <- list()

SRA_wrap <- function(...) {
  dots <- list(...)
  dots$max_F <- 10
  res <- do.call(SRA_scope, dots)
  dots$vul_par <- rbind(res@mean_fit$report$LFS, res@mean_fit$report$L5, res@mean_fit$report$Vmaxlen)
  dots$map_vul_par <- array(NA, dim(dots$vul_par))
  dots$resample <- TRUE
  do.call(SRA_scope, dots)
}

# High h, logistic sel
SRA[[1]] <- SRA_wrap(OM = OM, data = data, condition = "catch2", 
                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)

# High h, dome sel
SRA[[2]] <- SRA_wrap(OM = OM, data = data, condition = "catch2", 
                     selectivity = rep("dome", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)

# Low h, logistic sel
OM_low_h <- OM
OM_low_h@h <- c(0.7, 0.7)

SRA[[3]] <- SRA_wrap(OM = OM_low_h, data = data, condition = "catch2", 
                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)

# High M - age 3+ = 0.75
OM_high_M <- OM
OM_high_M@M[3:6] <- OM_high_M@M2[3:6] <- 0.75
SRA[[4]] <- SRA_wrap(OM = OM_high_M, data = data, condition = "catch2", 
                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)


# Lower t0, same maturity
OM_low_t0 <- OM
OM_low_t0@t0 <- c(-2, -2)
OM_low_t0@cpars$Mat_age <- SRA[[1]]@OM@cpars$Mat_age
SRA[[5]] <- SRA_wrap(OM_low_t0, data, condition = "catch2", 
                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)

# Lower t0, low maturity
OM_low2 <- OM
OM_low2@t0 <- c(-2, -2)
SRA[[6]] <- SRA_wrap(OM_low2, data, condition = "catch2", 
                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)


# 2 % increase in catchability of CPUE
# Project with hyperstable indices?

# The CPUE time series is very short so not included in set of OMs.

#data_qinc <- data
#data_qinc$Index <- data$Index / 1.02 ^ (c(1:OM@nyears)-1)
#data_qinc$Index <- data_qinc$Index / mean(data_qinc$Index, na.rm = TRUE) * mean(data$Index, na.rm = TRUE)
#
#SRA[[7]] <- SRA_wrap(OM = OM, data_qinc, condition = "catch2", 
#                     selectivity = rep("logistic", 2), s_selectivity = 1, map_log_rec_dev = map_log_rec_dev)



# Save SRA file
saveRDS(SRA, file = "OM/SRA_07.27.2020.rds")

# Save reports
scenario_names <- c("Base_h09", "Dome_sel", "Low_h07", "HighM_age3+", "Low_t0", "Low_t0_mat")
fleet_name <- c("Baitboat", "Handline")
survey_name <- "Baitboat_CPUE"

compare_SRA(SRA[[1]], SRA[[2]], SRA[[3]], SRA[[4]], SRA[[5]], SRA[[6]], dir = getwd(), 
            filename = "compare_SRA_07.27.2020",
            f_name = fleet_name, s_name = survey_name, 
            scenario = list(lty = 1, names = scenario_names))


lapply(1:length(SRA), function(x) plot(SRA[[x]], dir = getwd(), filename = paste0("SRA_", x),
                                       open_file = FALSE, f_name = fleet_name, s_name = survey_name))


# Surplus production model with r intrinsic prior
Data <- SRA[[1]]@OM@cpars$Data
Data@Mort <- 1
Data@CV_Mort <- 0.15

Data@steep <- 0.9
Data@CV_steep <- 0.1
Data@MaxAge <- 6

Data@wla <- mean(OM@a)
Data@wlb <- mean(OM@b)
Data@L50 <- mean(OM@L50)
Data@L95 <- mean(OM@L50 + OM@L50_95)

Data@vbLinf <- mean(OM@Linf)
Data@vbK <- mean(OM@K)
Data@vbt0 <- mean(OM@t0)

rp <- MSEtool:::r_prior_fn(1, Data, r_reps = 10000, SR_type = "BH")
hist(rp)
mean(rp)
sd(rp)

SP_mod <- SP(Data = Data, AddInd = 1, use_r_prior = TRUE, start = list(r_prior = c(1.33, 0.45)))
SP_mod <- SP_SS(Data = Data, AddInd = 1, use_r_prior = TRUE, start = list(n = 1, r_prior = c(1.33, 0.45)))
SP_mod <- SP_SS(Data = Data, AddInd = 1)

Data@Ind <- Data@AddInd[, 1, ]
Data@CV_Ind <- Data@CV_AddInd[, 1, ]

DLMextra()
res <- DLMextra::spict(Data = Data, fix_alpha = FALSE, fix_beta = FALSE)

