# Generate catch

library(dplyr)
library(reshape2)
library(readr)
library(ggplot2)

catch <- readxl::read_excel("raw_data/t1nc-ALL_20200115.xlsx", range = "A4:U68290") %>% 
  filter(Species == "SKJ" & Stock == "ATW")


######### Explore the data
# Summarize by country and gear
#country_gear <- catch %>% 
#  group_by(YearC, Flag, GearGrp) %>% 
#  summarise(Catch = sum(Qty_t)) 
#ggplot(country_gear, aes(YearC, Catch, colour = Flag)) + facet_wrap(~ GearGrp) + geom_line()

# Summarize by country
#country <- catch %>% 
#  group_by(YearC, Flag) %>% 
#  summarise(Catch = sum(Qty_t)) 
#ggplot(country, aes(YearC, Catch, colour = Flag)) + facet_wrap(~ Flag) + geom_line()
#
## Summarize by SampAreaCode and Area columns
#areas <- catch %>% group_by(YearC, GearGrp, SampAreaCode, Area) %>% 
#  summarise(Catch = sum(Qty_t)) 
#ggplot(areas, aes(YearC, Catch)) + facet_grid(SampAreaCode ~ Area) + geom_line()
#
#areas <- catch %>% filter(Flag == "Brazil") %>% group_by(YearC, GearGrp, Area) %>% 
#  summarise(Catch = sum(Qty_t)) 
#ggplot(areas, aes(YearC, Catch)) + facet_grid(GearGrp ~ Area) + geom_line()
#
#areas <- catch %>% filter(Flag == "Venezuela") %>% group_by(YearC, GearGrp, Area) %>% 
#  summarise(Catch = sum(Qty_t)) 
#ggplot(areas, aes(YearC, Catch)) + facet_grid(GearGrp ~ Area) + geom_line()


######### Generate catch for fleets
#1. Brazil BB
#2. Brazil HL

Year <- data.frame(YearC = 1958:2018)
catch_matrix <- matrix(NA, nrow = length(Year), ncol = 4)

bra_bb <- filter(catch, Area == "SW" & Flag == "Brazil" & GearGrp == "BB") %>% group_by(YearC) %>%
  summarise(Catch = sum(Qty_t)) %>% right_join(Year)
bra_hl <- filter(catch, Area == "SW" & Flag == "Brazil" & GearGrp == "HL") %>% group_by(YearC) %>%
  summarise(Catch = sum(Qty_t)) %>% right_join(Year)

catch_matrix <- cbind(bra_bb$Catch, bra_hl$Catch)

# Catch substitutions for missing catch data
catch_matrix[2:3, 1] <- 200
catch_matrix[15:18, 1] <- 100

catch_matrix[1:49, 2] <- 0.01
catch_matrix[51, 2] <- 4.878
catch_matrix[61, 2] <- 5300

#matplot(Year, catch_matrix, typ = 'l')

######### Generate Index
Ind <- read_csv("raw_data/index_ver00.csv") %>% right_join(data.frame(Year = Year$YearC), by = c("year" = "Year"))

Index <- Ind$predicted.average
SD <- Ind$predicted.sd


######### Length comps
cas <- read.csv("raw_data/t2sz_20131210_SKJ.csv") %>% filter(FreqTypeCod == "FL")

# Use 2017-2018 for Brazil BB
cas2 <- read_csv("raw_data/size_SKJ_BRA_BB_HL_2012_2018.csv")

# Use 2010 - 2013 for Brazil BB
cas3 <- readxl::read_excel("raw_data/BRA_SC_UNIVALI-EMCT-LEMA_Skipjack_Length_1995-2013.xlsx") %>% 
  reshape2::melt(id.var = "Length.cm")
cas3$Gear <- "BB"
cas3$YearC <- as.numeric(cas3[, 2]) + 1994
cas3 <- filter(cas3, YearC >= 2010)
cas3 <- cas3[, c(5, 4, 1, 3)]

names(cas2) <- names(cas3) <- c("YearC", "Gear", "ClassFrq", "N")

bra_bb <- filter(cas, Flag == "Brasil" & GearGrpCod == "BB") %>% group_by(YearC, ClassFrq) %>% 
  summarise(N = sum(Nr)) %>% rbind(filter(cas2, Gear == "BB")) %>% rbind(cas3)
#ggplot(bra_bb, aes(ClassFrq, N)) + facet_wrap(~YearC, scales = "free_y") + geom_line()

bra_hl <- filter(cas2, Gear == "HL")
#ggplot(bra_hl, aes(ClassFrq, N)) + facet_wrap(~YearC, scales = "free_y") + geom_line()





bins <- seq(30, 110, 2) # changed to accomadate new Linf and Linfsd (10%)

len <- list(bra_bb, bra_hl) %>% 
  lapply(function(x, Year, bins) {
    xx <- acast(x, list("YearC", "ClassFrq"), value.var = "N", fill = 0)
    yr_ind <- vapply(rownames(xx) %>% as.numeric(), match, numeric(1), table = Year)
    
    output <- matrix(NA, nrow = length(Year), ncol = length(bins))
    b_ind <- 2:(ncol(output)-1)
    
    minus_group <- as.numeric(colnames(xx)) <= bins[1]
    if(any(minus_group)) {
      output[yr_ind, 1] <- rowSums(xx[, minus_group, drop = FALSE], na.rm = TRUE)
    } else {
      output[yr_ind, 1] <- 0
    }
    
    plus_group <- as.numeric(colnames(xx)) >= bins[length(bins)]
    if(any(plus_group)) {
      output[yr_ind, length(bins)] <- rowSums(xx[, plus_group, drop = FALSE], na.rm = TRUE)
    } else {
      output[yr_ind, length(bins)] <- 0
    }
    
    match_columns <- lapply(bins[b_ind], function(y) {
      z <- match(y, as.numeric(colnames(xx)))
      if(is.na(z)) {
        return(rep(0, length(yr_ind)))
      } else {
        return(xx[, z])
      }
    })
    
    output[yr_ind, b_ind] <- do.call(cbind, match_columns)
    
    return(output)
  }, bins = bins, Year = Year$YearC)


dim_CAL <- c(dim(len[[1]]), length(len))
CAL <- len %>% unlist() %>% array(dim_CAL)
MSEtool::plot_composition(Year$YearC, CAL[,,1], CAL_bins = bins)
MSEtool::plot_composition(Year$YearC, CAL[,,2], CAL_bins = bins)

############# Save data
data <- list(Chist = catch_matrix, Index = Index, I_sd = SD, CAL = CAL, 
            length_bin = bins + 0.5)

saveRDS(data, file = "OM/SKJ_data_07.27.2020.rds")
