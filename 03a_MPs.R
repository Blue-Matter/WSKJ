
library(MSEtool)
library(DLMextra)

# MP wrapper that puts the BB fleet in Data@Ind
MP_wrapper <- function(MP, delay = 2, ...) {
  MP <- substitute(MP)
  MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data), reps = quote(reps), list(...)))
  force(delay)
  
  MP_body <- bquote({
    Data@Ind <- Data@AddInd[, 1, ]
    Data@CV_Ind <- Data@CV_AddInd[, 1, ]
    if(delay > 0) {
      delay_ind <- seq(ncol(Data@Cat) - delay + 1, ncol(Data@Cat), 1)
      Data@Ind <- Data@Ind[, -delay_ind]
      Data@CV_Ind <- Data@CV_Ind[, -delay_ind]
      Data@Cat <- Data@Cat[, -delay_ind, drop = FALSE]
      Data@CV_Cat <- Data@CV_Cat[, -delay_ind, drop = FALSE]
      Data@Year <- Data@Year[-delay_ind]
    }
    Rec <- .(MP_call)
    return(Rec)
  })
  
  MP_out <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  class(MP_out) <- "MP"
  return(MP_out)
}

# Constant catch MPs (mean 2015-2018 catch is about 20,000 t with 75%-25% baitboat to handline)
CC_30kt <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(3e4, reps)
  return(Rec)
}
class(CC_30kt) <- "MP"

CC_20kt <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(2e4, reps)
  return(Rec)
}
class(CC_20kt) <- "MP"

CC_15kt <- function(x, Data, reps) {
  Rec <- new("Rec")
  Rec@TAC <- rep(1.5e4, reps)
  return(Rec)
}
class(CC_15kt) <- "MP"

Iratio_ <- Iratio_3u <- MP_wrapper(Iratio)
Islope <- Islope_3u <- MP_wrapper(Islope1, xx = 0, yrsmth = 3, lambda = 0.6)
GBslope <- GBslope_3u <- MP_wrapper(GB_slope, yrsmth = 3)

# spict
spict_4010 <- make_MP(spict, HCR40_10, fix_alpha = FALSE, fix_beta = FALSE, diagnostic = "min")
spict_4010_ <- spict_4010_3u <- MP_wrapper(spict_4010)
