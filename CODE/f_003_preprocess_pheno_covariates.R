################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F003 READ RDS AND CLEAN PHENOTYPES AND COVARIATES:
#    + save rds files with column specification
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(data.table)
library(dplyr)
library(lubridate)
library(missRanger)
library(stringr)

source("./CODE/f_utils.r")

f_lag_date <- function(met, demog, vitals, medhist, mmse, adsp, meds) {
  
  ng <- readRDS(met) %>% group_by(RID) %>% arrange(EXAMDATE) %>%
    mutate(SEQUENCE = 1:n(),
           MET.DIFF.DAYS = ifelse(is.na(ymd(EXAMDATE) - ymd(lag(EXAMDATE))), 0,
                                  ymd(EXAMDATE) - ymd(lag(EXAMDATE)))
    ) %>%
    mutate(MET.CUMUL.DAYS = cumsum(MET.DIFF.DAYS))
  # -------------------------------------------------------------------------- #
  
  # VITALS
  vitals <- readRDS(vitals)[, c(1,2,5:10)]
  vitals <- vitals[vitals$RID %in% unique(ng$RID),]
  vitals <- cbind.data.frame(vitals[, 1:2], apply(vitals[, 3:ncol(vitals)], 2, 
                                                  function(x) as.numeric(x)))
  vitals <- cbind.data.frame(vitals[, 1:2], apply(vitals[, 3:ncol(vitals)], 2, 
                                                  function(x) ifelse(x < 0, NA, x)))
  vitals$USERDATE <- as.Date(vitals$USERDATE, "%Y-%m-%d")
  # imputation carry forward 1
  vitals <- vitals %>%
    group_by(RID) %>%
    arrange(USERDATE) %>%
    do(na.locf2(.)) %>%
    do(na.locf3(.))
  # 1 pound (.lb) = 0.453592 kg
  # 2 --> kg; 1 --> lb
  # 1 inch = 2.54 cm
  # 1 --> inch; 2 --> cm
  
  vitals <- vitals %>%
    mutate(VSWEIGHT = ifelse(VSWTUNIT == 1, VSWEIGHT * 0.453592, VSWEIGHT),
           VSHEIGHT = ifelse(VSHTUNIT == 1, VSHEIGHT * 2.54, VSHEIGHT))
  vitals <- vitals %>%
    select(-c(VSWTUNIT,VSHTUNIT)) %>%
    mutate(VSWEIGHT = ifelse(VSWEIGHT > 140, NA,ifelse(VSWEIGHT < 40, NA, VSWEIGHT)),
           VSHEIGHT = ifelse(VSHEIGHT > 200, NA, ifelse(VSHEIGHT < 140, NA, VSHEIGHT)))
  # imputation carry forward 2
  vitals <- vitals %>%
    group_by(RID) %>%
    arrange(USERDATE) %>%
    do(na.locf2(.)) %>%
    do(na.locf3(.))
  
  # imputation with RF
  vitals <-  missRanger::missRanger(vitals, seed = 1234, num.trees = 20)
  print(anyNA(vitals))
  
  # BMI = Weight (kg) / Height (m)^2 
  vitals <- vitals %>% mutate(BMI = VSWEIGHT/(VSHEIGHT/100)^2) %>%
    select(-VSWEIGHT, -VSHEIGHT)
  # -------------------------------------------------------------------------- #
  
  # DEMOG
  demog <- readRDS(demog)
  demog <- demog[demog$RID %in% unique(ng$RID),]
  demog <- demog[!duplicated(demog$RID),]
  demog$USERDATE <- as.Date(demog$USERDATE, "%Y-%m-%d")
  demog$PTGENDER <- ifelse(demog$PTGENDER == 2, 0, 1)
  
  demog$BIRTH <- paste0(demog$PTDOBYY, "-", 
                        ifelse(demog$PTDOBMM < 10, 
                               paste("0", demog$PTDOBMM, sep=""), 
                               demog$PTDOBMM), "-", "01")
  demog$BIRTH <- as.Date(demog$BIRTH, "%Y-%m-%d")
  demog <- demog %>% select(-PTDOBMM, -PTDOBYY)
  print(anyNA(demog))
  # -------------------------------------------------------------------------- #
  
  # MEDHIST
  medhist = readRDS(medhist)
  medhist <- medhist[medhist$RID %in% unique(ng$RID),]
  na.medhist = colSums(is.na(medhist)/nrow(medhist)*100)
  medhist = medhist[colnames(medhist) %in% 
                      names(na.medhist)[na.medhist == 0][-c(4,5,20:22,24:25,27:29,32:33)]]
  medhist$USERDATE <- as.Date(medhist$USERDATE, "%Y-%m-%d")
  print(anyNA(medhist))
  # -------------------------------------------------------------------------- #
  
  # MMSE
  mmse = readRDS(mmse)[, -c(4)]
  mmse = mmse[mmse$RID %in% unique(ng$RID), ]
  mmse$USERDATE <- as.Date(mmse$USERDATE, "%Y-%m-%d")
  mmse$MMSCORE <- ifelse(mmse$MMSCORE == -1, NA, mmse$MMSCORE)
  print(anyNA(mmse))
  print(table(is.na(mmse$MMSCORE)))
  
  # imputation carry forward 1
  mmse <- mmse %>%
    group_by(RID) %>%
    arrange(USERDATE) %>%
    do(na.locf2(.)) %>%
    do(na.locf3(.))
  
  # -------------------------------------------------------------------------- #
  
  # ADSP
  adsp = readRDS(adsp)
  adsp = adsp[adsp$RID %in% unique(ng$RID), ]
  adsp = as.data.frame(apply(adsp, 2, function(x) ifelse(x == "NA", NA, x)))
  
  # imputation
  adsp <- adsp %>%
    group_by(RID) %>%
    arrange(EXAMDATE) %>%
    do(na.locf2(.)) %>%
    do(na.locf3(.))
  
  print(anyNA(adsp))
  # -------------------------------------------------------------------------- #
  
  # MEDICATIONS
  meds = readRDS(meds)[, c(1,2,6)]
  meds = meds[!is.na(meds$final.class.code), ]
  colnames(meds) <- c("RID", "USERDATE", "ATC")
  meds$USERDATE <- gsub("\\/", "-", meds$USERDATE)
  meds$ATC <- gsub("\\,.*", "", meds$ATC)
  meds$USERDATE <- as.Date(meds$USERDATE, "%m-%d-%Y")
  meds$ATC.LVL1 <- str_extract(meds$ATC, "^[A-Z]")
  
  setDT(meds, key = c("RID", "USERDATE"))
  meds = dcast(meds, RID + USERDATE ~ ATC.LVL1, 
               fun = length, value.var = "ATC.LVL1")
  setDF(meds)
  meds = cbind.data.frame(meds[, c(1:2)], ifelse(meds[, c(3:ncol(meds))] >= 1, 1, 0))
  print(anyNA(meds))
  # -------------------------------------------------------------------------- #
  
  # MERGE BY RID AND DATE
  ls.vars = list(ng, vitals, demog, medhist, mmse, adsp, meds)
  names(ls.vars) = c("met", "vitals", "demog", "medhist", "mmse", "adsp", "meds")
  ls.vars = lapply(ls.vars, function(x) as.data.frame(x))
  
  ls.vars = lapply(ls.vars, function(x) 
    data.frame(RID = as.numeric(trimws(x[, c("RID")], 
                                       which = "both")), x[, 2:length(x)]))
  
  res = lapply(ls.vars, function(x) data.frame(x[, 1:length(x)], JOINDATE = if(
    isTRUE(colnames(x)[colnames(x) %in% c("EXAMDATE")] == "EXAMDATE")) {
    JOINDATE = as.Date(x[, c("EXAMDATE")], "%Y-%m-%d")} else {
      JOINDATE = as.Date(x[, c("USERDATE")], "%Y-%m-%d")}, EVENTDATE = if(
        isTRUE(colnames(x)[colnames(x) %in% c("EXAMDATE")] == "EXAMDATE")) {
        EVENTDATE = as.Date(x[, c("EXAMDATE")], "%Y-%m-%d")} else {
          EVENTDATE = as.Date(x[, c("USERDATE")], "%Y-%m-%d")}
  )
  )
  
  dt = lapply(res, function(x) setDT(x, key = c("RID","JOINDATE")))
  dt.met = dt$met
  dt.join = lapply(dt, function(x) x[dt.met, on = .(RID, JOINDATE), 
                                     roll = "nearest"][order(RID)])
  print(paste0("Before.......................................................")) 
  print(lapply(dt.join, function(x) table(duplicated(x[, c("RID", "JOINDATE")]))))
  dt.join = lapply(dt.join, function(x) unique(x, by = c("RID", "JOINDATE")))
  print(paste0("After........................................................")) 
  print(lapply(dt.join, function(x) table(duplicated(x[, c("RID", "JOINDATE")]))))
  print(paste0("RID all equal?........................................................"))
  print(lapply(dt.join[2:length(dt.join)], function(x) all.equal(x$RID, dt.join$met$RID)))
  #dt.join = Reduce(function(...) merge(..., all = TRUE), dt)
  dt.join = lapply(dt.join, function(x) x[, 
                                          DIFF.DAYS := as.numeric(difftime(EVENTDATE, JOINDATE, units = "days"))])
  dt.join = lapply(dt.join, function(x) x[, ("i.EVENTDATE"):=NULL])
  drop.cols <- grep("i.", colnames(dt.join$met))
  dt.join$met <- dt.join$met[, !..drop.cols]
  
  mincol <- unlist(Reduce(intersect, lapply(dt.join, colnames)))
  mincol <- mincol[-which(mincol == "DIFF.DAYS")]
  out = lapply(dt.join[2:length(dt.join)], function(x) x[, !..mincol])
  dt.out = cbind(dt.join$met, do.call(cbind, out))
  
  # AGE @blood draw-metabolites date
  dt.out$MET.AGE <- (dt.out$EXAMDATE - dt.out$demog.BIRTH)/365.25
  
  setDF(dt.out)
  dt.out <- dt.out[, colnames(dt.out) %nin% 
                     grep("i.viscode", colnames(dt.out), ignore.case = T, value = T)]
  saveRDS(dt.out, "./CLEAN_DATA/ng_data_matrix.rds")
  return(dt.out)
}

ng = f_lag_date(met = "./CLEAN_DATA/imputed_ng_longitudinal.rds", 
                demog = "./CLEAN_DATA/PTDEMOG.rds",
                vitals = "./CLEAN_DATA/VITALS.rds",
                medhist = "./CLEAN_DATA/MEDHIST.rds",
                mmse = "./CLEAN_DATA/MMSE.rds",
                adsp = "./CLEAN_DATA/ADSP_PHC_COGN.rds",
                meds = "./CLEAN_DATA/RECCMED_ATC.rds")
# ---------------------------------------------------------------------------- #