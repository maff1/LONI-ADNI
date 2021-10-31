################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F002 COVARIATES RDS FROM CSV:
#    + save rds files with column specification
# ---------------------------------------------------------------------------- #

rm(list = ls())
library(data.table)
library(dplyr)
library(lubridate)
library(zoo)
library(stringr)
library(missRanger)

source("./CODE/f_utils.r")

# DATADIC.csv # ADNI1  ADNI2 ADNIGO
dict <- read.csv("./DATA/DATADIC.csv")
dict <- dict[dict$Phase %in% c("ADNI1","ADNI2","ADNIGO"),]

# RID FROM NG
ng.rid <- unique(
  data.table::fread("./DATA/ADNINIGHTINGALELONG_05_24_21.csv", 
                    sep = ",", 
                    select = c("RID"), 
                    data.table = FALSE
  )$RID
)
print(length(ng.rid))
# 1524

# list all csv and count RID from NG
lsfiles = list.files("./DATA", pattern = ".csv", recursive = TRUE, full.names = TRUE)
names(lsfiles) <- toupper(gsub(".csv|_2020-11_20", "", basename(lsfiles)))
f_read_count <- function(files){
  tryCatch({
    dat <- fread(files, sep = ",", 
                 select = c("RID"), data.table = FALSE, nThread = 4)
    out <- length(unique(dat[dat$RID %in% ng.rid,]))
  },
  error = function(e) NULL)
}

`%nin%` <- Negate(`%in%`)
lsfiles <- lsfiles[lsfiles %nin% grep("DICT", lsfiles, ignore.case = T, value = T)]

counts = lapply(lsfiles, f_read_count)
filtered.counts = Filter(function(x) all(c(1524) %in% x), counts)
# 26
tbls <- lsfiles[names(lsfiles) %in% unique(names(filtered.counts))]
# ---------------------------------------------------------------------------- #

# phenotypes(outcomes) and covariates
f_extract_data <- function(X) {
  # outcomes
  cols.adsp <- c("RID","EXAMDATE","VISCODE2","DX","PHC_MEM","PHC_EXF","PHC_LAN")
  # NA == -1
  cols.mmse <- c("RID","USERDATE","VISCODE2","SITEID","MMSCORE")
  
  # covariates
  cols.medhist <- c("RID","USERDATE","VISCODE2","SITEID",
              grep("^MH",colnames(read.csv("./DATA/MEDHIST.csv")), value = TRUE))
  cols.neuroexa <- c("RID","USERDATE","VISCODE2","SITEID","NXVISUAL","NXAUDITO","NXTREMOR",
                     "NXCONSCI","NXNERVE","NXMOTOR","NXFINGER","NXHEEL","NXSENSOR",
                     "NXTENDON","NXPLANTA","NXGAIT","NXABNORM")
  cols.physical <- c("RID","USERDATE","VISCODE2","SITEID","PXGENAPP","PXHEADEY","PXNECK",
                     "PXCHEST","PXHEART","PXABDOM","PXEXTREM","PXPERIPH","PXSKIN",
                     "PXMUSCUL","PXBACK","PXABNORM")
  cols.demog <- c("RID","USERDATE","VISCODE2","SITEID","PTGENDER","PTDOBMM","PTDOBYY","PTMARRY",
                  "PTEDUCAT","PTETHCAT","PTRACCAT")
  cols.vitals <- c("RID","USERDATE","VISCODE2","SITEID","VSWEIGHT","VSWTUNIT","VSHEIGHT",
                   "VSHTUNIT","VSBPSYS","VSBPDIA","VSPULSE","VSRESP","VSTEMP")
  cols.meds <- c("RID", "USERDATE", "VISCODE2", "CMMED", "PROPERNAME", "final.class.code", "final.class.name")
  
  adniMerge <- data.table::fread("./DATA/ADNIMERGE.csv", sep = ",", 
                                 na.strings = c("", "NA", "NaN"))
  cols.adnimerge <- colnames(adniMerge)[colnames(adniMerge) %nin% 
                                    c(grep("_bl", colnames(adniMerge), value = T),
                                      "Month", "M", "update_stamp")]
  
  lsCols <- list(ADSP_PHC_COGN=cols.adsp,
                 MMSE=cols.mmse,
                 MEDHIST=cols.medhist,
                 NEUROEXM=cols.neuroexa,
                 PHYSICAL=cols.physical,
                 PTDEMOG=cols.demog,
                 VITALS=cols.vitals,
                 RECCMED_ATC=cols.meds,
                 ADNIMERGE=cols.adnimerge)
  tokeep = X[names(X) %in% names(lsCols)]
  keys <- unique(c(names(tokeep), names(lsCols)))
  lsdat = setNames(mapply(c, tokeep[keys], lsCols[keys]), keys)
  
  if(!dir.exists("./CLEAN_DATA")) {dir.create("./CLEAN_DATA")}
  for(i in 1:length(lsdat)) {
    dat <- data.table::fread(lsdat[[i]][1], sep = ",", na.strings = c("", "NA", "NaN", -4), 
                             select = as.character(lsdat[[i]][2:length(lsdat[[i]])]), 
                             data.table = FALSE, nThread = 12)
    saveRDS(dat, paste0("./CLEAN_DATA/", names(lsdat[i]), ".rds"))
  }
}
f_extract_data(X=lsfiles)
# ---------------------------------------------------------------------------- #