################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F011 BOOTSTRAP MODEL COEFFICIENTS
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(afex)
library(optimx)
library(dplyr)
library(parallel)
library(boot)

setwd("/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI")
source("./CODE/f_utils.r")

args <- commandArgs(trailingOnly=TRUE)
pheno = args[1]
sMet = args[2]
ncores = as.integer(args[3])
nsim = as.integer(args[4])
# ---------------------------------------------------------------------------- #

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")
selVars <- data.frame(RID = as.factor(metDat$RID),
                      VISIT = as.factor(metDat$VISCODE2),
                      SITE = as.factor(metDat$SITEID),
                      BPSYS = as.numeric(metDat$VSBPSYS),
                      BPDIA = as.numeric(metDat$VSBPDIA),
                      DX = as.factor(metDat$DX),
                      MEM = as.numeric(metDat$PHC_MEM),
                      LAN = as.numeric(metDat$PHC_LAN),
                      EXF = as.numeric(metDat$PHC_EXF),
                      MMSE = as.numeric(metDat$MMSCORE),
                      GENDER = as.factor(metDat$GENDER),
                      AGE = as.numeric(metDat$AGE),
                      RACE = as.factor(metDat$RACCAT),
                      ETHNICITY = as.factor(metDat$ETHCAT),
                      EDUCATION = as.numeric(metDat$EDUCAT),
                      BMI = as.numeric(metDat$BMI),
                      metDat[, c(298:311)],
                      metDat[, c(272:289)],
                      metDat[, c(4:253)]
)

scCols <- c("BPSYS", "BPDIA", "MMSE", "AGE", "BMI", "EDUCATION",
            colnames(selVars)[49:ncol(selVars)])
selVars[scCols] <- lapply(selVars[scCols], zscaleFunc)

fixVar = colnames(selVars)[11:48]
refVar = "(1|RID) + (1|SITE)"
varMet = colnames(selVars)[49:ncol(selVars)]

sFormula <- as.formula(
    paste(
      paste(sMet,
            paste(pheno, "VISIT", sep = "*"), sep = "~"),
      paste(c(fixVar, refVar), collapse = "+", sep = "+"),
      sep = "+"
  )
)
# ---------------------------------------------------------------------------- #

mod1 = lmer(formula = sFormula,
            data = selVars,
            REML = TRUE,
            control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb')
            )
)
b_par <- bootMer(x = mod1, 
                 FUN = fixef, 
                 nsim = nsim, 
                 seed = 2021,
        type = "parametric",
        parallel = "multicore",
        ncpus = 8)

# CI based on quantile approach ---------------------------------------------- #

tblOut <- data.frame(
  coef = colnames(b_par$t),
  pheno = paste0(pheno),
  beta.median = apply(b_par$t, 2, median),
  ci.lower = apply(b_par$t, 2, function(x) as.numeric(quantile(x, probs=0.025, na.rm=TRUE))),
  ci.upper = apply(b_par$t, 2, function(x) as.numeric(quantile(x, probs=0.975, na.rm=TRUE))),
  row.names = NULL
)
# ---------------------------------------------------------------------------- #
saveRDS(tblOut, paste0("./RESULTS/BOOT/", sMet, "_", pheno, "_boot_", nsim, ".rds")
        )
