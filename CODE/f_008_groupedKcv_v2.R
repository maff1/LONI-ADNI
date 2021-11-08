################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F008 10-FOLD CROSS-testing(CV):
#    + BY GROUPED K-FOLDS
#    + MAE
#    + MSE
#    + RMSE
#    + R2
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(caret)
library(Metrics)
library(afex)
library(optimx)
library(dplyr)
library(parallel)

source("./CODE/f_utils.r")

args <- commandArgs(trailingOnly=TRUE)
pheno = args[1]
ncores = as.integer(args[2])
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
# ---------------------------------------------------------------------------- #

fixVar = colnames(selVars)[11:48]
refVar = "(1|RID) + (1|SITE)"
varMet = colnames(selVars)[49:ncol(selVars)]

ls10cv = mclapply(varMet, mc.cores = ncores, function(sMet){
  sFormula <- as.formula(
    paste(
      paste(sMet,
            paste(pheno, "VISIT", sep = "*"), sep = "~"),
      paste(c(fixVar, refVar), collapse = "+", sep = "+"),
      sep = "+"
    )
  )
  
  set.seed(2021)
  fold <- fold_cv(selVars, k = 10)
  
  #vMet <- deparse(substitute(sMet))
  
  # createDF <- function(empty) {dat <- data.frame(
  #   RID = selVars$RID,
  #   Fold = empty,
  #   value = selVars[, colnames(selVars) %in% sMet],
  #   holdoutpred = empty,
  #   MSE = empty,
  #   RMSE = empty,
  #   MAE = empty,
  #   R2 = empty,
  #   AIC = empty,
  #   BIC = empty
  #   )
  # }
  # dat <- createDF(empty = paste0(rep(0, nrow(selVars))))
  
  dat <- selVars %>%
    mutate(Fold=rep(0,nrow(selVars)),
           holdoutpred=rep(0,nrow(selVars)),
           MSE=rep(0,nrow(.)),
           RMSE=rep(0,nrow(.)),
           MAE=rep(0,nrow(.)),
           R2=rep(0,nrow(.)),
           AIC=rep(0,nrow(.)),
           BIC=rep(0,nrow(.)))
  
  for(i in 1:10){
    train=dat[fold$subsets[fold$which != i], ]
    testing=dat[fold$subsets[fold$which == i], ]
    newlmx=lmer(formula = sFormula,
                data = train,
                REML = TRUE,
                control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb')
                )
    )
    newpred=predict(newlmx, newdata=testing, allow.new.levels = TRUE)
    trueValue = unlist(testing[sMet])
    error=(trueValue-newpred)
    rmse=sqrt(mean(error^2))
    mse=mean((newpred-trueValue)^2)
    R2=1-(sum((trueValue-newpred)^2)/sum((trueValue-mean(trueValue))^2))
    mae=mean(abs(error))
    dat[fold$subsets[fold$which == i], ]$holdoutpred <- newpred
    dat[fold$subsets[fold$which == i], ]$RMSE=rmse
    dat[fold$subsets[fold$which == i], ]$MSE=mse
    dat[fold$subsets[fold$which == i], ]$MAE=mae
    dat[fold$subsets[fold$which == i], ]$R2=R2
    dat[fold$subsets[fold$which == i], ]$AIC=AIC(newlmx)
    dat[fold$subsets[fold$which == i], ]$BIC=BIC(newlmx)
    dat[fold$subsets[fold$which == i], ]$Fold=i
  }
  return(dat)
})
names(ls10cv) <- varMet
cv10F = data.frame(do.call(rbind, ls10cv), 
                   metabolites = rep(varMet, each = length(ls10cv$TOTAL_C$RID)))
saveRDS(cv10F, paste0("./RESULTS/", pheno, "_cv10_lmx.rds"))
# ---------------------------------------------------------------------------- #