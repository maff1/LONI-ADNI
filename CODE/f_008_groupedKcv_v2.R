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

source("./CODE/f_utils.r")
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

pheno <- "BPSYS*VISIT"
fixVar = c(pheno,
            colnames(selVars)[11:48])
refVar = "(1|RID) + (1|SITE)"
varMet = colnames(selVars)[49:ncol(selVars)]

sFormula <- as.formula(
  paste("VAL",
        paste(c(fixVar, refVar), collapse = "+"), sep = "~")
)

set.seed(2021)
fold <- fold_cv(selVars, k = 10)

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
             REML = FALSE,
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='nlminb')
             )
  )
  newpred=predict(newlmx,newdata=testing, allow.new.levels = TRUE)
  true=testing$VAL
  error=(true-newpred)
  rmse=sqrt(mean(error^2))
  mse=mean((newpred-true)^2)
  R2=1-(sum((true-newpred)^2)/sum((true-mean(true))^2))
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
# ---------------------------------------------------------------------------- #