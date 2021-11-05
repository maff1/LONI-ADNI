################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F007 10-FOLD CROSS-VALIDATION(CV):
#    + BY GROUPED K-FOLDS
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(caret)
library(Metrics)
library(afex)
library(optimx)

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")

# 10% out
subjects = unique(metDat[["RID"]])
n_folds = 10
folds = createFolds(y = subjects, k = n_folds, list = TRUE, returnTrain = TRUE)

mae_id <- c()
mse_id <- c()
rmse_id <- c()

for(i in 1:n_folds) {
  
  train_ids = subjects[folds[[i]]]
  test_ids = subjects[!(subjects %in% train_ids)]
  
  train_dat = droplevels(subset(metDat, RID %in% train_ids))
  test_dat = droplevels(subset(metDat, RID %in% test_ids))
  
  pheno <- "MMSCORE"
  fixVar = c(pheno,
             "GENDER","AGE","EDUCAT","BMI","RACCAT","ALCH","SMOK","VISIT_NUM",
             grep("^meds", colnames(metDat), value = T)[-15],
             colnames(metDat)[272:289])
  refVar = "(1 + VISIT_NUM | RID) + (1 + VISIT_NUM | SITEID)"
  varMet = colnames(metDat)[4:253]
  
  sFormula <- as.formula(
    paste("VAL",
          paste(c(fixVar, refVar), collapse = "+"), sep = "~")
  )
  fit <- lmer(formula = sFormula,
              data = metDat,
              REML = FALSE,
              control = lmerControl(
                optimizer ='optimx', optCtrl=list(method='nlminb')
              )
  )
  
  test_subjects = unique(test_dat[["RID"]])
  
  for(subject in test_subjects){
    pred_subset = subset(test_dat, RID == subject)
    
    y_pred = predict(fit, newdata = pred_subset, allow.new.levels = TRUE)
    y_true = pred_subset[["MMSCORE"]]
    
    mae_id <- c(mae_id, Metrics::mae(y_true, y_pred))
    mse_id <- c(mse_id, Metrics::mse(y_true, y_pred))
    rmse_id <- c(mse_id, Metrics::rmse(y_true, y_pred))
  }
}