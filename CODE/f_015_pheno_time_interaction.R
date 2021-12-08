################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F015 PHENO*TIME
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(afex)
library(optimx)
library(dplyr)
library(parallel)

setwd("/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI")
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

fixVar = colnames(selVars)[11:48]
refVar = "(1|RID) + (1|SITE)"
varMet = colnames(selVars)[49:ncol(selVars)]
# ---------------------------------------------------------------------------- #

phenotypes = list("BPSYS", "BPDIA", "DX", "MMSE", "MEM", "EXF", "LAN")
names(phenotypes) = c("BPSYS", "BPDIA", "DX", "MMSE", "MEM", "EXF", "LAN")


for(pheno in 1:length(phenotypes)) {

lsFit <- mclapply(varMet,
                  mc.cores = 6,
                  function(sMet) {
                    sFormula <- as.formula(
                      paste(
                        paste(sMet,
                              paste(phenotypes[[pheno]], "VISIT", sep = "*"), sep = "~"),
                        paste(c(fixVar, refVar), collapse = "+", sep = "+"),
                        sep = "+"
                      )
                    )
                    fit <- lmer(formula = sFormula,
                                data = selVars,
                                REML = TRUE,
                                control = lmerControl(
                                  optimizer ='optimx', optCtrl=list(method='nlminb')
                                )
                    )
                    dfCoef <- as.data.frame(summary(fit)$coefficients)
                    tbldfCoef <- data.frame(coef.term = rownames(summary(fit)$coefficients),
                                            beta = dfCoef$Estimate,
                                            se = dfCoef$`Std. Error`,
                                            ci95lo = dfCoef$Estimate - dfCoef$`Std. Error`*1.96,
                                            ci95up = dfCoef$Estimate + dfCoef$`Std. Error`*1.96,
                                            dfreedom = dfCoef$df,
                                            t.value = dfCoef$`t value`,
                                            p.value = dfCoef$`Pr(>|t|)`,
                                            aic = AIC(fit), 
                                            bic = BIC(fit), 
                                            r2 = rsq::rsq(fit)$model, 
                                            sample.size = nrow(metDat),
                                            metabolite = paste0(sMet),
                                            phenotype = phenotypes[[pheno]]
                    )
                    
                    return(tbldfCoef)
                  } )
names(lsFit) <- varMet
modCoef = data.frame(do.call(rbind, lsFit), row.names = NULL)

if(!dir.exists("./RESULTS/TIME_INT")) {dir.create("./RESULTS/TIME_INT")}

saveRDS(modCoef, paste0("./RESULTS/TIME_INT/lxm_coef_time_interaction_", names(phenotypes[pheno]), ".rds"))
}
# ---------------------------------------------------------------------------- #