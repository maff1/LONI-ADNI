################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F005 LINEAR MIXED MODELS:
#    + test multi-phenotypes
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(lme4)
library(afex)
library(optimx)
library(parallel)

source("./CODE/f_utils.r")

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")
metDat$PHC_EXF <- as.numeric(metDat$PHC_EXF)
metDat$PHC_MEM <- as.numeric(metDat$PHC_MEM)
metDat$PHC_LAN <- as.numeric(metDat$PHC_LAN)

# PARALLEL IMPLEMENTATION LMX MODELS ----------------------------------------- #
# phenotypes = list("VSBPSYS", "VSBPDIA", "DX", "MMSCORE", "PHC_MEM", "PHC_EXF", "PHC_LAN")
# names(phenotypes) = c("BPSYS", "BPDIA", "DX", "MMSE", "MEMORY", "EXECUTIVE", "LANGUAGE")
phenotypes = list("PHC_MEM", "PHC_EXF", "PHC_LAN")
names(phenotypes) = c("MEMORY", "EXECUTIVE", "LANGUAGE")

for(pheno in 1:length(phenotypes)) {
  
ncores = 8
fixVar = c(phenotypes[[pheno]],
           "GENDER","AGE","EDUCAT","BMI","RACCAT","ALCH","SMOK","VISIT_NUM",
           grep("^meds", colnames(metDat), value = T)[-15],
           colnames(metDat)[272:289])
refVar = "(1 + VISIT_NUM | RID) + (1 + VISIT_NUM | SITEID)"
varMet = colnames(metDat)[4:253]

lsFit <- mclapply(varMet,
                  mc.cores = ncores,
                  function(sMet) {
                    sFormula <- as.formula(
                      paste(sMet,
                            paste(c(fixVar, refVar), collapse = "+"), sep = "~")
                    )
                    fit <- lmer(formula = sFormula,
                                data = metDat,
                                REML = FALSE,
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

saveRDS(modCoef, paste0("./RESULTS/lxm_coef_", names(phenotypes[pheno]), ".rds"))
}
# ---------------------------------------------------------------------------- #