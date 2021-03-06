################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F004 LINEAR MIXED MODELS:
#    + mod1 with all vars
#    + mod2 prunned mod1
#    + mod3 prunned mod2
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(lme4)
library(afex)
library(optimx)
library(parallel)

source("./CODE/f_utils.r")

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")

# PARALLEL IMPLEMENTATION LMX MODELS ----------------------------------------- #
ncores = 8
fixVar = c("GENDER*DX","AGE","EDUCAT","BMI","RACCAT","ALCH","SMOK","VISIT_NUM",
           grep("^meds", colnames(metDat), value = T)[-15],
           colnames(metDat)[272:289])
refVar = "(1 + VISIT_NUM | RID) + (1 + VISIT_NUM | SITEID)"
varMet <- colnames(metDat)[4:253]

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
                                            metabolite = paste0(sMet)
                    )
                    raef <- subset(as.data.frame(ranef(fit)), 
                             term == "VISIT_NUM" & grpvar == "RID")
                    
                    lsRes <- list(raef = raef, resCoef = tbldfCoef)
                    
                    return(lsRes)
                  } )
names(lsFit) <- varMet
modCoef = data.frame(do.call(rbind, lapply(lsFit, function(x) x$resCoef)),
  row.names = NULL
)
modRanef = data.frame(
  metabolite = gsub("\\..*", "", rownames(
    do.call(rbind, lapply(lsFit, function(x) x$raef)))),
  do.call(rbind, lapply(lsFit, function(x) x$raef)),
  row.names = NULL
)

saveRDS(modCoef, "./RESULTS/lxm_ng_longitudinal_gender_coef.rds")
saveRDS(modRanef, "./RESULTS/lxm_ng_longitudinal_gender_ranef.rds")
# ---------------------------------------------------------------------------- #