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
library(parallel)

source("./CODE/f_utils.r")

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")

# PARALLEL IMPLEMENTATION LMX MODELS ----------------------------------------- #
ncores = 8
fixVar = c("GENDER*DX","AGE","EDUCAT","BMI","RACCAT","ALCH","SMOK","VISIT_NUM")
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
                                control = lmerControl()
                    )
                    resCoef <- data.frame(summary(fit)$coefficients,
                                          AIC(fit), 
                                          BIC(fit), 
                                          rsq::rsq(fit)$model, 
                                          nrow(metDat),
                                          metabolite = paste0(sMet))
                    colnames(resCoef) <- c("beta", "se", "df", "t.value", "pvalue",
                                           "aic", "bic", "r2", "sample.size", "metabolite")
                    raef <- subset(as.data.frame(ranef(fit)), 
                             term == "VISIT_NUM" & grpvar == "RID")
                    
                    lsRes <- list(raef = raef, resCoef = resCoef)
                    
                    return(lsRes)
                  } )
names(lsFit) <- varMet
modCoef = data.frame(
  coefTerm = unlist(lapply(strsplit(rownames(do.call(rbind, lapply(lsFit, function(x) x$resCoef))), split = "\\."),
                           function(s) s[2])),
  do.call(rbind, lapply(lsFit, function(x) x$resCoef)),
  row.names = NULL
)
modRanef = data.frame(
  metabolite = gsub("\\..*", "", rownames(
    do.call(rbind, lapply(lsFit, function(x) x$raef)))),
  do.call(rbind, lapply(lsFit, function(x) x$raef)),
  row.names = NULL
)

saveRDS(modCoef, "./RESULTS/lxm_ng_longitudinal_coef.rds")
saveRDS(modRanef, "./RESULTS/lxm_ng_longitudinal_ranef.rds")
# --------------------------------------------------------------------------
