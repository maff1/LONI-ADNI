################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F014 PLOT MARGIN EFFECTS
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(lme4)
library(afex)
library(optimx)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggeffects)
source("./CODE/f_utils.r")

metDat <- readRDS("./CLEAN_DATA/ng_data_matrix.rds")
selVars <- data.frame(RID = as.factor(metDat$RID),
                      VISIT = as.factor(metDat$VISCODE2),
                      VISIT_NUM = as.numeric(metDat$VISCODE2),
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
# ---------------------------------------------------------------------------- #

fixVar = colnames(selVars)[c(12:49)]
refVar = "(1|RID) + (1|SITE)"
varMet = colnames(selVars)[50:ncol(selVars)]

sFormula <- as.formula(
  paste(
    paste("CREATININE",
          paste("DX", "VISIT_NUM", sep = "*"), sep = "~"),
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

# ---------------------------------------------------------------------------- #

pdf("./RESULTS/creatinine_marginal_plot.pdf")
  ggplot(selVars, aes(x=VISIT_NUM, y=CREATININE, group=DX)) +
  geom_point(aes(color=DX),size=1,alpha=0.5) +
    scale_x_continuous(breaks=seq(1,13,1), labels = levels(selVars$VISIT),
                       name="VISIT") +
  stat_summary(fun="median", geom="line", aes(color=DX), size=1) +
  facet_wrap(~DX, ncol=1, scales="free") +
  theme_bw()
dev.off()
# ---------------------------------------------------------------------------- #

preddf<-selVars%>%mutate(PredCREATININE=predict(fit,.))

pdf("./RESULTS/predPHE_marginal_plot.pdf")
ggplot(preddf, aes(x=VISIT_NUM, y=PredCREATININE, group=DX))+
  geom_point(aes(color=DX), size=1, alpha=0.5) +
  scale_x_continuous(breaks=seq(1,13,1), labels = levels(selVars$VISIT),
                     name="VISIT") +
  stat_summary(fun="median", geom="line", aes(color=DX), size=1) +
  stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max,aes(color=DX),size=1)+
  stat_summary(aes(fill=DX),fun.ymin=min,fun.ymax=max,fun.y="median",shape=21,color="black",size=5,geom ="point")+
  facet_wrap(~DX, ncol=1, scales="free") +
  theme_bw()
dev.off()
# ---------------------------------------------------------------------------- #

creat <- ggpredict(fit, terms = c("VISIT_NUM", "DX"), type = "random")
pdf("./RESULTS/CREAT_marginal_plot.pdf")
p1 = plot(creat, ci = FALSE)
p1 + scale_x_continuous(breaks=seq(1,13,1), labels = levels(selVars$VISIT),
                        name="VISIT")
dev.off()