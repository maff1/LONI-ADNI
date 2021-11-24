################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F013 PLOT HEATMAP BETAS
# ---------------------------------------------------------------------------- #

rm(list = ls())

library(data.table)
library(ComplexHeatmap)
library(ggplot2)
source("./CODE/f_utils.r")

# MAKE MATRIX ---------------------------------------------------------------- #

resF = readRDS("./RESULTS/bootBetas.rds")
setDT(resF)
mres <- dcast(resF, metabolite~coef, value.var = "beta.median", fun.aggregate = mean)

lsCoef = list.files("./RESULTS", pattern = "lxm_coef_", full.names = TRUE)
names(lsCoef) <- gsub("lxm_coef_|.rds", "", basename(lsCoef))


f_get_pval <- function(p) {
  dat <- readRDS(p)
  dat <- subset(dat, 
                coef.term == "VSBPDIA"|
                  coef.term == "VSBPSYS"|
                  coef.term == "DXALZ"|
                  coef.term == "DXMCI"|
                  coef.term == "DXALZ"|
                  coef.term == "PHC_LAN"|
                  coef.term == "PHC_EXF"|
                  coef.term == "PHC_MEM"|
                  coef.term == "MMSCORE")[, c(1,2,8,13)]
  return(dat)
}

lsPval <- lapply(lsCoef, f_get_pval)
resCoef <- do.call(rbind, lsPval)
resCoef$coef.term <- factor(resCoef$coef.term, 
                            levels = unique(resCoef$coef.term),
                            labels = c("BPDIA","BPSYS","DXMCI","DXALZ",
                                       "EXF","LAN","MEM","MMSE")
                            )
setDT(resCoef)
mBetas <- dcast(resCoef, metabolite~coef.term, value.var = "beta", fun.aggregate = mean)
setDF(mBetas)
mBetas <- data.frame(mBetas, row.names = mBetas$metabolite); mBetas <- mBetas[, -1]

mPvalue <- dcast(resCoef, metabolite~coef.term, value.var = "p.value", fun.aggregate = mean)
setDF(mPvalue)
mPvalue <- data.frame(mPvalue, row.names = mPvalue$metabolite); mPvalue <- mPvalue[, -1]
mPvalue <- ifelse(mPvalue < 0.001, "***",
                  ifelse(mPvalue < 0.01, "**",
                         ifelse(mPvalue < 0.05, "*",
                                ifelse(mPvalue < 0.1, ".", "")
                                )
                         )
                  )
# ---------------------------------------------------------------------------- #

pdf("./RESULTS/heatmap_lxm.pdf", height = 30, width = 15)
Heatmap(mBetas, 
        name = "betas", cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%s", mPvalue[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 3),)
dev.off()
# ---------------------------------------------------------------------------- #

filMetBetas <- mBetas[rownames(mBetas) %nin% grep("_", rownames(mBetas), value = TRUE),]
filMetPvalue <- mPvalue[rownames(mPvalue) %nin% grep("_", rownames(mPvalue), value = TRUE),]

pdf("./RESULTS/heatmap_lxm_selected_MET.pdf")
Heatmap(filMetBetas, 
        name = "betas", cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%s", filMetPvalue[i, j]), x, y, gp = gpar(fontsize = 12))
        },
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 10),)
dev.off()
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #