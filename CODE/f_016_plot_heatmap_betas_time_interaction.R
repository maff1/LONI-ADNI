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

filesList = list.files("./RESULTS/TIME_INT", pattern = ".rds", full.names = TRUE)
names(filesList) <- gsub("lxm_coef_time_interaction_|.rds", "", basename(filesList))

f_get_coef <- function(p) {
  dat <- readRDS(paste0("./RESULTS/TIME_INT/lxm_coef_time_interaction_", p, ".rds"))
  dat <- subset(dat, 
                coef.term == paste(p, "VISITm12", sep = ":")|
                  coef.term == paste(p, "VISITm24", sep = ":")|
                    coef.term == paste(p, "VISITm36", sep = ":")
  )[, c(1,2,8,13)]
  
  return(dat)
}

lsPval <- lapply(names(filesList), f_get_coef)
resCoef <- do.call(rbind, lsPval)

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

pdf("./RESULTS/heatmap_lxm_time_interaction.pdf", height = 30, width = 15)
Heatmap(mBetas, 
        name = "betas", cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%s", mPvalue[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 3),
        column_title = "NMR LONGITUDINAL|raw pvalues|<0.001 ***, <0.01**, <0.05*, <0.1.")
dev.off()
# ---------------------------------------------------------------------------- #

filMetBetas <- mBetas[rownames(mBetas) %nin% grep("_", rownames(mBetas), value = TRUE),]
filMetPvalue <- mPvalue[rownames(mPvalue) %nin% grep("_", rownames(mPvalue), value = TRUE),]

timeLabels = gsub("^[A-Z]{1,5}\\.VISIT", "", colnames(filMetBetas))
timeLabels_col = colorRamp2(breaks = unique(timeLabels), colors = c("yellow", "pink", "brown"))


col_fun = colorRamp2(breaks = unique(timeLabels), colors = c("yellow", "pink", "brown"))
ha = HeatmapAnnotation(
  time = timeLabels,
  col = list(
             time = c("m12" = "yellow", "m24" = "pink", "m36" = "brown")
  ),
  gp = gpar(col = "black")
)

# ---------------------------------------------------------------------------- #

pdf("./RESULTS/heatmap_lxm_time_interaction_selected_MET_LABELS.pdf")
Heatmap(filMetBetas, 
        name = "betas", cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%s", filMetPvalue[i, j]), x, y, gp = gpar(fontsize = 12))
        },
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_title = "NMR LONGITUDINAL|raw pvalues|<0.001 ***, <0.01**, <0.05*, <0.1.",
        top_annotation = ha)
dev.off()
# ---------------------------------------------------------------------------- #