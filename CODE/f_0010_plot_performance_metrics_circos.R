################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F009 DISTRIBUTION PERFORMANCE METRICS:
#    + BY PHENOTYPE
# ---------------------------------------------------------------------------- #

rm(list = ls())
library(ggplot2)
library(data.table)
library(tidyverse)
library(dendsort)
library(circlize)
library(ComplexHeatmap)
library(gridBase)

lsPheno <- list.files(path = "./RESULTS", pattern = "_cv10", full.names = TRUE)
names(lsPheno) <- gsub("_cv10_lmx.rds", "", basename(lsPheno))
lsTemp = lapply(lsPheno, function(x) {
  df <- readRDS(x)[, c(301:307)] %>%
    group_by(metabolites) %>%
    dplyr::summarize(AIC = mean(AIC, na.rm=TRUE),
                     R2 = mean(R2, na.rm=TRUE), 
                     RMSE = mean(RMSE, na.rm=TRUE)
    )
  return(as.data.frame(df))
}
)

lsAIC = lapply(lsTemp, function(x) data.frame(x["metabolites"], x["AIC"]))
lsR2 = lapply(lsTemp, function(x) data.frame(x["metabolites"], x["R2"]))
lsRMSE = lapply(lsTemp, function(x) data.frame(x["metabolites"], x["RMSE"]))
# ---------------------------------------------------------------------------- #

f_rename_merge <- function(X) {
  for(i in 1:length(lsAIC)) {
    colnames(lsAIC[[i]])[2] <- paste0(colnames(lsAIC[[i]])[2], "_", names(lsAIC)[i])
  }
}

mergedDF = Reduce(function(...) merge(..., by = "metabolites", sort = FALSE), lsAIC)
colnames(mergedDF) <- gsub("AIC_", "", colnames(mergedDF))
rownames(mergedDF) <- mergedDF$metabolites; mergedDF <- mergedDF[, -1]
# ---------------------------------------------------------------------------- #

metAnn <- read.csv("./DATA/NG-annotations.csv")
met <- metAnn[match(toupper(unique(temp$metabolites)), toupper(metAnn$CSV_column_name)),]
met$sector <- factor(met$Group, 
                     levels = unique(met$Group), 
                     abbreviate(unique(met$Group), 3, named = FALSE)
)
# ---------------------------------------------------------------------------- #

split = factor(letters[1:5], levels = letters[1:5])
col_funAIC = colorRamp2(breaks = c(0, 7500, 15000), 
                        colors = c("white", "yellow", "red"), transparency = 0.25)
col_funR2 = colorRamp2(breaks = c(0, 0.5, 1), 
                       colors = c("white", "yellow", "red"), transparency = 0.25)
col_funRMSE = colorRamp2(breaks = c(0, 0.5, 1), 
                         colors = c("white", "yellow", "red"), transparency = 0.25)

lgd_aic = Legend(title = "AIC", col_fun = col_funAIC)
lgd_r2 = Legend(title = "r2", col_fun = col_funR2)
lgd_rmse = Legend(title = "RMSE", col_fun = col_funRMSE)
# ---------------------------------------------------------------------------- #

pdf(file = "./RESULTS/circosMetrics_AIC_split.pdf", width = 9, height = 8)
circos.par(gap.after = c(2, 15))
circos.heatmap(mergedDF, col = col_funAIC, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               rownames.side = "outside",
               track.height = 0.1,
               split = split2,
               show.sector.labels = TRUE,
               dend.side = "inside",
               dend.callback = function(dend, m, si) {
                 dendsort(dend)
               }
)

circos.track(track.index = 2, 
             panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 2) {
    cn = rev(colnames(mergedDF))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 0.25, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)
circos.clear()
lgd = Legend(title = "AIC", col_fun = col_funAIC)
grid.draw(lgd)
dev.off()
# ---------------------------------------------------------------------------- #

for(i in 1:length(lsR2)) {
  colnames(lsR2[[i]])[2] <- paste0(colnames(lsR2[[i]])[2], "_", names(lsR2)[i])
}
mergedR2 = Reduce(function(...) merge(..., by = "metabolites", sort = FALSE), lsR2)
colnames(mergedR2) <- gsub("R2_", "", colnames(mergedR2))
rownames(mergedR2) <- mergedR2$metabolites; mergedR2 <- mergedR2[, -1]

split2 = as.factor(ifelse(
  rownames(mergedR2) %in% grep("_", rownames(mergedR2), value = TRUE), "b", "a"
  )
)
# ---------------------------------------------------------------------------- #

pdf(file = "./RESULTS/circosMetrics_R2_split.pdf", width = 9, height = 8)
circos.par(gap.after = c(2, 2, 2, 15, 2))
circos.heatmap(mergedR2, col = col_funR2, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               rownames.side = "outside",
               track.height = 0.1,
               split = split2,
               show.sector.labels = TRUE,
               dend.side = "inside",
               dend.callback = function(dend, m, si) {
                 dendsort(dend)
               }
)

circos.track(track.index = 2, 
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 4) {
                 cn = rev(colnames(mergedR2))
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                             1:n - 0.5, cn,
                             cex = 0.25, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)
circos.clear()
lgd_r2 = Legend(title = "AIC", col_fun = col_funR2)
grid.draw(lgd_r2)
dev.off()
# ---------------------------------------------------------------------------- #

for(i in 1:length(lsR2)) {
  colnames(lsR2[[i]])[2] <- paste0(colnames(lsR2[[i]])[2], "_", names(lsR2)[i])
}
mergedR2 = Reduce(function(...) merge(..., by = "metabolites", sort = FALSE), lsR2)
colnames(mergedR2) <- gsub("R2_", "", colnames(mergedR2))
rownames(mergedR2) <- mergedR2$metabolites; mergedR2 <- mergedR2[, -1]

# ---------------------------------------------------------------------------- #

pdf(file = "./RESULTS/circosMetrics_R2_split.pdf", width = 9, height = 8)
circos.par(gap.after = c(2, 15))
circos.heatmap(mergedR2, col = col_funR2, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               rownames.side = "outside",
               track.height = 0.1,
               split = split2,
               show.sector.labels = TRUE,
               dend.side = "inside",
               dend.callback = function(dend, m, si) {
                 dendsort(dend)
               }
)

circos.track(track.index = 2, 
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 2) {
                 cn = rev(colnames(mergedR2))
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                             1:n - 0.5, cn,
                             cex = 0.25, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)
circos.clear()
lgd_r2 = Legend(title = "r2", col_fun = col_funR2)
grid.draw(lgd_r2)
dev.off()
# ---------------------------------------------------------------------------- #

for(i in 1:length(lsRMSE)) {
  colnames(lsRMSE[[i]])[2] <- paste0(colnames(lsRMSE[[i]])[2], "_", names(lsRMSE)[i])
}
mergedRMSE = Reduce(function(...) merge(..., by = "metabolites", sort = FALSE), lsRMSE)
colnames(mergedRMSE) <- gsub("RMSE_", "", colnames(mergedRMSE))
rownames(mergedRMSE) <- mergedRMSE$metabolites; mergedRMSE <- mergedRMSE[, -1]

pdf(file = "./RESULTS/circosMetrics_RMSE_split.pdf", width = 9, height = 8)
circos.par(gap.after = c(2, 15))
circos.heatmap(mergedRMSE, col = col_funRMSE, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               rownames.side = "outside",
               track.height = 0.1,
               split = split2,
               show.sector.labels = TRUE,
               dend.side = "inside",
               dend.callback = function(dend, m, si) {
                 dendsort(dend)
               }
)

circos.track(track.index = 2, 
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 2) {
                 cn = rev(colnames(mergedRMSE))
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                             1:n - 0.5, cn,
                             cex = 0.25, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)
circos.clear()
lgd_rmse = Legend(title = "RMSE", col_fun = col_funRMSE)
grid.draw(lgd_rmse)
dev.off()