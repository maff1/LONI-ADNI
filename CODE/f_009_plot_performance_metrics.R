################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F009 DISTRIBUTION PERFORMANCE METRICS:
#    + BY PHENOTYPE
# ---------------------------------------------------------------------------- #

library(ggplot2)
library(data.table)
library(tidyverse)
library(dendsort)
library(circlize)
library(ComplexHeatmap)
library(gridBase)

lsPheno <- list.files(path = "./RESULTS", pattern = "_cv10", full.names = TRUE)
names(lsPheno) <- gsub("_cv10_lmx.rds", "", basename(lsPheno))
temp = readRDS(lsPheno[1])[, c(301:307)]
# ---------------------------------------------------------------------------- #

dat1 <- temp %>%
  group_by(metabolites) %>%
  dplyr::summarize(AIC = mean(AIC, na.rm=TRUE),
                   R2 = mean(R2, na.rm=TRUE), 
                   RMSE = mean(RMSE, na.rm=TRUE)
  )
datAIC <- data.frame(AIC = dat1[, 2], 
                     row.names = tolower(
                       gsub("_", ".", dat1$metabolites)
                     )
)

datR2 <- data.frame(R2 = dat1[, 3], 
                    row.names = tolower(
                      gsub("_", ".", dat1$metabolites)
                    )
)

datRMSE <- data.frame(RMSE = dat1[, 4], 
                      row.names = tolower(
                        gsub("_", ".", dat1$metabolites)
                      )
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

pdf(file = "./RESULTS/circosMetrics_BPDIA.pdf", width = 9, height = 8)
plot.new()
circle_size = unit(1, "snpc")
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
# ---------------------------------------------------------------------------- #
circos.par(gap.after = c(2, 2, 2, 15, 2))
circos.heatmap(datAIC, col = col_funAIC, 
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               rownames.side = "outside",
               track.height = 0.1,
               split = split,
               show.sector.labels = TRUE
)
circos.heatmap(datR2, col = col_funR2,
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               track.height = 0.1,
               split = split
)
circos.heatmap(datRMSE, col = col_funRMSE,
               bg.border = "black", 
               bg.lwd = 1, 
               bg.lty = 1,
               cell.border = "black", 
               cell.lty = 1, 
               cell.lwd = 0.25,
               track.height = 0.1,
               split = split,
               dend.side = "inside",
               dend.callback = function(dend, m, si) {
                 dendsort(dend)
               }
)
circos.track(track.index = 4, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 4) {
    cn = c("RMSE", "r2", "AIC")
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)
text(0,0,"BPDIA")
#circos.text(0,0, labels = "BPDIA")
# ---------------------------------------------------------------------------- #
upViewport()

h = dev.size()[2]
lgd_list = packLegend(lgd_aic, lgd_r2, lgd_rmse, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")
circos.clear()
dev.off()
# ---------------------------------------------------------------------------- #