################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F006 COHORT DESCRIPTIVE STATISTICS:
#    + BY VISIT TIME | DIAGNOSIS
# ---------------------------------------------------------------------------- #

library(arsenal)

tblTitle <- "ADNI - Nightingale (NMR) longitudinal"
selCols <- as.formula(
  paste("VISCODE2", paste(colnames(metDat)[c(260:262,265:269,272:289,291,293:296,298:311,313)],
                 collapse = "+"), sep = "~"
  )
)

tbl <- tableby(formula = selCols, 
               data = metDat) 
tbl1 = summary(tbl, title = tblTitle)
write2html(tbl1, 
file = "/t1-data/project/psychdementia/mfernand/PROJECTS/LONI-ADNI/RESULTS/tbl_cohort.html", 
title = tblTitle)
# ---------------------------------------------------------------------------- #