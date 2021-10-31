################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F001 PREPROCESS NG:
#    + MERGE DUPLICATES
#    + IF OUTLIER ADD NA
#    + IMPUTE NA's WITH DOUBLE PASS
#    + LOG2(X+1) TRANSFORM
#    + SAVE CLEAN DATA
# periods of time aligned using VISCODE2
# ---------------------------------------------------------------------------- #

rm(list=ls())
library(data.table)
library(missRanger)
source("./CODE/f_utils.r")

# longitudinal NG-NMR metabolite data
f_clean_ng_longitudinal <- function(metFile) {
  ng <- data.table::fread(metFile, sep = ",", na.strings = c("", "NA", "NaN"), 
                         data.table = FALSE, nThread = 4, select = c(1,3,30:279))
  ng <- as.data.frame(apply(ng, 2, function(x) ifelse(x == "TAG", NA, x)))
  ng <- cbind.data.frame(ng[, 1:2], apply(ng[,3:ncol(ng)], 2, 
                                          function(x) as.numeric(x, na.rm = TRUE)))
  setDT(ng)
  ng[, RID_VISIT := trimws(paste0(RID, "_", VISCODE2), which = "both")][, c("RID", "VISCODE2"):=NULL]
  setkey(ng, "RID_VISIT")
  
  # remove variables or individuals-visit with 100% NA ----------------------- #
  colNA <- data.frame(metName = colnames(ng[, 1:250]), 
                      percNA = round(colSums(is.na(ng[, 1:250]))/nrow(ng)*100, digits = 2),
                      row.names = NULL
  )
  rowNA <- data.frame(rid_visit = ng$RID_VISIT, 
                      percNA = round(rowSums(is.na(ng[, 1:250]))/ncol(ng[, 1:250])*100, digits = 2),
                      row.names = NULL
  )
  `%nin%` <- Negate(`%in%`)
  ng <- ng[ng$RID_VISIT %nin% subset(rowNA, percNA == 100)$rid_visit, ]
  print(paste0("removed...................", 
               length(subset(rowNA, percNA == 100)$rid_visit), " data points")
  )
  anyDuplicated(ng$RID_VISIT) # 40
  ng <- ng[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                                       by="RID_VISIT", .SDcols = 1:(ncol(ng)-1)]
  # -------------------------------------------------------------------------- #
  setDF(ng)
  ng <- cbind(data.frame(
     RID = as.numeric(trimws(unlist(lapply(
       strsplit(ng$RID_VISIT, split = "_"), function(x) x[1])
     ), which = "both")),
    VISCODE2 = trimws(unlist(lapply(
      strsplit(ng$RID_VISIT, split = "_"), function(x) x[2])
    ), which = "both"),
    row.names = NULL), ng[, 2:ncol(ng)]
  )
  
  # UNIVARIATE OUTLIERS ------------------------------------------------------ #
  mad <- apply(ng[, 3:ncol(ng)], 2, isnt_out_mad)
  ng <- cbind.data.frame(ng[,1:2], ifelse(mad == FALSE, NA, TRUE) * ng[,3:ncol(ng)])
  print(
    all(apply(ng[, 3:ncol(ng)], 2, is.numeric))
  )
  
  # double pass imputation per visit
  # imputation can take a while...depends on tree number --------------------- #

  lsVisits <- lapply(names(table(ng$VISCODE2)), function(x) subset(ng, VISCODE2 == x))
  names(lsVisits) <- names(table(ng$VISCODE2))
  
  lsNgImp_pass_1 <- lapply(lsVisits, 
                    function(X) cbind.data.frame(X[, 1:2],
                                                   missRanger::missRanger(
                                                     data = X[, 3:ncol(X)],
                                                     formula = .~.,
                                                     seed = 1234,
                                                     num.trees = 20,
                                                     maxiter = 10
        )
    )
)
  ###
  lsNgImp_pass_2 <- lapply(lsNgImp_pass_1, 
                           function(X) cbind.data.frame(X[, 1:2],
                                                        log2(
                                                          missRanger::missRanger(
                                                            data = X[, 3:ncol(X)],
                                                            formula = .~.,
                                                            seed = 1234,
                                                            num.trees = 20,
                                                            maxiter = 10L) + 1
                                                        )
                           )
  )
  ###
  ngImp <- do.call(rbind, lsNgImp_pass_2)
  ngImp$RID_VISIT <- paste(ngImp$RID, ngImp$VISCODE2, sep = "_")
  ###
  # add EXAMDATE
  ngExamDate <- data.table::fread(metFile, 
                                  sep = ",", na.strings = c("", "NA", "NaN"), 
                                  data.table = FALSE, nThread = 4, select = c(1,3,4))
  ngExamDate$EXAMDATE <- as.Date(ngExamDate$EXAMDATE, "%Y-%m-%d")
  ngExamDate$RID_VISIT <- paste(ngExamDate$RID, ngExamDate$VISCODE2, sep = "_")
  ngExamDate <- ngExamDate[, c("RID_VISIT", "EXAMDATE")]
  ngExamDate <- ngExamDate[!duplicated(ngExamDate$RID_VISIT),]
  
  out <- merge(ngImp, ngExamDate, by = "RID_VISIT", all = FALSE, all.x = TRUE)
  out <- data.frame(RID = out[, "RID"], 
                    EXAMDATE = out[, "EXAMDATE"], 
                    VISCODE2 = out[, "VISCODE2"],
                    out[, 4:253],
                    row.names = NULL)
  
  ###
  if(!dir.exists("./CLEAN_DATA")) {
    dir.create("./CLEAN_DATA")
  }
  saveRDS(out, "./CLEAN_DATA/imputed_ng_longitudinal.rds")
  return(out)
}
# ---------------------------------------------------------------------------- #
ng = f_clean_ng_longitudinal(metFile = "./DATA/ADNINIGHTINGALELONG_05_24_21.csv")
# ---------------------------------------------------------------------------- #